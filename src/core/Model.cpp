#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#define MODEL_PI 3.141592653589793238462643

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Spectra/SymEigsSolver.h>

#include "Elements/Elements.h"
#include "Components.h"
#include "LoadComponent.h"

#include "Model.h"
#include "SparseMatrixUtils.h"

#include "analysis/FEDynamic.h"

FEModel::FEModel()
    : RigidLinkData(std::make_shared<RigidLinks>())
{
}

namespace {

// 振動固有値問題 K φ = ω^2 M φ を、質量あり自由度上の標準固有値問題
// A y = μ y, A = L^T P^T K^{-1} P L, μ = 1/ω^2 として解くための Spectra 用作用素。
// P は質量あり自由度の選択、L は質量ブロック M_mm = L L^T のコレスキー因子
// （元の質量行列は対角だが、剛体リンク縮約後のマスタ自由度には連成が入るため
// 対角とは限らない）。シューア補元（静的縮約）を陽に作らないため K の疎性が
// 保たれ、K の分解は1回で済む。
class VibrationInverseOp {
public:
    using Scalar = double;

    VibrationInverseOp(ISparseSolver& solver, const std::vector<int>& massed_indices,
        const Eigen::SparseMatrix<double>& mass_cholL, int full_size)
        : solver_(solver), massed_indices_(massed_indices),
        mass_cholL_(mass_cholL), full_size_(full_size) {}

    Eigen::Index rows() const { return (Eigen::Index)massed_indices_.size(); }
    Eigen::Index cols() const { return (Eigen::Index)massed_indices_.size(); }

    // y_out = L^T P^T K^{-1} P L x_in
    void perform_op(const double* x_in, double* y_out) const
    {
        int n = (int)massed_indices_.size();
        Eigen::VectorXd v = mass_cholL_ * Eigen::Map<const Eigen::VectorXd>(x_in, n);

        Eigen::VectorXd full_rhs = Eigen::VectorXd::Zero(full_size_);
        for (int i = 0; i < n; i++)
            full_rhs(massed_indices_[i]) = v(i);

        Eigen::VectorXd sol = solver_.solve(full_rhs);

        Eigen::VectorXd g(n);
        for (int i = 0; i < n; i++)
            g(i) = sol(massed_indices_[i]);

        Eigen::Map<Eigen::VectorXd>(y_out, n) = mass_cholL_.transpose() * g;
    }

private:
    ISparseSolver& solver_;
    const std::vector<int>& massed_indices_;
    const Eigen::SparseMatrix<double>& mass_cholL_;
    int full_size_;
};

} // namespace

int FEModel::SolveVibration(const int nev, std::vector<double>& eigen_values,
    std::vector<std::vector<Displacement>>& mode_vectors)
{
    int computed_num = nev;

    // インデックスの取得（RigidLinkを考慮）
    std::vector<int> slave_indices = RigidLinkData->SlaveDOFIndices();
    std::vector<int> free_indices = FreeIndices(true);  // rigid_link=true
    std::vector<int> fixed_indices = FixIndices();

    // 変換行列の取得
    Eigen::SparseMatrix<double> linkTransMat =
        RigidLinkData->TransformationMatrix().sparseView(1e-10);
    int master_dof_num = linkTransMat.cols();

    Eigen::SparseMatrix<double> ka, ma;

    if (master_dof_num > 0) {
        // RigidLinkがある場合: 3x3ブロックに分割して縮小
        Eigen::SparseMatrix<double> k11, k12, k13, k22, k23, k33;
        SparseMatrixUtils::splitMatrix3x3(AssembleStiffnessMatrix(), slave_indices, free_indices,
            k11, k12, k13, k22, k23, k33);

        Eigen::SparseMatrix<double> m11, m12, m13, m22, m23, m33;
        SparseMatrixUtils::splitMatrix3x3(AssembleMassMatrix(), slave_indices, free_indices,
            m11, m12, m13, m22, m23, m33);

        // 剛性行列の縮小
        Eigen::SparseMatrix<double> kaa, kab;
        kaa = (linkTransMat.transpose() * k11.selfadjointView<Eigen::Upper>() * linkTransMat)
              .triangularView<Eigen::Upper>();
        kab = (linkTransMat.transpose() * k12);
        SparseMatrixUtils::mergeMatrixWithResize(kaa, kab, k22, ka);

        // 質量行列の縮小
        Eigen::SparseMatrix<double> maa, mab;
        maa = (linkTransMat.transpose() * m11.selfadjointView<Eigen::Upper>() * linkTransMat)
              .triangularView<Eigen::Upper>();
        mab = (linkTransMat.transpose() * m12);
        SparseMatrixUtils::mergeMatrixWithResize(maa, mab, m22, ma);
    }
    else {
        // RigidLinkがない場合: 従来通り2x2分割
        SparseMatrixUtils::splitMatrixWithResize(AssembleStiffnessMatrix(), fixed_indices, ka);
        SparseMatrixUtils::splitMatrixWithResize(AssembleMassMatrix(), fixed_indices, ma);
    }

    // 質量あり自由度の抽出
    // 質量ゼロ自由度は行列を縮約せず、VibrationInverseOp が暗黙に静的縮約と
    // 同じ固有値問題を解く（対角質量では対角ゼロ⇔行・列全体ゼロなので厳密）
    std::vector<int> massed_indices, massless_indices;
    Eigen::VectorXd mdiag = ma.diagonal();
    for (int i = 0; i < mdiag.size(); i++)
    {
        if (mdiag(i) >= 0.0000001)
            massed_indices.push_back(i);
        else
            massless_indices.push_back(i);
    }

    int mat_size = (int)massed_indices.size();
    if (computed_num > mat_size - 1)
        computed_num = mat_size - 1;
    if (computed_num < 1 || mat_size - 2 < computed_num)
        return -1;

    // 質量あり自由度ブロック M_mm のコレスキー分解 M_mm = L L^T
    // （剛体リンクのマスタ自由度は質量が連成するため対角とは限らない。
    //   並べ替えを行わない NaturalOrdering で置換の扱いを不要にする）
    Eigen::SparseMatrix<double> m_sh;
    SparseMatrixUtils::splitMatrixWithResize(ma, massless_indices, m_sh);
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper,
        Eigen::NaturalOrdering<int>> mass_llt(m_sh);
    if (mass_llt.info() != Eigen::Success)
        return -1;
    Eigen::SparseMatrix<double> mass_cholL = mass_llt.matrixL();

    // 剛性が全く付かない自由度（トラス節点の回転、板のドリリング等）で
    // K が特異になるのを防ぐ。PSD の組立行列では対角ゼロ⇔行・列全体ゼロ
    // （完全非連成）なので、質量ゼロの死自由度の対角に正値を置いても
    // 他自由度の解は変わらず、当該モード成分は 0 になる。
    {
        Eigen::VectorXd kdiag = ka.diagonal();
        std::vector<Eigen::Triplet<double>> reg;
        for (int i = 0; i < kdiag.size(); i++)
        {
            if (kdiag(i) > 0.0)
                continue;
            if (mdiag(i) >= 0.0000001)
                return -1; // 質量があるのに剛性ゼロの自由度は解けない
            reg.emplace_back(i, i, 1.0);
        }
        if (!reg.empty()) {
            Eigen::SparseMatrix<double> kreg(ka.rows(), ka.cols());
            kreg.setFromTriplets(reg.begin(), reg.end());
            ka += kreg;
        }
    }

    // K の分解は全体でこの1回のみ
    auto solver_vib = createSolver();
    if (!solver_vib->compute(ka))
        return -1;

    int ncv = 2 * computed_num + 1; // Recommended value
    if (ncv > mat_size) ncv = mat_size;

    // A = L^T P^T K^{-1} P L の最大固有値 μ = 1/ω^2 を求める（ω 昇順で得られる）
    VibrationInverseOp op(*solver_vib, massed_indices, mass_cholL, (int)ka.rows());
    Spectra::SymEigsSolver<VibrationInverseOp> geigs(op, computed_num, ncv);

    geigs.init();
    int nconv = geigs.compute();

    if (geigs.info() != Spectra::CompInfo::Successful)
        return -1;

    Eigen::VectorXd mu = geigs.eigenvalues();   // μ 降順 = ω 昇順
    Eigen::MatrixXd u1s = geigs.eigenvectors(); // 正規直交 (y^T y = 1)

    // 固有値を元の固有値問題に戻す（ω = 角振動数）
    for (int i = 0; i < nconv; i++)
        eigen_values.push_back(1.0 / sqrt(mu(i)));

    // 縮小空間（master + free）の固有ベクトルを復元: φ = K^{-1} P L y / μ
    // 質量ゼロ自由度の成分も静的縮約関係を満たす形で同時に得られる。
    // φ^T M φ = y^T y = 1 となり質量正規化が構成上厳密に成立する。
    Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(ka.rows(), nconv);
    for (int j = 0; j < nconv; j++)
    {
        Eigen::VectorXd v = mass_cholL * u1s.col(j);
        for (int i = 0; i < mat_size; i++)
            rhs(massed_indices[i], j) = v(i);
    }

    Eigen::MatrixXd reduced_vectors = solver_vib->solveMulti(rhs);
    for (int j = 0; j < nconv; j++)
        reduced_vectors.col(j) /= mu(j);

    // 全体DOFへの固有ベクトルを構築
    Eigen::MatrixXd eigs_vector = Eigen::MatrixXd::Zero(DOFNum(), nconv);

    if (master_dof_num > 0) {
        // RigidLinkがある場合: master DOFをslave DOFに展開
        for (int i = 0; i < nconv; i++) {
            Eigen::VectorXd part_vec = reduced_vectors.col(i);
            Eigen::VectorXd d_master = part_vec.head(master_dof_num);
            Eigen::VectorXd d_free = part_vec.tail(free_indices.size());
            Eigen::VectorXd d_slave = linkTransMat * d_master;

            for (size_t j = 0; j < slave_indices.size(); j++)
                eigs_vector(slave_indices[j], i) = d_slave(j);
            for (size_t j = 0; j < free_indices.size(); j++)
                eigs_vector(free_indices[j], i) = d_free(j);
        }
    }
    else {
        // RigidLinkがない場合: free_indicesに直接配置
        std::vector<int> all_free = FreeIndices(false);
        for (int j = 0; j < nconv; j++)
            for (size_t i = 0; i < all_free.size(); i++)
                eigs_vector(all_free[i], j) = reduced_vectors(i, j);
    }

    // 固有ベクトルを格納（φ^T M φ = 1 の質量正規化済み）
    for (int i = 0; i < nconv; i++)
    {
        std::vector<Displacement> v(NodeNum());
        for (size_t j = 0; j < NodeNum(); j++)
        {
            int p = j * 6;
            v[j] = Displacement(
                eigs_vector(p, i), eigs_vector(p + 1, i), eigs_vector(p + 2, i),
                eigs_vector(p + 3, i), eigs_vector(p + 4, i), eigs_vector(p + 5, i));
        }
        mode_vectors.push_back(v);
    }

    return nconv;
}

std::vector<int> FEModel::FreeIndices(bool rigid_link)
{
    std::vector<int> slave_indices = RigidLinkData->SlaveDOFIndices();
    std::unordered_set<int> slaveid_set(slave_indices.begin(), slave_indices.end());
    std::vector<int> indices;
    int idx = 0;
    for (Node n : Nodes)
        for (const bool fixed : n.Fix.isdof_fixed()) {
            
            bool is_slave = (slaveid_set.find(idx) != slaveid_set.end());
            if (rigid_link){
                if (!fixed && !is_slave)
                    indices.push_back(idx);
            }
            else if (!fixed)
                indices.push_back(idx);
            
            idx++;
        }
    return indices;
}

std::vector<int> FEModel::SlaveIndices()
{
    return RigidLinkData->SlaveDOFIndices();
}

/// <summary>
/// 固定された自由度について全体自由度におけるインデックスを格納
/// </summary>
/// <returns>長さが固定自由度数で全体自由度インデックスが格納されたint型vector</returns>
std::vector<int> FEModel::FixIndices()
{
    std::vector<int> indices;
    int idx = 0;
    for (Node n : Nodes) {
        auto fixed = n.Fix.isdof_fixed();
        for (bool f : fixed) {
            if (f) indices.push_back(idx);
            idx++;
        }
    }
    return indices;
}

void FEModel::add_element(BeamElement data)
{
    std::shared_ptr<BeamElement> ptr = std::make_shared<BeamElement>(data);
    Elements.push_back(ptr);
    
    Nodes[data.Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[data.Nodes[1]->id].Fix.UnlockAllRot();
}

void FEModel::add_element(ComplexBeamElement data)
{
    std::shared_ptr<ComplexBeamElement> ptr = std::make_shared<ComplexBeamElement>(data);
    Elements.push_back(ptr);

    Nodes[data.Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[data.Nodes[1]->id].Fix.UnlockAllRot();
}

void FEModel::add_element(TrussElement data)
{
    std::shared_ptr<TrussElement> ptr = std::make_shared<TrussElement>(data);
    Elements.push_back(ptr);
}

void FEModel::add_element(TensionTrussElement data)
{
    std::shared_ptr<TensionTrussElement> ptr = std::make_shared<TensionTrussElement>(data);
    Elements.push_back(ptr);
}

void FEModel::add_element(TriPlaneElement data)
{
    std::shared_ptr<TriPlaneElement> ptr = std::make_shared<TriPlaneElement>(data);
    Elements.push_back(ptr);
}

void FEModel::add_element(TriPlateElement data)
{
    std::shared_ptr<TriPlateElement> ptr = std::make_shared<TriPlateElement>(data);
    Elements.push_back(ptr);

    Nodes[ptr->Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[1]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[2]->id].Fix.UnlockAllRot();
}

void FEModel::add_element(QuadPlaneElement data)
{
    std::shared_ptr<QuadPlaneElement> ptr = std::make_shared<QuadPlaneElement>(data);
    Elements.push_back(ptr);
}

void FEModel::add_element(QuadPlateElement data)
{
    std::shared_ptr<QuadPlateElement> ptr = std::make_shared<QuadPlateElement>(data);
    Elements.push_back(ptr);

    Nodes[ptr->Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[1]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[2]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[3]->id].Fix.UnlockAllRot();
}

void FEModel::add_truss_element(int id, int n1_id, int n2_id, int sec_id, int mat_id)
{
    Node* n1 = &Nodes[n1_id];
    Node* n2 = &Nodes[n2_id];
    Section* sec = &Sections[sec_id];
    Material mat = Materials[mat_id];

    std::shared_ptr<TrussElement> ptr = std::make_shared<TrussElement>(id, n1, n2, sec, mat);
    Elements.push_back(ptr);
}

void FEModel::add_beam_element(int id, int n1_id, int n2_id, int sec_id, int mat_id, double beta)
{
    Node* n1 = &Nodes[n1_id];
    Node* n2 = &Nodes[n2_id];
    Section* sec = &Sections[sec_id];
    Material mat = Materials[mat_id];

    std::shared_ptr<BeamElement> ptr = std::make_shared<BeamElement>(id, n1, n2, sec, mat, beta);
    Elements.push_back(ptr);

    Nodes[ptr->Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[1]->id].Fix.UnlockAllRot();
}

void FEModel::add_tri_plate_element(int id, int n1_id, int n2_id, int n3_id, double thickness, int mat_id)
{
    Node* n1 = &Nodes[n1_id];
    Node* n2 = &Nodes[n2_id];
    Node* n3 = &Nodes[n3_id];
    Material mat = Materials[mat_id];

    std::shared_ptr<TriPlateElement> ptr = std::make_shared<TriPlateElement>(id, n1, n2, n3, thickness, mat);
    Elements.push_back(ptr);

    Nodes[ptr->Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[1]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[2]->id].Fix.UnlockAllRot();

}

void FEModel::add_quad_plate_element(int id, int n1_id, int n2_id, int n3_id, int n4_id, double thickness, int mat_id)
{
    Node* n1 = &Nodes[n1_id];
    Node* n2 = &Nodes[n2_id];
    Node* n3 = &Nodes[n3_id];
    Node* n4 = &Nodes[n4_id];
    Material mat = Materials[mat_id];

    std::shared_ptr<QuadPlateElement> ptr = 
        std::make_shared<QuadPlateElement>(id, n1, n2, n3, n4, thickness, mat);
    Elements.push_back(ptr);

    Nodes[ptr->Nodes[0]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[1]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[2]->id].Fix.UnlockAllRot();
    Nodes[ptr->Nodes[3]->id].Fix.UnlockAllRot();
}

BarElementBase* FEModel::GetBarElement(int id)
{
    BarElementBase* be = dynamic_cast<BarElementBase*>(Elements[id].get());
    return be;
}

BeamElement* FEModel::GetBeamElement(int id)
{
    BeamElement* be = dynamic_cast<BeamElement*>(Elements[id].get());
    return be;
}

TrussElement* FEModel::GetTrussElement(int id)
{
    TrussElement* be = dynamic_cast<TrussElement*>(Elements[id].get());
    return be;
}

QuadPlateElement* FEModel::GetQuadPlateElement(int id)
{
    QuadPlateElement* be = dynamic_cast<QuadPlateElement*>(Elements[id].get());
    return be;
}

TriPlateElement* FEModel::GetTriPlateElement(int id)
{
    TriPlateElement* be = dynamic_cast<TriPlateElement*>(Elements[id].get());
    return be;
}

std::vector<NodeLoadData> FEModel::InnertialForceToNodeLoads(const InertialForce inertial_force)
{
    std::vector<NodeLoadData> node_loads;
    //if (std::shared_ptr<InertialForce> inertial = std::dynamic_pointer_cast<InertialForce>(load)) {
    for (auto& e : Elements)
    {
        std::vector<NodeLoadData> elem_loads =
            e->InertialForceToNodeLoadData(
                Eigen::Vector3d(
                    inertial_force.accels.x, 
                    inertial_force.accels.y, 
                    inertial_force.accels.z));
        node_loads.insert(node_loads.end(), elem_loads.begin(), elem_loads.end());
    }
    return node_loads;
}

double FEModel::SumNodeMass()
{
    Vector mass;
    double sum = 0;
    for (Node n : Nodes)
        sum += n.MassData.SumMass();

    return sum;
}

Eigen::SparseMatrix<double> FEModel::AssembleStiffnessMatrix(bool applyTensionOnly)
{
    int mat_size = Nodes.size() * 6;

    // 非零要素数を推定
    size_t estimated_nnz = Elements.size() * 100;
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(estimated_nnz);

    // 各要素からTripletを収集
    for (const std::shared_ptr<ElementBase>& eh : Elements) {
        if (applyTensionOnly) {
            // 状態依存要素は現在状態に応じた接線剛性を組立てる
            if (auto sde = std::dynamic_pointer_cast<IStateDependentElement>(eh)) {
                sde->GetTangentStiffnessTriplets(tripletList);
                continue;
            }
        }
        eh->GetStiffnessTriplets(tripletList);
    }

    // Tripletから疎行列を一括構築
    Eigen::SparseMatrix<double> mat(mat_size, mat_size);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
}

Eigen::SparseMatrix<double> FEModel::AssembleMassMatrix()
{
    std::vector<Eigen::Triplet<double>> tripletList;
    for (int i = 0; i < NodeNum(); ++i) {
        tripletList.push_back(Eigen::Triplet<double>(i * 6, i * 6, 0.0));
        tripletList.push_back(Eigen::Triplet<double>(i * 6 + 1, i * 6 + 1, 0.0));
        tripletList.push_back(Eigen::Triplet<double>(i * 6 + 2, i * 6 + 2, 0.0));
    }
    Eigen::SparseMatrix<double> mass_mat(NodeNum() * 6, NodeNum() * 6);
    mass_mat.setFromTriplets(tripletList.begin(), tripletList.end());

    // Construct MassMatrix
    //for each(std::shared_ptr<ElementBase> eh in Elements)
    //    eh->AssembleMassMatrix(mass_mat);

    for (const Node& n : Nodes)
    {
        int i = n.id * 6;
        mass_mat.coeffRef(i, i) += n.MassData.SumMass();
        mass_mat.coeffRef(i + 1, i + 1) += n.MassData.SumMass();
        mass_mat.coeffRef(i + 2, i + 2) += n.MassData.SumMass();
    }

    mass_mat *= (1.0 / GraityAccel);
    return mass_mat;
}

Eigen::SparseMatrix<double> FEModel::AssembleGeometricStiffnessMatrix(
    const std::vector<Displacement>& displacements)
{
    int mat_size = Nodes.size() * 6;

    size_t estimated_nnz = Elements.size() * 100;
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(estimated_nnz);

    // 各要素の節点変位を準備して Triplet を収集
    for (const std::shared_ptr<ElementBase>& eh : Elements) {
        std::vector<Displacement> disp_vec(eh->NodeNum());
        std::vector<Node*> nodes = eh->NodesList();
        for (size_t i = 0; i < eh->NodeNum(); i++)
            disp_vec[i] = displacements[nodes[i]->id];

        eh->GetGeometricStiffnessTriplets(disp_vec, tripletList);
    }

    Eigen::SparseMatrix<double> mat(mat_size, mat_size);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
}

void FEModel::ComputeElementNodeMass()
{
    for (size_t i = 0; i < Nodes.size(); i++)
        Nodes[i].MassData.ElementMass = 0.0;

    for (std::shared_ptr<ElementBase> eh : Elements) {
        Eigen::VectorXd masses = eh->NodeLumpedMass();
        for (size_t i = 0; i < eh->NodeNum(); i++) {
            int node_id = eh->NodesList()[i]->id;
            Nodes[node_id].MassData.ElementMass += masses(i);
        }
    }
}

void FEModel::SolveLinearStatic(std::vector<std::shared_ptr<LoadBase>>& loads, 
    std::vector<Displacement>& disp, std::vector<NodeLoad>& react)
{
    Eigen::VectorXd f(Nodes.size() * 6);
    f.setZero();
    for (auto& load : loads)
    {

        std::vector<NodeLoadData> node_loads;
        if (std::shared_ptr<InertialForce> inertial = std::dynamic_pointer_cast<InertialForce>(load)) {
            for (auto& e : Elements)
            {
                std::vector<NodeLoadData> elem_loads = 
                    e->InertialForceToNodeLoadData(Eigen::Vector3d(inertial->accels.x, inertial->accels.y, inertial->accels.z));
                node_loads.insert(node_loads.end(), elem_loads.begin(), elem_loads.end());
            }
        }
        else {
            node_loads = load->NodeLoads();
        }

        for (auto &nl : node_loads)
        {
            if (nl.id < 0) continue;
            int pos = nl.id * 6;
            f[pos] += nl.Px();
            f[pos + 1] += nl.Py();
            f[pos + 2] += nl.Pz();
            f[pos + 3] += nl.Mx();
            f[pos + 4] += nl.My();
            f[pos + 5] += nl.Mz();
        }
    }

    std::vector<int> slave_indices = RigidLinkData->SlaveDOFIndices();
    std::vector<int> free_indices = FreeIndices(true);
    std::vector<int> fixed_indices = FixIndices();

    Eigen::SparseMatrix<double> m11, m12, m13, m22, m23, m33;
    SparseMatrixUtils::splitMatrix3x3(AssembleStiffnessMatrix(), slave_indices, free_indices,
        m11, m12, m13, m22, m23, m33);

    Eigen::SparseMatrix<double> maa, mab, mac;
    Eigen::SparseMatrix<double> linkTransMat = RigidLinkData->TransformationMatrix().sparseView(1e-10);
    maa = (linkTransMat.transpose() * m11.selfadjointView<Eigen::Upper>() * linkTransMat).triangularView<Eigen::Upper>();
    mab = (linkTransMat.transpose() * m12);
    mac = (linkTransMat.transpose() * m13);

    Eigen::SparseMatrix<double> mii, mij, mjj;
    if (maa.rows() > 0){
        SparseMatrixUtils::mergeMatrixWithResize(maa, mab, m22, mii);
        mij = SparseMatrixUtils::vstack(mac, m23);
    }
    else{
        mii = m22;
        mij = m23;
    }
    mjj = m33;

    Eigen::VectorXd f_slave(slave_indices.size());
    Eigen::VectorXd f_free(free_indices.size());
    Eigen::VectorXd f_fix(fixed_indices.size());
    for (size_t i = 0; i < slave_indices.size(); i++)
        f_slave(i) = f(slave_indices[i]);
    for (size_t i = 0; i < free_indices.size(); i++)
        f_free(i) = f(free_indices[i]);
    for (size_t i = 0; i < fixed_indices.size(); i++)
        f_fix(i) = f(fixed_indices[i]);

    Eigen::VectorXd f_master = linkTransMat.transpose() * f_slave;
    Eigen::VectorXd f_input(f_master.size() + f_free.size());
    f_input << f_master, f_free;

    // Solve
    auto solver_static = createSolver();
    solver_static->compute(mii);
    Eigen::VectorXd d_result = solver_static->solve(f_input);
    Eigen::VectorXd r_fix = mij.transpose() * d_result - f_fix;

    // 反力データ整理
    Eigen::VectorXd r = Eigen::VectorXd::Zero(Nodes.size() * 6);
    for (size_t i = 0; i < fixed_indices.size(); i++)
        r(fixed_indices[i]) = r_fix(i);
    for (size_t i = 0; i < Nodes.size(); i++)
    {
        if (!Nodes[i].Fix.IsAnyFix()) continue;
        int pos = i * 6;
        react.push_back(NodeLoad(i, r[pos], r[pos + 1], r[pos + 2], r[pos + 3], r[pos + 4], r[pos + 5]));
    }

    // 変形データ整理
    Eigen::VectorXd d = Eigen::VectorXd::Zero(Nodes.size() * 6);

    // Eigen::VectorXd d_master = d_result.head(linkTransMat.cols());
    Eigen::VectorXd d_slave = linkTransMat * d_result.head(linkTransMat.cols());
    Eigen::VectorXd d_free = d_result.tail(free_indices.size());

    for (size_t i = 0; i < slave_indices.size(); i++)
        d(slave_indices[i]) = d_slave(i);

        for (size_t i = 0; i < free_indices.size(); i++)
        d(free_indices[i]) = d_free(i);

    for (size_t i = 0; i < Nodes.size(); i++)
    {
        int pos = i * 6;
        disp.push_back(Displacement(d[pos], d[pos + 1], d[pos + 2], d[pos + 3], d[pos + 4], d[pos + 5]));
    }
}

double IResponseSpectrum::DampingCorrectionFactor(const double t)
{
    if (DampInitializer == nullptr) return 1.0;

    double h = DampInitializer->DampRateAtPeriod(t);
    // 算定不能(t<=0のZPAや初期化前など、h<0)は補正なし。
    // これにより零周期応答(t=0)は減衰非依存となり、Fhが負に化けるのも防ぐ。
    if (h < 0.0) return 1.0;
    return 1.5 / (1.0 + 10.0 * h);
}

double IResponseSpectrum::acceleration_factored(double t)
{
    if (enable_damp_factor)
        return DampingCorrectionFactor(t) * Acceleration(t);
    else
        return Acceleration(t);
}

double IResponseSpectrum::velocity_factored(double t)
{
    if (enable_damp_factor)
        return DampingCorrectionFactor(t) * Velocity(t);
    else
        return Velocity(t);
}

double IResponseSpectrum::displacement_factored(double t)
{
    if (enable_damp_factor)
        return DampingCorrectionFactor(t) * Displacement(t);
    else
        return Displacement(t);
}

void FEModel::SolveLinearStaticIter(std::vector<std::shared_ptr<LoadBase>> &loads, std::vector<Displacement> &disp, std::vector<NodeLoad> &react)
{
    int max_iter = 100;

    Eigen::VectorXd force_vec(Nodes.size() * 6);
    Eigen::VectorXd residual_vec(Nodes.size() * 6);
    Eigen::VectorXd disp_vec = Eigen::VectorXd::Zero(Nodes.size() * 6);
    force_vec.setZero();
    for (auto &load : loads)
    {

        std::vector<NodeLoadData> node_loads;
        if (std::shared_ptr<InertialForce> inertial = std::dynamic_pointer_cast<InertialForce>(load))
        {
            for (auto &e : Elements)
            {
                std::vector<NodeLoadData> elem_loads =
                    e->InertialForceToNodeLoadData(Eigen::Vector3d(inertial->accels.x, inertial->accels.y, inertial->accels.z));
                node_loads.insert(node_loads.end(), elem_loads.begin(), elem_loads.end());
            }
        }
        else
        {
            node_loads = load->NodeLoads();
        }

        for (auto &nl : node_loads)
        {
            if (nl.id < 0)
                continue;
            int pos = nl.id * 6;
            force_vec[pos] += nl.Px();
            force_vec[pos + 1] += nl.Py();
            force_vec[pos + 2] += nl.Pz();
            force_vec[pos + 3] += nl.Mx();
            force_vec[pos + 4] += nl.My();
            force_vec[pos + 5] += nl.Mz();
        }
    }
    residual_vec = force_vec;

    std::vector<int> slave_indices = RigidLinkData->SlaveDOFIndices();
    std::vector<int> free_indices = FreeIndices(true);
    std::vector<int> fixed_indices = FixIndices();

    // 状態依存要素の状態を初期化 (規定剛性で機能する状態へ)
    for (auto &e : Elements)
        if (auto sde = std::dynamic_pointer_cast<IStateDependentElement>(e))
            sde->IsActive = true;
    
    Eigen::SparseMatrix<double> full_stiffmat = AssembleStiffnessMatrix(true);


    for (size_t iter = 0; iter < max_iter; iter++)
    {
        Eigen::SparseMatrix<double> m11, m12, m13, m22, m23, m33;
        SparseMatrixUtils::splitMatrix3x3(full_stiffmat, slave_indices, free_indices,
                                        m11, m12, m13, m22, m23, m33);

        Eigen::SparseMatrix<double> maa, mab, mac;
        Eigen::SparseMatrix<double> linkTransMat = RigidLinkData->TransformationMatrix().sparseView(1e-10);
        maa = (linkTransMat.transpose() * m11.selfadjointView<Eigen::Upper>() * linkTransMat).triangularView<Eigen::Upper>();
        mab = (linkTransMat.transpose() * m12);
        mac = (linkTransMat.transpose() * m13);

        Eigen::SparseMatrix<double> mii, mij, mjj;
        if (maa.rows() > 0)
        {
            SparseMatrixUtils::mergeMatrixWithResize(maa, mab, m22, mii);
            mij = SparseMatrixUtils::vstack(mac, m23);
        }
        else
        {
            mii = m22;
            mij = m23;
        }
        mjj = m33;

        Eigen::VectorXd f_slave(slave_indices.size());
        Eigen::VectorXd f_free(free_indices.size());
        Eigen::VectorXd f_fix(fixed_indices.size());
        for (size_t i = 0; i < slave_indices.size(); i++)
            f_slave(i) = residual_vec(slave_indices[i]);
        for (size_t i = 0; i < free_indices.size(); i++)
            f_free(i) = residual_vec(free_indices[i]);
        for (size_t i = 0; i < fixed_indices.size(); i++)
            f_fix(i) = residual_vec(fixed_indices[i]);

        Eigen::VectorXd f_master = linkTransMat.transpose() * f_slave;
        Eigen::VectorXd f_input(f_master.size() + f_free.size());
        f_input << f_master, f_free;

        // Solve
        auto solver_static = createSolver();
        solver_static->compute(mii);
        Eigen::VectorXd d_result = solver_static->solve(f_input);
        Eigen::VectorXd r_fix = mij.transpose() * d_result - f_fix;

        // 反力データ整理
        Eigen::VectorXd r = Eigen::VectorXd::Zero(Nodes.size() * 6);
        for (size_t i = 0; i < fixed_indices.size(); i++)
            r(fixed_indices[i]) = r_fix(i);
        react.clear();
        for (size_t i = 0; i < Nodes.size(); i++)
        {
            if (!Nodes[i].Fix.IsAnyFix())
                continue;
            int pos = i * 6;
            react.push_back(NodeLoad(i, r[pos], r[pos + 1], r[pos + 2], r[pos + 3], r[pos + 4], r[pos + 5]));
        }

        // 変形データ整理
        Eigen::VectorXd delta_disp_vec = Eigen::VectorXd::Zero(Nodes.size() * 6);

        // Eigen::VectorXd d_master = d_result.head(linkTransMat.cols());
        Eigen::VectorXd d_slave = linkTransMat * d_result.head(linkTransMat.cols());
        Eigen::VectorXd d_free = d_result.tail(free_indices.size());

        for (size_t i = 0; i < slave_indices.size(); i++)
            delta_disp_vec(slave_indices[i]) = d_slave(i);

        for (size_t i = 0; i < free_indices.size(); i++)
            delta_disp_vec(free_indices[i]) = d_free(i);
        
        disp_vec += delta_disp_vec;
        disp.clear();
        for (size_t i = 0; i < Nodes.size(); i++)
        {
            int pos = i * 6;
            disp.push_back(Displacement(disp_vec[pos], disp_vec[pos + 1], disp_vec[pos + 2], disp_vec[pos + 3], disp_vec[pos + 4], disp_vec[pos + 5]));
        }

        // 判定と更新
        // 状態依存要素の update を呼び、状態変化があれば剛性を再構築する
        bool any_change = false;
        for (const auto& elem : Elements) {
            if (auto sde = std::dynamic_pointer_cast<IStateDependentElement>(elem)) {
                if (sde->update(disp))
                    any_change = true;
            }
        }

        if (!any_change){
            break; // 収束判定: 状態変化なしなら終了
        }
        else{
            full_stiffmat = AssembleStiffnessMatrix(true); // 状態変化あり: 剛性行列を再構築
            residual_vec = force_vec - full_stiffmat.selfadjointView<Eigen::Upper>() * disp_vec; // 内力を再計算
        }
    }
    return; // 最大反復数に達して終了
}
