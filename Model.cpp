#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#define MODEL_PI 3.141592653589793238462643

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/SymGEigsShiftSolver.h>

#include <Eigen/Eigenvalues>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/DenseCholesky.h>

#include "Elements/Elements.h"
#include "Components.h"
#include "LoadComponent.h"

#include "Model.h"
#include "SparseMatrixUtils.h"

FEModel::FEModel()
    : RigidLinkData(std::make_shared<RigidLinks>())
{
}

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

    // 質量ゼロの自由度を縮約
    std::vector<int> shrink_indices, other_indices;
    Eigen::Diagonal mdiag = ma.diagonal();
    for (size_t i = 0; i < mdiag.size(); i++)
    {
        if (mdiag.coeffRef(i) < 0.0000001)
            shrink_indices.push_back(i);
        else
            other_indices.push_back(i);
    }

    Eigen::SparseMatrix<double> k_sha, k_shb, k_shc, k_shd, m_sh;
    SparseMatrixUtils::splitMatrixWithResize(ma, shrink_indices, m_sh);
    SparseMatrixUtils::splitMatrixWithResize(ka, shrink_indices, k_sha, k_shb, k_shd);

    k_shc = k_shb.transpose();
    auto solver_vib = createSolver();
    solver_vib->compute(k_shd);
    Eigen::MatrixXd k_shc_dense(k_shc);
    Eigen::SparseMatrix<double> tmp_mat = (k_shb * solver_vib->solveMulti(k_shc_dense)).sparseView();
    k_sha -= tmp_mat.triangularView<Eigen::Upper>();

    // A_op: 行列 A に対する作用素
    Spectra::SparseSymMatProd<double, Eigen::Upper> A_op(m_sh);
    // B_op: 行列 B に対する作用素
    Spectra::SparseCholesky<double, Eigen::Upper> B_op(k_sha);

    // --- 一般固有値問題の設定 ---
    int mat_size = other_indices.size();
    if (computed_num > mat_size - 1)
        computed_num = mat_size - 1;
    if (computed_num < 1 || mat_size - 2 < computed_num)
        return -1;

    int ncv = 2 * computed_num + 1; // Recommended value
    if (ncv > mat_size) ncv = mat_size;
    Spectra::SymGEigsSolver<Spectra::SparseSymMatProd<double, Eigen::Upper>,
        Spectra::SparseCholesky<double, Eigen::Upper>, Spectra::GEigsMode::Cholesky>
        geigs(A_op, B_op, computed_num, ncv);

    geigs.init();
    int nconv = geigs.compute();

    if (geigs.info() == Spectra::CompInfo::Successful)
    {
        Eigen::MatrixXd u1s = geigs.eigenvectors();
        Eigen::MatrixXd tmp_mat2 = -k_shc * u1s;
        Eigen::MatrixXd u2s = solver_vib->solveMulti(tmp_mat2);

        // 縮小空間（master + free）での固有ベクトルを復元
        int reduced_size = master_dof_num + free_indices.size();
        Eigen::MatrixXd reduced_vectors = Eigen::MatrixXd::Zero(reduced_size, nconv);
        for (size_t j = 0; j < nconv; j++) {
            for (size_t i = 0; i < other_indices.size(); i++)
                reduced_vectors(other_indices[i], j) = u1s(i, j);
            for (size_t i = 0; i < shrink_indices.size(); i++)
                reduced_vectors(shrink_indices[i], j) = u2s(i, j);
        }

        // 全体DOFへの固有ベクトルを構築
        Eigen::MatrixXd eigs_vector = Eigen::MatrixXd::Zero(DOFNum(), nev);

        if (master_dof_num > 0) {
            // RigidLinkがある場合: master DOFをslave DOFに展開
            for (size_t i = 0; i < nconv; i++) {
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
            for (size_t j = 0; j < nconv; j++) {
                for (size_t i = 0; i < other_indices.size(); i++)
                    eigs_vector(all_free[other_indices[i]], j) = u1s(i, j);
                for (size_t i = 0; i < shrink_indices.size(); i++)
                    eigs_vector(all_free[shrink_indices[i]], j) = u2s(i, j);
            }
        }

        for (size_t i = 0; i < nconv; i++)
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

        // 固有値を元の固有値問題に戻す
        for (double v : geigs.eigenvalues())
            eigen_values.push_back(1.0 / sqrt(v));

    }
    else {
        return -1;
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

Eigen::SparseMatrix<double> FEModel::AssembleStiffnessMatrix()
{
    int mat_size = Nodes.size() * 6;

    // 非零要素数を推定
    size_t estimated_nnz = Elements.size() * 100;
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(estimated_nnz);

    // 各要素からTripletを収集
    for (const std::shared_ptr<ElementBase>& eh : Elements) {
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

