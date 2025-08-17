#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#define MODEL_PI 3.141592653589793238462643

#ifdef EIGEN_USE_MKL_ALL
    //#define EIGEN_USE_MKL_ALL
    #include <Eigen/Sparse>
    #include <Eigen/PardisoSupport>
    #include <Eigen/SparseCholesky>
    #include <Spectra/MatOp/SparseSymMatProd.h>
    #include <Spectra/MatOp/SparseCholesky.h>
    #include <Spectra/MatOp/SparseSymShiftSolve.h>
    #include <Spectra/SymGEigsSolver.h>
    #include <Spectra/SymGEigsShiftSolver.h>
#else
    #include <Eigen/Sparse>
    #include <Spectra/MatOp/SparseSymMatProd.h>
    #include <Spectra/MatOp/SparseCholesky.h>
    #include <Spectra/MatOp/SparseSymShiftSolve.h>
    #include <Spectra/SymGEigsSolver.h>
    #include <Spectra/SymGEigsShiftSolver.h>
    // #include <Spectra/Util/CompInfo.h>
#endif

#include <Eigen/Eigenvalues>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/DenseCholesky.h>

#include "Elements/Elements.h"
#include "Components.h"
#include "LoadComponent.h"

#include "Model.h"


// 入力: 対称なSparseMatrix、行インデックス配列、列インデックス配列
void FEModel::splitMatrixWithResize(
    const Eigen::SparseMatrix<double>& A,
    const std::vector<int>& fixed_indices,
    Eigen::SparseMatrix<double>& free_matrix,
    Eigen::SparseMatrix<double>& free_fixed_matrix,
    Eigen::SparseMatrix<double>& fixed_matrix)
{
    typedef std::pair<bool, int> idx_attr;
    // インデックスセット（高速検索用）
    std::unordered_set<int> fixed_set(fixed_indices.begin(), fixed_indices.end());

    // 出力行列のTripletを準備
    std::vector<Eigen::Triplet<double>> free_triplets;
    std::vector<Eigen::Triplet<double>> free_fixed_triplets;
    std::vector<Eigen::Triplet<double>> fixed_triplets;

    int fixed_size = fixed_indices.size();
    int free_size = A.cols() - fixed_size;
    std::vector<idx_attr> modified_indices(A.cols());
    int ifixed = 0, ifree = 0;
    for (int i = 0; i < A.cols(); ++i) {
        bool is_fixed = fixed_set.find(i) != fixed_set.end();
        //bool is_fixed = fixed_set.contains(i);
        modified_indices[i] = std::make_pair(is_fixed, is_fixed ? ifixed++ : ifree++);
    }

    // 行列Aを走査
    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            int i = it.row();    // 行
            int j = it.col();    // 列
            double value = it.value();

            idx_attr idx_attr_i = modified_indices[i];
            idx_attr idx_attr_j = modified_indices[j];
            if (!idx_attr_i.first && !idx_attr_j.first)
                // 両方自由
                free_triplets.emplace_back(idx_attr_i.second, idx_attr_j.second, value);
            else if (!idx_attr_i.first && idx_attr_j.first)
                // 行自由, 列固定
                free_fixed_triplets.emplace_back(idx_attr_i.second, idx_attr_j.second, value);
            else if (idx_attr_i.first && !idx_attr_j.first)
                // 行固定, 列自由
                free_fixed_triplets.emplace_back(idx_attr_j.second, idx_attr_i.second, value);
            else if (idx_attr_i.first && idx_attr_j.first)
                // 両方固定
                fixed_triplets.emplace_back(idx_attr_i.second, idx_attr_j.second, value);
        }
    }

    // リサイズと構築
    free_matrix.resize(free_size, free_size);
    free_fixed_matrix.resize(free_size, fixed_size);
    fixed_matrix.resize(fixed_size, fixed_size);

    free_matrix.setFromTriplets(free_triplets.begin(), free_triplets.end());
    free_fixed_matrix.setFromTriplets(free_fixed_triplets.begin(), free_fixed_triplets.end());
    fixed_matrix.setFromTriplets(fixed_triplets.begin(), fixed_triplets.end());
}

// 入力: 対称なSparseMatrix、行インデックス配列、列インデックス配列
void FEModel::splitMatrixWithResize(
    const Eigen::SparseMatrix<double>& A,
    const std::vector<int>& fixed_indices,
    Eigen::SparseMatrix<double>& free_matrix)
{
    typedef std::pair<bool, int> idx_attr;
    // インデックスセット（高速検索用）
    std::unordered_set<int> fixed_set(fixed_indices.begin(), fixed_indices.end());

    // 出力行列のTripletを準備
    std::vector<Eigen::Triplet<double>> free_triplets;

    int fixed_size = fixed_indices.size();
    int free_size = A.cols() - fixed_size;
    std::vector<idx_attr> modified_indices(A.cols());
    int ifixed = 0, ifree = 0;
    for (int i = 0; i < A.cols(); ++i) {
        bool is_fixed = fixed_set.find(i) != fixed_set.end();
        //bool is_fixed = fixed_set.contains(i);
        modified_indices[i] = std::make_pair(is_fixed, is_fixed ? ifixed++ : ifree++);
    }

    // 行列Aを走査
    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            int i = it.row();    // 行
            int j = it.col();    // 列
            double value = it.value();

            idx_attr idx_attr_i = modified_indices[i];
            idx_attr idx_attr_j = modified_indices[j];
            if (!idx_attr_i.first && !idx_attr_j.first)
                // 両方自由
                free_triplets.emplace_back(idx_attr_i.second, idx_attr_j.second, value);
        }
    }

    // リサイズと構築
    free_matrix.resize(free_size, free_size);
    free_matrix.setFromTriplets(free_triplets.begin(), free_triplets.end());
}

Eigen::SparseMatrix<double> FEModel::extractSubMatrix(
    const Eigen::SparseMatrix<double>& mat,
    const std::vector<int>& rowIndices,
    const std::vector<int>& colIndices) {

    // 新しい疎行列を作成
    Eigen::SparseMatrix<double> subMat(rowIndices.size(), colIndices.size());
    

    // インデックス変換用のマップ
    std::unordered_map<int, int> rowMap;
    for (int i = 0; i < rowIndices.size(); ++i)
        rowMap[rowIndices[i]] = i;

    std::unordered_map<int, int> colMap;
    for (int j = 0; j < colIndices.size(); ++j)
        colMap[colIndices[j]] = j;

    // 指定された行と列の要素を新しい疎行列にコピー
    for (int k = 0; k < mat.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
            auto rowMapItem = rowMap.find(it.row());
            auto colMapItem = colMap.find(it.col());

            if (rowMapItem != rowMap.end() && colMapItem != colMap.end())
                subMat.insert(rowMapItem->second, colMapItem->second) = it.value();
        }
    }

    return subMat;
}

int FEModel::SolveVibration(const int nev, std::vector<double>& eigen_values, 
    std::vector<std::vector<Displacement>>& mode_vectors)
{
    int computed_num = nev;

    std::vector<int> free_indices = FreeIndices();
    std::vector<int> fixed_indices = FixIndices();
    Eigen::SparseMatrix<double> ka; //, kb, kc;
    FEModel::splitMatrixWithResize(AssembleStiffnessMatrix(), fixed_indices, ka);
    Eigen::SparseMatrix<double> ma;
    FEModel::splitMatrixWithResize(AssembleMassMatrix(), fixed_indices, ma);

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
    FEModel::splitMatrixWithResize(ma, shrink_indices, m_sh);
    FEModel::splitMatrixWithResize(ka, shrink_indices, k_sha, k_shb, k_shd);

#ifdef EIGEN_USE_MKL_ALL
    Eigen::PardisoLLT<Eigen::SparseMatrix<double>> solver;
#else
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
#endif

    k_shc = k_shb.transpose();
    solver.compute(k_shd);
    Eigen::SparseMatrix<double> tmp_mat = k_shb * solver.solve(k_shc);
    k_sha -= tmp_mat.triangularView<Eigen::Upper>();

    // A_op: 行列 A に対する作用素
    Spectra::SparseSymMatProd<double, Eigen::Upper> A_op(m_sh);
    // B_op: 行列 B に対する作用素
    Spectra::SparseCholesky<double, Eigen::Upper> B_op(k_sha);

    // --- 一般固有値問題の設定 ---
    // 求める固有値の個数 (nev) と、アルゴリズム内部で使用する次元 (ncv) を指定します
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
        Eigen::MatrixXd u2s = solver.solve(tmp_mat2);
        Eigen::MatrixXd eigs_vector = Eigen::MatrixXd::Zero(DOFNum(), nev);
        for (size_t i = 0; i < free_indices.size(); i++)
        {
            for (size_t i = 0; i < other_indices.size(); i++)
                eigs_vector.row(free_indices[other_indices[i]]) = u1s.row(i);
            for (size_t i = 0; i < shrink_indices.size(); i++)
                eigs_vector.row(free_indices[shrink_indices[i]]) = u2s.row(i);
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

std::vector<int> FEModel::FreeIndices()
{
    std::vector<int> indices;
    int idx = 0;
    for (Node n : Nodes)
        for (const bool f : n.Fix.isdof_fixed()) {
            if (!f) indices.push_back(idx);
            idx++;
        }
    return indices;
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
    Eigen::SparseMatrix<double> mat(mat_size, mat_size);
    for (const std::shared_ptr<ElementBase>& eh : Elements)
        eh->AssembleStiffMatrix(mat);

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
    Eigen::SparseMatrix<double> mat(mat_size, mat_size);
    for (const std::shared_ptr<ElementBase>& eh : Elements) {

        std::vector<Displacement> disp_vec(eh->NodeNum());
        std::vector<Node*> nodes = eh->NodesList();
        for (size_t i = 0; i < eh->NodeNum(); i++)
            disp_vec[i] = displacements[nodes[i]->id];
        eh->AssembleGeometricStiffMatrix(mat, disp_vec);
    }
        
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

        for (auto& nl : node_loads) {
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

    std::vector<int> free_indices = FreeIndices();
    std::vector<int> fixed_indices = FixIndices();
    Eigen::SparseMatrix<double> m, mb, mc;
    FEModel::splitMatrixWithResize(AssembleStiffnessMatrix(), fixed_indices, m, mb, mc);

    Eigen::VectorXd f_free(free_indices.size());
    Eigen::VectorXd f_fix(fixed_indices.size());
    for (size_t i = 0; i < free_indices.size(); i++)
        f_free(i) = f(free_indices[i]);
    for (size_t i = 0; i < fixed_indices.size(); i++)
        f_fix(i) = f(fixed_indices[i]);

#ifdef EIGEN_USE_MKL_ALL
    Eigen::PardisoLLT<Eigen::SparseMatrix<double>> solver;
#else
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
#endif

    // Solve
    solver.compute(m);
    Eigen::VectorXd d_free = solver.solve(f_free);
    Eigen::VectorXd r_fix = mb.transpose() * d_free - f_fix;

    // 反力データ整理
    Eigen::VectorXd r = Eigen::VectorXd::Zero(Nodes.size() * 6);
    for (size_t i = 0; i < fixed_indices.size(); i++)
        r(fixed_indices[i]) = r_fix(i);

    // std::vector<NodeLoad> react;
    for (size_t i = 0; i < Nodes.size(); i++)
    {
        if (!Nodes[i].Fix.IsAnyFix()) continue;
        int pos = i * 6;
        react.push_back(NodeLoad(i, r[pos], r[pos + 1], r[pos + 2], r[pos + 3], r[pos + 4], r[pos + 5]));
    }

    // 変形データ整理
    Eigen::VectorXd d = Eigen::VectorXd::Zero(Nodes.size() * 6);
    //Eigen::VectorXd d(Nodes.size() * 6) = Eigen::VectorXd::;
    //d.setZero();
    for (size_t i = 0; i < free_indices.size(); i++)
        d(free_indices[i]) = d_free(i);

    // std::vector<Displacement> disp;
    for (size_t i = 0; i < Nodes.size(); i++)
    {
        int pos = i * 6;
        disp.push_back(Displacement(d[pos], d[pos + 1], d[pos + 2], d[pos + 3], d[pos + 4], d[pos + 5]));
    }
}

