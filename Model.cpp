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

#include "Element.h"
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

DynamicAnalysis::DynamicAnalysis(std::shared_ptr<FEModel> model, DynamicAccelLoad accel_load, FEDynamicDampInitializer* damp)
    : FEDeformOperator(model), accel_load(accel_load)
{
    if (!damp_initializer) {
        // デフォルトの減衰初期化子を用意
        static FEDynamicStiffDampInitializer defaultDamp;
        damp_initializer = &defaultDamp;
    }
}

bool DynamicAnalysis::Initialize()
{
    // 値を初期化
    current_step = 0;
    current_disp = Eigen::VectorXd::Zero(model->FreeDOFNum());
    current_vel = Eigen::VectorXd::Zero(model->FreeDOFNum());
    current_accel = Eigen::VectorXd::Zero(model->FreeDOFNum());

    // マトリクスの組み立て

	// StiffnessMatrixの組み立て
    free_indices = model->FreeIndices();
    fixed_indices = model->FixIndices();
    FEModel::splitMatrixWithResize(model->AssembleStiffnessMatrix(), fixed_indices, matK_aa, matK_ab, matK_bb);
    FEModel::splitMatrixWithResize(model->AssembleMassMatrix(), fixed_indices, matM_aa, matM_ab, matM_bb);

    // 減衰マトリクスの組み立て
    bool damp_init = damp_initializer->Initialize(this);
	if (!damp_init) {
		std::cerr << "Failed to initialize damping matrix." << std::endl;
		return false;
	}

    // 因数分解しておく
	double dt = accel_load.timestep;
    Eigen::SparseMatrix<double> compute_mat;
	compute_mat = matM_aa + 0.5 * dt * matC_aa + beta * dt * dt * matK_aa;
	solver.compute(compute_mat);

    // Recorder初期化
	energy_recorder.Initialize();

    return true;
}

void DynamicAnalysis::ComputeStep()
{
	// ステップ数が最大に達した場合は終了 
    // Note: 加速度をゼロとして続けるという選択肢もある,,,
	if (current_step >= accel_load.Accels.size()) {
		std::cout << "Dynamic analysis completed." << std::endl;
		return;
	}

	double dt = accel_load.timestep;
    Vector gacc = accel_load.Direction * accel_load.Accels[current_step];

	std::vector<int> free_indices = model->FreeIndices();
	Eigen::VectorXd post_accel0 = Eigen::VectorXd::Zero(free_indices.size());
    for (size_t i = 0; i < free_indices.size(); i++)
    {
		int fi = free_indices[i] % NODE_DOF;
		if (fi == 0) post_accel0[i] = gacc.x;
		else if (fi == 1) post_accel0[i] = gacc.y;
		else if (fi == 2) post_accel0[i] = gacc.z;
    }

    std::vector<int> fixed_indices = model->FixIndices();
    Eigen::VectorXd post_accel0_fixed = Eigen::VectorXd::Zero(fixed_indices.size());
    for (size_t i = 0; i < fixed_indices.size(); i++)
    {
        int fi = fixed_indices[i] % NODE_DOF;
        if (fi == 0) post_accel0_fixed[i] = gacc.x;
        else if (fi == 1) post_accel0_fixed[i] = gacc.y;
        else if (fi == 2) post_accel0_fixed[i] = gacc.z;
    }

    // 次ステップの変位、速度、加速度を取得	
	Eigen::VectorXd post_accel = matM_aa.selfadjointView<Eigen::Upper>() * (-post_accel0)
        - matC_aa.selfadjointView<Eigen::Upper>() * (current_vel + 0.5 * dt * current_accel)
        - matK_aa.selfadjointView<Eigen::Upper>() * (current_disp + dt * current_vel + (0.5 - beta) * dt * dt * current_accel);
    post_accel = solver.solve(post_accel);
	Eigen::VectorXd post_vel = current_vel + 0.5 * (current_accel + post_accel) * dt;
	Eigen::VectorXd post_disp = current_disp + dt * current_vel + (0.5 - beta) * dt * dt * current_accel + beta * dt * dt * post_accel;

    // Update
    current_step++;
	current_accel = post_accel;
	current_vel = post_vel;
	current_disp = post_disp;

    // Compute reaction force
	Eigen::VectorXd rf = matK_ab.transpose() * current_disp + matM_ab.transpose() * current_accel + 
        matC_ab.transpose() * current_vel + matM_bb * post_accel0_fixed;
	Eigen::VectorXd rf_full = Eigen::VectorXd::Zero(model->NodeNum() * 6);
    for (size_t i = 0; i < fixed_indices.size(); i++)
		rf_full(fixed_indices[i]) = rf(i);

	current_react_force.clear();
    for (size_t i = 0; i < model->Nodes.size(); i++)
    {
        if (!model->Nodes[i].Fix.IsAnyFix()) continue;
        int pos = i * 6;
        current_react_force.push_back(NodeLoad(i, rf_full[pos], rf_full[pos + 1], rf_full[pos + 2], rf_full[pos + 3], rf_full[pos + 4], rf_full[pos + 5]));
    }

    if (RecordEnabled) {
        // Sampling
        for (auto sampler : samplers) {
            sampler->Sampling(*this);
        }
        // Recording
        energy_recorder.Record(*this);
    }

}

void DynamicAnalysis::ComputeSteps(int steps)
{
	// ステップ数が最大に達した場合は終了
	if (current_step >= steps) {
		std::cout << "Dynamic analysis completed." << std::endl;
		return;
	}

	// ステップ数分計算
	for (int i = current_step; i < steps; i++)
		ComputeStep();
}

bool DynamicAnalysis::SetDisplacements(std::vector<Displacement> disps)
{
    for (size_t i = 0; i < free_indices.size(); i++) {
		size_t idx = free_indices[i];
		size_t pos = idx % NODE_DOF;
		size_t node_id = idx / NODE_DOF;
		
        if (node_id >= disps.size()) return false;

        current_disp(i) = disps[node_id].displace[pos];
    }
    return true;
}

bool DynamicAnalysis::SetVelocities(std::vector<Displacement> vels)
{
    for (size_t i = 0; i < free_indices.size(); i++) {
        size_t idx = free_indices[i];
        size_t pos = idx % NODE_DOF;
        size_t node_id = idx / NODE_DOF;

        if (node_id >= vels.size()) return false;

        current_vel(i) = vels[node_id].displace[pos];
    }
    return true;
}

bool DynamicAnalysis::SetAccelerations(std::vector<Displacement> accs)
{
    for (size_t i = 0; i < free_indices.size(); i++) {
        size_t idx = free_indices[i];
        size_t pos = idx % NODE_DOF;
        size_t node_id = idx / NODE_DOF;

        if (node_id >= accs.size()) return false;

        current_accel(i) = accs[node_id].displace[pos];
    }
    return true;
}

std::vector<Displacement> DynamicAnalysis::GetDisplacements()
{
    std::vector<Displacement> disp;

    // 変形データ整理
    Eigen::VectorXd d = Eigen::VectorXd::Zero(model->Nodes.size() * 6);
    //Eigen::VectorXd d(Nodes.size() * 6) = Eigen::VectorXd::;
    //d.setZero();
    for (size_t i = 0; i < free_indices.size(); i++) {
        d(free_indices[i]) = current_disp(i);
        //std::cout << current_disp(i);
    }

    for (size_t i = 0; i < model->Nodes.size(); i++)
    {
        int pos = i * 6;
        disp.push_back(Displacement(d[pos], d[pos + 1], d[pos + 2], d[pos + 3], d[pos + 4], d[pos + 5]));
    }

    return disp;
}

std::vector<Displacement> DynamicAnalysis::GetVelocities()
{
    std::vector<Displacement> vel;

    // 変形データ整理
    Eigen::VectorXd d = Eigen::VectorXd::Zero(model->Nodes.size() * 6);
    //Eigen::VectorXd d(Nodes.size() * 6) = Eigen::VectorXd::;
    //d.setZero();
    for (size_t i = 0; i < free_indices.size(); i++)
        d(free_indices[i]) = current_vel(i);

    for (size_t i = 0; i < model->Nodes.size(); i++)
    {
        int pos = i * 6;
        vel.push_back(Displacement(d[pos], d[pos + 1], d[pos + 2], d[pos + 3], d[pos + 4], d[pos + 5]));
    }

    return vel;
}

std::vector<Displacement> DynamicAnalysis::GetAccelerations()
{
    std::vector<Displacement> acc;

    // 変形データ整理
    Eigen::VectorXd d = Eigen::VectorXd::Zero(model->Nodes.size() * 6);
    //Eigen::VectorXd d(Nodes.size() * 6) = Eigen::VectorXd::;
    //d.setZero();
    for (size_t i = 0; i < free_indices.size(); i++)
        d(free_indices[i]) = current_accel(i);

    for (size_t i = 0; i < model->Nodes.size(); i++)
    {
        int pos = i * 6;
        acc.push_back(Displacement(d[pos], d[pos + 1], d[pos + 2], d[pos + 3], d[pos + 4], d[pos + 5]));
    }

    return acc;
}

BeamStressData DynamicAnalysis::GetBeamStress(int eid, double p)
{
    BarElementBase *be = dynamic_cast<BarElementBase *>(model->Elements[eid].get());

    BeamStress b_strs = be->stress(this->GetDisplacements()[be->Nodes[0]->id], this->GetDisplacements()[be->Nodes[1]->id]);
    BeamStressData strs = b_strs.Interpolate(p);
    return strs;
}

PlateStressData DynamicAnalysis::GetPlateStressData(int eid, double xi, double eta)
{
    PlateStressData data;
    std::vector<Displacement> displace = this->GetDisplacements();
    if (model->Elements[eid]->Type() == ElementType::DKT)
    {
        TriPlateElement *el = dynamic_cast<TriPlateElement *>(model->Elements[eid].get());
        data = el->stress(displace[el->Nodes[0]->id], displace[el->Nodes[1]->id],
                          displace[el->Nodes[2]->id], xi, eta);
    }
    else if (model->Elements[eid]->Type() == ElementType::DKQ)
    {
        QuadPlateElement *el = dynamic_cast<QuadPlateElement *>(model->Elements[eid].get());
        data = el->stress(displace[el->Nodes[0]->id], displace[el->Nodes[1]->id],
                          displace[el->Nodes[2]->id], displace[el->Nodes[3]->id], xi, eta);
    }
    return data;
}

Displacement DynamicAnalysis::GetBeamDisplace(int eid, double p)
{
    BeamElement *elm = model->GetBeamElement(eid);
    std::vector<Displacement> displace = this->GetDisplacements();
    Displacement disp = elm->DisplaceAt(
        displace[elm->Nodes[0]->id], displace[elm->Nodes[1]->id], p);

    return disp;
}

bool FEDynamicStiffDampInitializer::Initialize(DynamicAnalysis* analysis)
{
	// 解析モデルの固有振動数を計算
	std::vector<double> eigen_values;
	std::vector<std::vector<Displacement>> mode_vectors;
	int nconv = analysis->model->SolveVibration(1, eigen_values, mode_vectors);
	if (nconv < 0) {
		std::cout << "Eigenvalue calculations did not converge." << std::endl;
		return false;
	}

	natural_angle_velocity = eigen_values[0];
    analysis->matC_aa = analysis->matK_aa * (2.0 * damp_rate / natural_angle_velocity);
    analysis->matC_ab = analysis->matK_ab * (2.0 * damp_rate / natural_angle_velocity);
    analysis->matC_bb = analysis->matK_bb * (2.0 * damp_rate / natural_angle_velocity);

	return true;
}

std::vector<double> FEVibrateResult::ParticipationFactors()
{
    Eigen::SparseMatrix<double> mass_mat= model->AssembleMassMatrix();
    
	std::vector<double> participation_factors;

    for (size_t i = 0; i < ModeNum(); i++)
    {
		Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
		for (size_t j = 0; j < model->NodeNum(); j++)
		{
			int p = j * 6;
			v(p) = mode_vectors[i][j].Dx();
			v(p + 1) = mode_vectors[i][j].Dy();
			v(p + 2) = mode_vectors[i][j].Dz();
			v(p + 3) = mode_vectors[i][j].Rx();
			v(p + 4) = mode_vectors[i][j].Ry();
			v(p + 5) = mode_vectors[i][j].Rz();
		}

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        double beta_i = Eigen::VectorXd::Ones(v.size()).dot(vM) / v.dot(vM);

		participation_factors.push_back(beta_i);
    }

	return participation_factors;
}

std::vector<double> FEVibrateResult::ParticipationFactors(Vector direction)
{
    // MassMatrixの組み立て
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();

    std::vector<double> participation_factors;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();

            f(p) = direction.x;
            f(p + 1) = direction.y;
			f(p + 2) = direction.z;
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        double beta_i = f.dot(vM) / v.dot(vM);

        participation_factors.push_back(beta_i);
    }

    return participation_factors;
}

std::vector<Displacement> FEVibrateResult::ParticipationDirectedFactors()
{
    // MassMatrixの組み立て
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();
    std::vector<Displacement> participation_factors;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        vM /= v.dot(vM);
        Displacement facs;
		for (size_t j = 0; j < model->NodeNum(); j++)
		{
			int p = j * 6;
            // 方向別足し合わせ
            Displacement dj(
				vM(p), vM(p + 1), vM(p + 2),
				vM(p + 3), vM(p + 4), vM(p + 5));
			facs += dj;
		}
        participation_factors.push_back(facs);
    }

    return participation_factors;
}

std::vector<double> FEVibrateResult::EffectiveMassRates()
{
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();

	double total_mass = mass_mat.diagonal().sum();

    std::vector<double> mass_rates;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        double mi = v.dot(vM);
        double beta_i = Eigen::VectorXd::Ones(v.size()).dot(vM) / mi;

        mass_rates.push_back(beta_i * beta_i * mi / total_mass);
    }

    return mass_rates;
}

std::vector<Displacement> FEVibrateResult::EffectiveDirectedMassRates()
{
    // MassMatrixの組み立て
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();
    double total_mass = mass_mat.diagonal().sum();
    std::vector<Displacement> mass_rates;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        double vMv = v.dot(vM);
        Displacement facs;
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            // 方向別足し合わせ
            Displacement dj(
                vM(p), vM(p + 1), vM(p + 2),
                vM(p + 3), vM(p + 4), vM(p + 5));
            facs += dj;
        }

		facs = Displacement(
			facs.Dx() * facs.Dx() / vMv,
			facs.Dy() * facs.Dy() / vMv,
            facs.Dz() * facs.Dz() / vMv,
            facs.Rx() * facs.Rx() / vMv,
            facs.Ry() * facs.Ry() / vMv,
            facs.Rz() * facs.Rz() / vMv);

        facs = facs / (total_mass / 3.0);
        mass_rates.push_back(facs);
    }
    
    return mass_rates;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_responseCQC(ResponseValueType vt)
{
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    std::vector<double> spectrums(periods.size());
    for (size_t i = 0; i < periods.size(); i++){
        if (vt == ResponseValueType::Displacement)
            spectrums[i] = SpectrumFunction->Displacement(periods[i]);
        else if (vt == ResponseValueType::Velocity)
            spectrums[i] = SpectrumFunction->Velocity(periods[i]);
        else // (vt == ResponseValueType::Acceleration)
            spectrums[i] = SpectrumFunction->Acceleration(periods[i]);
    }

    std::vector<Displacement> responses(model->NodeNum());
    for (size_t j = 0; j < part_facs.size(); j++)
    {
        std::vector<Displacement> uj = VibrateResult.ModeVectors()[j];
        for (size_t k = 0; k < part_facs.size(); k++)
        {
            std::vector<Displacement> uk = VibrateResult.ModeVectors()[k];

            double rjk = periods[k] / periods[j];
            double correlation = 8.0 * damping_rate * damping_rate * (1.0 + rjk) * pow(rjk, 1.5) /
                (pow((1.0 - rjk * rjk), 2.0) + 4.0 * damping_rate * damping_rate * rjk * pow(1.0 + rjk, 2.0));
            double fac = spectrums[j] * spectrums[k] * part_facs[j] * part_facs[k] * correlation;

            for (size_t i = 0; i < model->NodeNum(); i++)
            {
                Displacement idj = uj[i];
                Displacement idk = uk[i];

                Displacement d = Displacement(
                    idj.Dx() * idk.Dx(), idj.Dy() * idk.Dy(), idj.Dz() * idk.Dz(),
                    idj.Rx() * idk.Rx(), idj.Ry() * idk.Ry(), idj.Rz() * idk.Rz());

                responses[i] += fac * d;
            }
        }
    }

    // 平方根を取る
    for (size_t j = 0; j < responses.size(); j++) {
        responses[j] = Displacement(
            sqrt(responses[j].Dx()), sqrt(responses[j].Dy()), sqrt(responses[j].Dz()),
            sqrt(responses[j].Rx()), sqrt(responses[j].Ry()), sqrt(responses[j].Rz()));
    }

    return responses;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_responseSRSS(ResponseValueType vt)
{
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    std::vector<double> spectrums(periods.size());
    for (size_t i = 0; i < periods.size(); i++){
        if (vt == ResponseValueType::Displacement)
            spectrums[i] = SpectrumFunction->Displacement(periods[i]);
        else if (vt == ResponseValueType::Velocity)
            spectrums[i] = SpectrumFunction->Velocity(periods[i]);
        else // (vt == ResponseValueType::Acceleration)
            spectrums[i] = SpectrumFunction->Acceleration(periods[i]);
    }

    std::vector<Displacement> responses(model->NodeNum());
    for (size_t i = 0; i < part_facs.size(); i++)
    {
        std::vector<Displacement> mode_vector = VibrateResult.ModeVectors()[i];

        // SRSS Method
        for (size_t j = 0; j < mode_vector.size(); j++)
        {
            Displacement d2 = spectrums[i] * part_facs[i] * mode_vector[j];
            // 変位の二乗和を計算
            d2 = Displacement(d2.Dx() * d2.Dx(), d2.Dy() * d2.Dy(), d2.Dz() * d2.Dz(),
                              d2.Rx() * d2.Rx(), d2.Ry() * d2.Ry(), d2.Rz() * d2.Rz());
            responses[j] += d2;
        }
    }

    // SRSS Methodの結果を平方根で正規化
    for (size_t j = 0; j < responses.size(); j++)
    {
        responses[j] = Displacement(
            sqrt(responses[j].Dx()), sqrt(responses[j].Dy()), sqrt(responses[j].Dz()),
            sqrt(responses[j].Rx()), sqrt(responses[j].Ry()), sqrt(responses[j].Rz()));
    }

    return responses;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_responseABS(ResponseValueType vt)
{
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    std::vector<double> spectrums(periods.size());
    for (size_t i = 0; i < periods.size(); i++){
        if (vt == ResponseValueType::Displacement)
            spectrums[i] = SpectrumFunction->Displacement(periods[i]);
        else if (vt == ResponseValueType::Velocity)
            spectrums[i] = SpectrumFunction->Velocity(periods[i]);
        else // (vt == ResponseValueType::Acceleration)
            spectrums[i] = SpectrumFunction->Acceleration(periods[i]);
    }

    std::vector<Displacement> responses(model->NodeNum());
    for (size_t i = 0; i < part_facs.size(); i++)
    {
        std::vector<Displacement> mode_vector = VibrateResult.ModeVectors()[i];
        // | sd x vector x beta_i |

        // ABS Method
        for (size_t j = 0; j < mode_vector.size(); j++)
        {
            // 変位の絶対値を計算
            Displacement d = spectrums[i] * part_facs[i] * mode_vector[j];
            responses[j] += Displacement(abs(d.Dx()), abs(d.Dy()), abs(d.Dz()),
                                         abs(d.Rx()), abs(d.Ry()), abs(d.Rz()));
        }
    }

    return responses;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_response(ResponseValueType vt)
{
    std::vector<Displacement> responses(model->NodeNum());
    if (MethodType == ResponseSpectrumMethodType::CQC)
        responses = calculate_responseCQC(vt);
    else if (MethodType == ResponseSpectrumMethodType::SRSS)
        responses = calculate_responseSRSS(vt);
    else // MethodType == ResponseSpectrumMethodType::ABS
        responses = calculate_responseABS(vt);

    return responses;
}

ResponseSpectrumMethod::ResponseSpectrumMethod(std::shared_ptr<FEModel> model, 
    FEVibrateResult vibrate_result, Vector direction, IResponseSpectrum* spectrum_function, ResponseSpectrumMethodType type)
    : FEDeformOperator(model), VibrateResult(vibrate_result), SpectrumFunction(spectrum_function), Direction(direction), MethodType(type) {

    Compute();
}

void ResponseSpectrumMethod::Compute()
{
    displacements = calculate_response(ResponseValueType::Displacement);
    velocities = calculate_response(ResponseValueType::Velocity);
    accelerations = calculate_response(ResponseValueType::Acceleration);
    m_computed = true;
}

std::vector<Displacement> ResponseSpectrumMethod::GetDisplacements()
{
    if (m_computed)
        return displacements;
    else
        return calculate_response(ResponseValueType::Displacement);
}

std::vector<Displacement> ResponseSpectrumMethod::GetVelocities()
{
    if (m_computed)
        return velocities;
    else
        return calculate_response(ResponseValueType::Velocity);
}

std::vector<Displacement> ResponseSpectrumMethod::GetAccelerations()
{
    if (m_computed)
        return accelerations;
    else
        return calculate_response(ResponseValueType::Acceleration);
}

BeamStressData ResponseSpectrumMethod::GetBeamStress(int eid, double p)
{
    if (!m_computed)
        throw std::runtime_error("ResponseSpectrumMethod: need to call Compute()");

    BarElementBase* be = dynamic_cast<BarElementBase*>(model->Elements[eid].get());

    BeamStress b_strs = be->stress(this->GetDisplacements()[be->Nodes[0]->id], this->GetDisplacements()[be->Nodes[1]->id]);
    BeamStressData strs = b_strs.Interpolate(p);
    return strs;
}

PlateStressData ResponseSpectrumMethod::GetPlateStressData(int eid, double xi, double eta)
{
    if (!m_computed)
        throw std::runtime_error("ResponseSpectrumMethod: need to call Compute()");
    
    PlateStressData data;
    std::vector<Displacement> displace = this->GetDisplacements();
    if (model->Elements[eid]->Type() == ElementType::DKT)
    {
        TriPlateElement* el = dynamic_cast<TriPlateElement*>(model->Elements[eid].get());
        data = el->stress(displace[el->Nodes[0]->id], displace[el->Nodes[1]->id],
            displace[el->Nodes[2]->id], xi, eta);
    }
    else if (model->Elements[eid]->Type() == ElementType::DKQ)
    {
        QuadPlateElement* el = dynamic_cast<QuadPlateElement*>(model->Elements[eid].get());
        data = el->stress(displace[el->Nodes[0]->id], displace[el->Nodes[1]->id],
            displace[el->Nodes[2]->id], displace[el->Nodes[3]->id], xi, eta);
    }
    return data;
}

Displacement ResponseSpectrumMethod::GetBeamDisplace(int eid, double p)
{
    if (!m_computed)
        throw std::runtime_error("ResponseSpectrumMethod: need to call Compute()");
    
    BeamElement* elm = model->GetBeamElement(eid);
    std::vector<Displacement> displace = this->GetDisplacements();
    Displacement disp = elm->DisplaceAt(
        displace[elm->Nodes[0]->id], displace[elm->Nodes[1]->id], p);

    return disp;
}

FELinearStaticOp ResponseSpectrumMethod::GetLinearStaticCase()
{
	std::vector<std::shared_ptr<LoadBase>> loads;
	std::vector<Displacement> accels = this->GetAccelerations();
    for (size_t i = 0; i < model->NodeNum(); i++)
    {
        // とりあえず並進だけ
		double x = accels[i].Dx() / model->GraityAccel;
		double y = accels[i].Dy() / model->GraityAccel;
		double z = accels[i].Dz() / model->GraityAccel;
        loads.push_back(std::make_shared<NodeBodyForce>( NodeBodyForce(&model->Nodes[i], x, y, z)));
    }
    return FELinearStaticOp(model,loads);
}

void DASampler_MaxDisplacement::Sampling(DynamicAnalysis &da)
{
    bool updated = false;
    std::vector<Displacement> disp = da.GetDisplacements();

    //任意�?�節点の変位が最大となるス�?�?プを記録
    for (int i = 0; i < da.model->Nodes.size(); i++) {
        double d_length = disp[i].Translation().norm();
        if (d_length > max_displacement) {
            max_displacement = d_length;
            step = da.current_step;
            updated = true;
        }
    }

    if (updated) {
        velocity = da.GetVelocities();
        displacement = disp;
        acceleration = da.GetAccelerations();
    }

}

void DAEnergyRecorder::Initialize()
{
    kinetic_energy.clear();
    potential_energy.clear();
    damping_energy.clear();
    input_energy.clear();

	kinetic_energy.push_back(0);
	potential_energy.push_back(0);
	damping_energy.push_back(0);
	input_energy.push_back(0);
}

void DAEnergyRecorder::RecordKineticEnergy(DynamicAnalysis& da)
{
    double energy = 0.5 * da.current_disp.dot(da.matK_aa.selfadjointView<Eigen::Upper>()* da.current_disp);
	kinetic_energy.push_back(energy);
}

void DAEnergyRecorder::RecordPotentialEnergy(DynamicAnalysis& da)
{
    double energy = 0.5 * da.current_vel.dot(da.matM_aa.selfadjointView<Eigen::Upper>() * da.current_vel);
	potential_energy.push_back(energy);
}

void DAEnergyRecorder::RecordDampingEnergy(DynamicAnalysis& da)
{
	double energy = da.current_vel.dot(da.matC_aa.selfadjointView<Eigen::Upper>() * da.current_vel);
	damping_energy.push_back(energy);
}

void DAEnergyRecorder::RecordInputEnergy(DynamicAnalysis& da)
{
    Vector gacc = da.accel_load.Direction * da.accel_load.Accels[da.current_step - 1];
    std::vector<int> free_indices = da.model->FreeIndices();
    Eigen::VectorXd post_accel0 = Eigen::VectorXd::Zero(free_indices.size());
    for (size_t i = 0; i < free_indices.size(); i++)
    {
        int fi = free_indices[i] % NODE_DOF;
        if (fi == 0) post_accel0[i] = gacc.x;
        else if (fi == 1) post_accel0[i] = gacc.y;
        else if (fi == 2) post_accel0[i] = gacc.z;
        else if (fi == 3) post_accel0[i] = 0.0;
        else if (fi == 4) post_accel0[i] = 0.0;
        else if (fi == 5) post_accel0[i] = 0.0;
    }

    Eigen::VectorXd post_accel = da.matM_aa.selfadjointView<Eigen::Upper>() * post_accel0;
    input_energy.push_back(post_accel.dot(da.current_vel));
}

void DAEnergyRecorder::Record(DynamicAnalysis& da)
{
	RecordKineticEnergy(da);
	RecordPotentialEnergy(da);
	RecordDampingEnergy(da);
	RecordInputEnergy(da);
}

int FEBucklingAnalysis::SolveBuckling()
{
    int computed_num = mode_num;
    bool solve_in_spectra = true;
	int solver_selection = 1; // 1: Spectra sparse, 2: Spectra dense, other: Eigen dense

    std::vector<int> free_indices = model->FreeIndices();
    std::vector<int> fixed_indices = model->FixIndices();
    Eigen::SparseMatrix<double> k_full = model->AssembleStiffnessMatrix();
    if (InitailDeformOp != nullptr)
        // 初期変形がある場合は、剛性行列を変形に応じて更新
		k_full += model->AssembleGeometricStiffnessMatrix(InitailDeformOp->GetDisplacements());

    Eigen::SparseMatrix<double> ka; //, kb, kc;
    FEModel::splitMatrixWithResize(k_full, fixed_indices, ka);

    Eigen::SparseMatrix<double> kg; //, kb, kc;
    FEModel::splitMatrixWithResize(
        model->AssembleGeometricStiffnessMatrix(deform_case->GetDisplacements()), 
        fixed_indices, kg);

    if (solver_selection == 1)
    {

        int ncv = 2 * computed_num + 1; // Recommended value

        //using OpType = Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse, Eigen::Upper>;
        //using BOpType = Spectra::SparseSymMatProd<double, Eigen::Upper>;
        //OpType A_op(ka, kg);
        //BOpType B_op(kg);
        //Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::Buckling>
        //   geigs(A_op, B_op, computed_num, ncv, 0.1);

        using OpType = Spectra::SparseSymMatProd<double, Eigen::Upper>;
        using BOpType = Spectra::SparseCholesky<double, Eigen::Upper>;
        OpType A_op(-kg); // Invert
        BOpType B_op(ka); // Invert
        Spectra::SymGEigsSolver<OpType, BOpType, Spectra::GEigsMode::Cholesky> 
             geigs(A_op, B_op, computed_num, ncv);
    
        //  using OpType = Spectra::SparseSymMatProd<double, Eigen::Upper>;
        //  using BOpType = Spectra::SparseRegularInverse<double, Eigen::Upper>;
        //  OpType A_op(kg);
        //  BOpType B_op(ka);
        //  Spectra::SymGEigsSolver<OpType, BOpType, Spectra::GEigsMode::RegularInverse>
        //   geigs(A_op, B_op, computed_num, ncv);

        geigs.init();
        //int nconv = geigs.compute(Spectra::SortRule::LargestMagn);
        //int nconv = geigs.compute(Spectra::SortRule::SmallestMagn);
        // int nconv = geigs.compute(Spectra::SortRule::SmallestAlge);
        int nconv = geigs.compute(Spectra::SortRule::LargestAlge);

        if (geigs.info() == Spectra::CompInfo::Successful)
        {
            Eigen::MatrixXd part_eigen_vectors = geigs.eigenvectors();
            // Eigen::MatrixXd tmp_mat2 = -kg * u1s;
            // Eigen::MatrixXd u2s = solver.solve(tmp_mat2);
            Eigen::MatrixXd eigs_vector = Eigen::MatrixXd::Zero(model->DOFNum(), computed_num);
            for (size_t i = 0; i < free_indices.size(); i++)
            {
                eigs_vector.row(free_indices[i]) = part_eigen_vectors.row(i);
                // for (size_t i = 0; i < other_indices.size(); i++)
                //     eigs_vector.row(free_indices[other_indices[i]]) = u1s.row(i);
                // for (size_t i = 0; i < shrink_indices.size(); i++)
                //     eigs_vector.row(free_indices[shrink_indices[i]]) = u2s.row(i);
            }

            for (size_t i = 0; i < nconv; i++)
            {
                std::vector<Displacement> v(model->NodeNum());
                for (size_t j = 0; j < model->NodeNum(); j++)
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
                eigs.push_back(1.0 / v);
            
        }
        else
        {
            return -1;
        }
        mode_num = nconv;
    }
    else if (solver_selection == 2) {
        // Dense Solver
        Eigen::MatrixXd kg_dense = Eigen::MatrixXd(kg);
        Eigen::MatrixXd k_dense = Eigen::MatrixXd(ka);

        // 行列演算子を作成
        Spectra::DenseSymMatProd<double, Eigen::Upper> kg_dense_op(-kg_dense);
        Spectra::DenseCholesky<double, Eigen::Upper> k_dense_op(k_dense);

        Spectra::SymGEigsSolver<Spectra::DenseSymMatProd<double, Eigen::Upper>,
            Spectra::DenseCholesky<double, Eigen::Upper>,
            Spectra::GEigsMode::Cholesky> solver(kg_dense_op, k_dense_op, computed_num, 2 * computed_num + 1);

        solver.init();
        int nconv = solver.compute(Spectra::SortRule::LargestAlge);

        if (solver.info() == Spectra::CompInfo::Successful)
        {
            Eigen::MatrixXd part_eigen_vectors = solver.eigenvectors();
            Eigen::MatrixXd eigs_vector = Eigen::MatrixXd::Zero(model->DOFNum(), computed_num);
            for (size_t i = 0; i < free_indices.size(); i++)
            {
                eigs_vector.row(free_indices[i]) = part_eigen_vectors.row(i);
                // for (size_t i = 0; i < other_indices.size(); i++)
                //     eigs_vector.row(free_indices[other_indices[i]]) = u1s.row(i);
                // for (size_t i = 0; i < shrink_indices.size(); i++)
                //     eigs_vector.row(free_indices[shrink_indices[i]]) = u2s.row(i);
            }

            for (size_t i = 0; i < nconv; i++)
            {
                std::vector<Displacement> v(model->NodeNum());
                for (size_t j = 0; j < model->NodeNum(); j++)
                {
                    int p = j * 6;
                    v[j] = Displacement(
                        eigs_vector(p, i), eigs_vector(p + 1, i), eigs_vector(p + 2, i),
                        eigs_vector(p + 3, i), eigs_vector(p + 4, i), eigs_vector(p + 5, i));
                }
                mode_vectors.push_back(v);
            }

            // 固有値を元の固有値問題に戻す
            for (double v : solver.eigenvalues())
                eigs.push_back(1.0 / v);

        }
        else
        {
            return -1;
        }
        mode_num = nconv;

    }
    else {
        // 密行列への変換
        Eigen::MatrixXd kg_dense = Eigen::MatrixXd(-kg);
        Eigen::MatrixXd k_dense = Eigen::MatrixXd(ka);

        // 対称行列の場合、下三角部分を補完
        kg_dense.triangularView<Eigen::Lower>() = kg_dense.triangularView<Eigen::Upper>().transpose();
        k_dense.triangularView<Eigen::Lower>() = k_dense.triangularView<Eigen::Upper>().transpose();

        // 一般固有値問題を解く
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(kg_dense, k_dense);

        if (solver.info() == Eigen::Success) {
            Eigen::VectorXd eigs_arr = solver.eigenvalues();
            for (size_t i = 0; i < mode_num; i++)
            {
				double eig = eigs_arr[eigs_arr.size() - i - 1];
                eigs.push_back(1.0 / eig);
            }

            //for (size_t i = 0; i < mode_num; i++)
            //{
            //    double eig = eigs_arr[eigs_arr.size() - i - 1];
            //    eigs.push_back(1.0 / eig);
            //}

   //         for (double v : solver.eigenvalues()) {
   //             eigs.push_back(1.0 / v);
			//}

            Eigen::MatrixXd part_eigen_vectors = solver.eigenvectors();
            // Eigen::MatrixXd tmp_mat2 = -kg * u1s;
            // Eigen::MatrixXd u2s = solver.solve(tmp_mat2);
            Eigen::MatrixXd eigs_vector = Eigen::MatrixXd::Zero(model->DOFNum(), mode_num);
            for (size_t j = 0; j < mode_num; j++)
            {
                for (size_t i = 0; i < free_indices.size(); i++)
                    eigs_vector(free_indices[i], j) = part_eigen_vectors(i, j);
            }
        }
        else {
            return -1;
        }
    }

    
    return mode_num;
}

BeamStressData LinearStaticCombinationOperator::GetBeamStress(int eid, double p)
{
    BeamStressData data;
    for (LinearStaticDeformFactor op : this->cases){
        if (!op.op->Computed())
            throw std::runtime_error("FELinearStaticOP: need to call compute()");

        data += op.factor * op.op->GetBeamStress(eid, p);
    }
    return data;
}

PlateStressData LinearStaticCombinationOperator::GetPlateStressData(int eid, double xi, double eta)
{
    PlateStressData data;
    for (LinearStaticDeformFactor op : this->cases) {
        if (!op.op->Computed())
            throw std::runtime_error("FELinearStaticOP: need to call compute()");

        data += op.factor * op.op->GetPlateStressData(eid, xi, eta);
    }

    return data;
}

Displacement LinearStaticCombinationOperator::GetBeamDisplace(int eid, double p)
{
	Displacement disp;
    for (LinearStaticDeformFactor op : this->cases) {
        if (!op.op->Computed())
            throw std::runtime_error("FELinearStaticOP: need to call compute()");

		disp += op.factor * op.op->GetBeamDisplace(eid, p);
    }
    return disp;
}

std::vector<Displacement> LinearStaticCombinationOperator::GetDisplacements()
{
	std::vector<Displacement> disp(model->NodeNum());
	for (LinearStaticDeformFactor op : this->cases)
	{
        if (!op.op->Computed())
            throw std::runtime_error("FELinearStaticOP: need to call compute()");

		std::vector<Displacement> d = op.op->GetDisplacements();
		for (size_t i = 0; i < d.size(); i++)
			disp[i] += d[i] * op.factor;
	}
    return disp;
}

std::vector<NodeLoad> LinearStaticCombinationOperator::GetReactForces()
{
    std::vector<NodeLoad> combined_react;
    for (const auto& case_factor : cases) {
        if (!case_factor.op->Computed())
            throw std::runtime_error("FELinearStaticOP: need to call compute()");

        auto react = case_factor.op->GetReactForces();
        for (size_t i = 0; i < react.size(); i++)
        {
            bool need_fallback = false;
            // compine_reactがi番目要素を含む場合はi番目要素をチェック
            if (i < combined_react.size()) {
                // i番目要素が対象のノードか？
                if (combined_react[i].id == react[i].id) {
                    combined_react[i].data += (case_factor.factor * react[i].data);
                    continue;
                }
            }

            for (auto& cr : combined_react) {
                if (cr.id == react[i].id) {
                    cr.data += (case_factor.factor * react[i].data);
                    break;
                }
                need_fallback = true;
            }

            if (need_fallback)
                combined_react.push_back(react[i]);
        }
    }
    return combined_react;
}

void FELinearStaticOp::Compute()
{
	this->model->SolveLinearStatic(this->loads, this->displace, this->react_force);
	m_computed = true;
}

BeamStressData FELinearStaticOp::GetBeamStress(int eid, double p)
{
    if (!m_computed)
		throw std::runtime_error("FELinearStaticOP: need to call compute()");

    BarElementBase* be = dynamic_cast<BarElementBase*>(model->Elements[eid].get());

    // BeamElement* elm = GetBeamElement(eid);
    BeamStress b_strs = be->stress(displace[be->Nodes[0]->id], displace[be->Nodes[1]->id]);
    BeamStressData strs = b_strs.Interpolate(p);

    // Add beam stresses due to beam loads
    BeamElement* beamElement = dynamic_cast<BeamElement*>(be);
    if (beamElement != nullptr)
    {
        for (const auto& l : loads)
        {
            // BeamLoadBase* bpl = dynamic_cast<BeamLoadBase*>(l);
            std::shared_ptr<BeamLoadBase> bpl = std::dynamic_pointer_cast<BeamLoadBase>(l);
            if (bpl == NULL)
                continue;
            if (bpl->element->id != beamElement->id)
                continue;
            BeamStressData bsd = bpl->GetBeamStress(p);
            strs.Nx += bsd.Nx;
            strs.My += bsd.My;
            strs.Mz += bsd.Mz;
            strs.Mx += bsd.Mx;
            strs.Qy += bsd.Qy;
            strs.Qz += bsd.Qz;
        }
    }
    return strs;
}

PlateStressData FELinearStaticOp::GetPlateStressData(int eid, double xi, double eta)
{
    if (!m_computed)
        throw std::runtime_error("FELinearStaticOP: need to call compute()");

    PlateStressData data;
    if (model->Elements[eid]->Type() == ElementType::DKT) {
        TriPlateElement* el = dynamic_cast<TriPlateElement*>(model->Elements[eid].get());
        data = el->stress(displace[el->Nodes[0]->id], displace[el->Nodes[1]->id],
            displace[el->Nodes[2]->id], xi, eta);
    }
    else if (model->Elements[eid]->Type() == ElementType::DKQ) {
        QuadPlateElement* el = dynamic_cast<QuadPlateElement*>(model->Elements[eid].get());
        data = el->stress(displace[el->Nodes[0]->id], displace[el->Nodes[1]->id],
            displace[el->Nodes[2]->id], displace[el->Nodes[3]->id], xi, eta);
    }
    return data;
}

Displacement FELinearStaticOp::GetBeamDisplace(int eid, double p)
{
    if (!m_computed)
        throw std::runtime_error("FELinearStaticOP: need to call compute()");

    BeamElement* elm = model->GetBeamElement(eid);
    Displacement disp = elm->DisplaceAt(
        displace[elm->Nodes[0]->id], displace[elm->Nodes[1]->id], p);

    //BeamStressData strs = b_strs.Interpolate(p);
    for (const auto& l : loads)
    {
        std::shared_ptr<BeamLoadBase> bpl = std::dynamic_pointer_cast<BeamLoadBase>(l);
        if (bpl == NULL) continue;
        if (bpl->element->id != elm->id) continue;

        // Not Implemented
        if (bpl->axis == BeamLoadAxis::XAxis) continue;

        Displacement disp_i = bpl->GetDisplacement(p);
        disp = Displacement(
            disp.Dx() + disp_i.Dx(),
            disp.Dy() + disp_i.Dy(),
            disp.Dz() + disp_i.Dz(),
            disp.Rx() + disp_i.Rx(),
            disp.Ry() + disp_i.Ry(),
            disp.Rz() + disp_i.Rz()
        );
    }
    return disp;
}

std::vector<Displacement> FELinearStaticOp::GetDisplacements()
{
    if (!m_computed)
        throw std::runtime_error("FELinearStaticOP: need to call compute()");

    return this->displace;
}

std::vector<NodeLoad> FELinearStaticOp::GetReactForces()
{
    if (!m_computed)
        throw std::runtime_error("FELinearStaticOP: need to call compute()");

    return this->react_force;
}
