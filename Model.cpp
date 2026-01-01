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

// 2x2ブロック行列を単一の対称疎行列に結合
void FEModel::mergeMatrixWithResize(
    const Eigen::SparseMatrix<double>& free_free,
    const Eigen::SparseMatrix<double>& free_fixed,
    const Eigen::SparseMatrix<double>& fixed_fixed,
    Eigen::SparseMatrix<double>& A)
{
    // 各ブロックのサイズから全体サイズを計算
    int size_free = free_free.rows();
    int size_fixed = fixed_fixed.rows();
    int total_size = size_free + size_fixed;

    // インデックスオフセット
    int offset_free = 0;
    int offset_fixed = size_free;

    // Tripletベクトルの準備
    std::vector<Eigen::Triplet<double>> triplets;

    // メモリ最適化: 全ブロックの非零要素数を事前計算
    size_t estimated_nnz = free_free.nonZeros() + fixed_fixed.nonZeros()
                         + 2 * free_fixed.nonZeros();
    triplets.reserve(estimated_nnz);

    // ブロック11 (free_free): オフセット(0, 0)
    // 対角ブロックなので、上三角のみの場合は対称化が必要
    for (int k = 0; k < free_free.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(free_free, k); it; ++it) {
            int row = offset_free + it.row();
            int col = offset_free + it.col();
            triplets.emplace_back(row, col, it.value());

            // 非対角要素は転置も追加
            if (row != col) {
                triplets.emplace_back(col, row, it.value());
            }
        }
    }

    // ブロック12 (free_fixed): オフセット(0, size_free)
    for (int k = 0; k < free_fixed.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(free_fixed, k); it; ++it) {
            int global_row = offset_free + it.row();
            int global_col = offset_fixed + it.col();

            // 上三角要素を追加
            triplets.emplace_back(global_row, global_col, it.value());

            // 対称性: 下三角要素（ブロック21 = fixed_free）も追加
            triplets.emplace_back(global_col, global_row, it.value());
        }
    }

    // ブロック22 (fixed_fixed): オフセット(size_free, size_free)
    // 対角ブロックなので、上三角のみの場合は対称化が必要
    for (int k = 0; k < fixed_fixed.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(fixed_fixed, k); it; ++it) {
            int row = offset_fixed + it.row();
            int col = offset_fixed + it.col();
            triplets.emplace_back(row, col, it.value());

            // 非対角要素は転置も追加
            if (row != col) {
                triplets.emplace_back(col, row, it.value());
            }
        }
    }

    // 疎行列の構築
    A.resize(total_size, total_size);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

// 対称疎行列を3x3ブロックに分割
void FEModel::splitMatrix3x3(
    const Eigen::SparseMatrix<double>& A,
    const std::vector<int>& indices_group1,
    const std::vector<int>& indices_group2,
    Eigen::SparseMatrix<double>& mat_11,
    Eigen::SparseMatrix<double>& mat_12,
    Eigen::SparseMatrix<double>& mat_13,
    Eigen::SparseMatrix<double>& mat_22,
    Eigen::SparseMatrix<double>& mat_23,
    Eigen::SparseMatrix<double>& mat_33)
{
    // グループ分類用の列挙型
    enum class GroupType : uint8_t {
        GROUP_1 = 0,
        GROUP_2 = 1,
        GROUP_3 = 2
    };

    // インデックス属性構造体
    struct idx_attr_three {
        GroupType group;
        int new_index;
    };

    // 高速検索用のunordered_set
    std::unordered_set<int> set1(indices_group1.begin(), indices_group1.end());
    std::unordered_set<int> set2(indices_group2.begin(), indices_group2.end());

    // 各グループのサイズ計算（重複除去済み）
    int size_1 = static_cast<int>(set1.size());
    int size_2 = static_cast<int>(set2.size());
    int size_3 = A.cols() - size_1 - size_2;

    // インデックスマッピングテーブル構築
    std::vector<idx_attr_three> index_map(A.cols());
    int counter_1 = 0, counter_2 = 0, counter_3 = 0;

    for (int i = 0; i < A.cols(); ++i) {
        if (set1.find(i) != set1.end()) {
            index_map[i] = {GroupType::GROUP_1, counter_1++};
        } else if (set2.find(i) != set2.end()) {
            index_map[i] = {GroupType::GROUP_2, counter_2++};
        } else {
            index_map[i] = {GroupType::GROUP_3, counter_3++};
        }
    }

    // 6つのTripletベクトル（上三角6ブロック用）
    std::vector<Eigen::Triplet<double>> trip_11, trip_12, trip_13;
    std::vector<Eigen::Triplet<double>> trip_22, trip_23, trip_33;

    // メモリ最適化: 事前容量確保
    size_t estimated_nnz = A.nonZeros();
    trip_11.reserve(estimated_nnz / 9);
    trip_12.reserve(estimated_nnz / 9);
    trip_13.reserve(estimated_nnz / 9);
    trip_22.reserve(estimated_nnz / 9);
    trip_23.reserve(estimated_nnz / 9);
    trip_33.reserve(estimated_nnz / 9);

    // 疎行列の非零要素を1回だけ走査
    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            double value = it.value();

            const idx_attr_three& attr_i = index_map[i];
            const idx_attr_three& attr_j = index_map[j];

            // 9パターンを判定（対称性を考慮して6パターンのみ処理）
            if (attr_i.group == GroupType::GROUP_1) {
                if (attr_j.group == GroupType::GROUP_1) {
                    trip_11.emplace_back(attr_i.new_index, attr_j.new_index, value);
                } else if (attr_j.group == GroupType::GROUP_2) {
                    trip_12.emplace_back(attr_i.new_index, attr_j.new_index, value);
                } else {
                    trip_13.emplace_back(attr_i.new_index, attr_j.new_index, value);
                }
            } else if (attr_i.group == GroupType::GROUP_2) {
                if (attr_j.group == GroupType::GROUP_1) {
                    // 転置: ブロック21 → ブロック12
                    trip_12.emplace_back(attr_j.new_index, attr_i.new_index, value);
                } else if (attr_j.group == GroupType::GROUP_2) {
                    trip_22.emplace_back(attr_i.new_index, attr_j.new_index, value);
                } else {
                    trip_23.emplace_back(attr_i.new_index, attr_j.new_index, value);
                }
            } else {  // GROUP_3
                if (attr_j.group == GroupType::GROUP_1) {
                    // 転置: ブロック31 → ブロック13
                    trip_13.emplace_back(attr_j.new_index, attr_i.new_index, value);
                } else if (attr_j.group == GroupType::GROUP_2) {
                    // 転置: ブロック32 → ブロック23
                    trip_23.emplace_back(attr_j.new_index, attr_i.new_index, value);
                } else {
                    trip_33.emplace_back(attr_i.new_index, attr_j.new_index, value);
                }
            }
        }
    }

    // 疎行列の再構築
    mat_11.resize(size_1, size_1);
    mat_12.resize(size_1, size_2);
    mat_13.resize(size_1, size_3);
    mat_22.resize(size_2, size_2);
    mat_23.resize(size_2, size_3);
    mat_33.resize(size_3, size_3);

    mat_11.setFromTriplets(trip_11.begin(), trip_11.end());
    mat_12.setFromTriplets(trip_12.begin(), trip_12.end());
    mat_13.setFromTriplets(trip_13.begin(), trip_13.end());
    mat_22.setFromTriplets(trip_22.begin(), trip_22.end());
    mat_23.setFromTriplets(trip_23.begin(), trip_23.end());
    mat_33.setFromTriplets(trip_33.begin(), trip_33.end());
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
    std::vector<int> slave_indices = RigidLinkData.SlaveDOFIndices();
    std::unordered_set<int> slaveid_set(slave_indices.begin(), slave_indices.end());
    std::vector<int> indices;
    int idx = 0;
    for (Node n : Nodes)
        for (const bool f : n.Fix.isdof_fixed()) {
            bool is_slave = (slaveid_set.find(idx) != slaveid_set.end());
            if (!f && !is_slave) indices.push_back(idx);
            idx++;
        }
    return indices;
}

std::vector<int> FEModel::SlaveIndices()
{
    return RigidLinkData.SlaveDOFIndices();
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

    std::vector<int> slave_indices = RigidLinkData.SlaveDOFIndices();
    std::vector<int> free_indices = FreeIndices();
    std::vector<int> fixed_indices = FixIndices();

    Eigen::SparseMatrix<double> m11, m12, m13, m22, m23, m33;
    FEModel::splitMatrix3x3(AssembleStiffnessMatrix(), slave_indices, free_indices,
        m11, m12, m13, m22, m23, m33);

    Eigen::MatrixXd linkTransMat = RigidLinkData.TransformationMatrix();
    Eigen::SparseMatrix<double> m, mb, mc;
    m = linkTransMat.transpose() * m11.selfadjointView<Eigen::Upper>() * linkTransMat;
    mb = linkTransMat.transpose() * m12;
    mc = linkTransMat.transpose() * m13;

    FEModel::mergeMatrixWithResize(m, mb, m22, m);
    

    // FEModel::splitMatrixWithResize(AssembleStiffnessMatrix(), fixed_indices, m, mb, mc);

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

