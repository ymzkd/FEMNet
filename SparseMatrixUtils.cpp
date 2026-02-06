#include "SparseMatrixUtils.h"
#include <unordered_set>
#include <stdexcept>
#include <iostream>

// 3-parameter版: 4-parameter版を呼び出すラッパー（リファクタリング版）
void SparseMatrixUtils::splitMatrixWithResize(
    const Eigen::SparseMatrix<double>& A,
    const std::vector<int>& fixed_indices,
    Eigen::SparseMatrix<double>& free_matrix)
{
    // 4-parameter版を呼び出し、不要な出力は破棄
    Eigen::SparseMatrix<double> free_fixed_matrix, fixed_matrix;
    splitMatrixWithResize(A, fixed_indices, free_matrix, free_fixed_matrix, fixed_matrix);
}

// 4-parameter版: Model.cpp lines 40-98からコピー
void SparseMatrixUtils::splitMatrixWithResize(
    const Eigen::SparseMatrix<double>& A,
    const std::vector<int>& fixed_indices,
    Eigen::SparseMatrix<double>& free_matrix,
    Eigen::SparseMatrix<double>& free_fixed_matrix,
    Eigen::SparseMatrix<double>& fixed_matrix)
{
    typedef std::pair<bool, int> idx_attr;

    // 出力行列のTripletを準備
    std::vector<Eigen::Triplet<double>> free_triplets;
    std::vector<Eigen::Triplet<double>> free_fixed_triplets;
    std::vector<Eigen::Triplet<double>> fixed_triplets;

    // インデックスマッピングテーブル構築
    // 入力順序を保持するため、fixed_indicesを先に処理
    std::vector<idx_attr> modified_indices(A.cols());
    std::unordered_set<int> used_indices;
    int ifixed = 0, ifree = 0;

    // Fixed indices: fixed_indicesの順序を保持
    for (int idx : fixed_indices) {
        if (idx >= 0 && idx < A.cols() && used_indices.find(idx) == used_indices.end()) {
            modified_indices[idx] = std::make_pair(true, ifixed++);
            used_indices.insert(idx);
        }
    }

    // Free indices: 残りのインデックス（昇順で処理）
    for (int i = 0; i < A.cols(); ++i) {
        if (used_indices.find(i) == used_indices.end()) {
            modified_indices[i] = std::make_pair(false, ifree++);
        }
    }

    // サイズ（重複除去済み）
    int fixed_size = ifixed;
    int free_size = ifree;

    // 行列Aを走査
    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            int i = it.row();    // 行
            int j = it.col();    // 列
            double value = it.value();

            idx_attr idx_attr_i = modified_indices[i];
            idx_attr idx_attr_j = modified_indices[j];
            if (!idx_attr_i.first && !idx_attr_j.first) {
                // 両方自由 - 対角ブロックなので上三角のみ保持
                if (idx_attr_i.second <= idx_attr_j.second) {
                    free_triplets.emplace_back(idx_attr_i.second, idx_attr_j.second, value);
                } else {
                    free_triplets.emplace_back(idx_attr_j.second, idx_attr_i.second, value);
                }
            } else if (!idx_attr_i.first && idx_attr_j.first) {
                // 行自由, 列固定 - 非対角ブロック
                free_fixed_triplets.emplace_back(idx_attr_i.second, idx_attr_j.second, value);
            } else if (idx_attr_i.first && !idx_attr_j.first) {
                // 行固定, 列自由 - ブロック転置して非対角ブロックへ
                free_fixed_triplets.emplace_back(idx_attr_j.second, idx_attr_i.second, value);
            } else if (idx_attr_i.first && idx_attr_j.first) {
                // 両方固定 - 対角ブロックなので上三角のみ保持
                if (idx_attr_i.second <= idx_attr_j.second) {
                    fixed_triplets.emplace_back(idx_attr_i.second, idx_attr_j.second, value);
                } else {
                    fixed_triplets.emplace_back(idx_attr_j.second, idx_attr_i.second, value);
                }
            }
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

// splitMatrix3x3: Model.cpp lines 217-333からコピー
void SparseMatrixUtils::splitMatrix3x3(
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

    // インデックスマッピングテーブル構築
    // 入力順序を保持するため、入力ベクトルを順に処理
    std::vector<idx_attr_three> index_map(A.cols());
    std::unordered_set<int> used_indices;
    int counter_1 = 0, counter_2 = 0, counter_3 = 0;

    // Group 1: indices_group1の順序を保持
    for (int idx : indices_group1) {
        if (idx >= 0 && idx < A.cols() && used_indices.find(idx) == used_indices.end()) {
            index_map[idx] = {GroupType::GROUP_1, counter_1++};
            used_indices.insert(idx);
        }
    }

    // Group 2: indices_group2の順序を保持
    for (int idx : indices_group2) {
        if (idx >= 0 && idx < A.cols() && used_indices.find(idx) == used_indices.end()) {
            index_map[idx] = {GroupType::GROUP_2, counter_2++};
            used_indices.insert(idx);
        }
    }

    // Group 3: 残りのインデックス（昇順で処理）
    for (int i = 0; i < A.cols(); ++i) {
        if (used_indices.find(i) == used_indices.end()) {
            index_map[i] = {GroupType::GROUP_3, counter_3++};
        }
    }

    // 各グループのサイズ（重複除去済み）
    int size_1 = counter_1;
    int size_2 = counter_2;
    int size_3 = counter_3;

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
            // 対角ブロックは上三角のみ保持、非対角ブロックはブロック間転置のみ
            if (attr_i.group == GroupType::GROUP_1) {
                if (attr_j.group == GroupType::GROUP_1) {
                    // 対角ブロック11: 上三角のみ保持
                    if (attr_i.new_index <= attr_j.new_index) {
                        trip_11.emplace_back(attr_i.new_index, attr_j.new_index, value);
                    } else {
                        trip_11.emplace_back(attr_j.new_index, attr_i.new_index, value);
                    }
                } else if (attr_j.group == GroupType::GROUP_2) {
                    // 非対角ブロック12: そのまま配置
                    trip_12.emplace_back(attr_i.new_index, attr_j.new_index, value);
                } else {
                    // 非対角ブロック13: そのまま配置
                    trip_13.emplace_back(attr_i.new_index, attr_j.new_index, value);
                }
            } else if (attr_i.group == GroupType::GROUP_2) {
                if (attr_j.group == GroupType::GROUP_1) {
                    // 転置: ブロック21 → ブロック12
                    trip_12.emplace_back(attr_j.new_index, attr_i.new_index, value);
                } else if (attr_j.group == GroupType::GROUP_2) {
                    // 対角ブロック22: 上三角のみ保持
                    if (attr_i.new_index <= attr_j.new_index) {
                        trip_22.emplace_back(attr_i.new_index, attr_j.new_index, value);
                    } else {
                        trip_22.emplace_back(attr_j.new_index, attr_i.new_index, value);
                    }
                } else {
                    // 非対角ブロック23: そのまま配置
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
                    // 対角ブロック33: 上三角のみ保持
                    if (attr_i.new_index <= attr_j.new_index) {
                        trip_33.emplace_back(attr_i.new_index, attr_j.new_index, value);
                    } else {
                        trip_33.emplace_back(attr_j.new_index, attr_i.new_index, value);
                    }
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

// mergeMatrixWithResize: Model.cpp lines 144-214からコピー
void SparseMatrixUtils::mergeMatrixWithResize(
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

    // メモリ最適化: 全ブロックの非零要素数を事前計算（上三角のみ）
    size_t estimated_nnz = free_free.nonZeros() + fixed_fixed.nonZeros()
                         + free_fixed.nonZeros();
    triplets.reserve(estimated_nnz);

    // ブロック11 (free_free): オフセット(0, 0)
    // 上三角のみ保持（対称行列として扱う）
    for (int k = 0; k < free_free.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(free_free, k); it; ++it) {
            int row = offset_free + it.row();
            int col = offset_free + it.col();
            triplets.emplace_back(row, col, it.value());
        }
    }

    // ブロック12 (free_fixed): オフセット(0, size_free)
    // 上三角要素のみ追加
    for (int k = 0; k < free_fixed.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(free_fixed, k); it; ++it) {
            int global_row = offset_free + it.row();
            int global_col = offset_fixed + it.col();
            triplets.emplace_back(global_row, global_col, it.value());
        }
    }

    // ブロック22 (fixed_fixed): オフセット(size_free, size_free)
    // 上三角のみ保持（対称行列として扱う）
    for (int k = 0; k < fixed_fixed.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(fixed_fixed, k); it; ++it) {
            int row = offset_fixed + it.row();
            int col = offset_fixed + it.col();
            triplets.emplace_back(row, col, it.value());
        }
    }

    // 疎行列の構築
    A.resize(total_size, total_size);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

// vstack: 新規実装（縦方向結合）
Eigen::SparseMatrix<double> SparseMatrixUtils::vstack(
    const Eigen::SparseMatrix<double>& A,
    const Eigen::SparseMatrix<double>& B)
{
    // サイズチェック（列数が一致する必要がある）
    if (A.cols() != B.cols()) {
        throw std::invalid_argument("vstack: column sizes must match");
    }

    int rows_A = A.rows();
    int rows_B = B.rows();
    int cols = A.cols();
    int total_rows = rows_A + rows_B;

    // Tripletベクトルの準備
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(A.nonZeros() + B.nonZeros());

    // 上側ブロック (A): オフセット行 = 0
    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            triplets.emplace_back(it.row(), it.col(), it.value());
        }
    }

    // 下側ブロック (B): オフセット行 = rows_A
    for (int k = 0; k < B.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(B, k); it; ++it) {
            triplets.emplace_back(rows_A + it.row(), it.col(), it.value());
        }
    }

    // 疎行列の構築
    Eigen::SparseMatrix<double> result(total_rows, cols);
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

// hstack: 新規実装（横方向結合）
Eigen::SparseMatrix<double> SparseMatrixUtils::hstack(
    const Eigen::SparseMatrix<double>& A,
    const Eigen::SparseMatrix<double>& B)
{
    // サイズチェック（行数が一致する必要がある）
    if (A.rows() != B.rows()) {
        throw std::invalid_argument("hstack: row sizes must match");
    }

    int rows = A.rows();
    int cols_A = A.cols();
    int cols_B = B.cols();
    int total_cols = cols_A + cols_B;

    // Tripletベクトルの準備
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(A.nonZeros() + B.nonZeros());

    // 左側ブロック (A): オフセット列 = 0
    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            triplets.emplace_back(it.row(), it.col(), it.value());
        }
    }

    // 右側ブロック (B): オフセット列 = cols_A
    for (int k = 0; k < B.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(B, k); it; ++it) {
            triplets.emplace_back(it.row(), cols_A + it.col(), it.value());
        }
    }

    // 疎行列の構築
    Eigen::SparseMatrix<double> result(rows, total_cols);
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}
