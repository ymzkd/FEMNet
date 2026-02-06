#ifndef _SPARSE_MATRIX_UTILS_
#define _SPARSE_MATRIX_UTILS_

#include <vector>
#include <Eigen/Sparse>

/// <summary>
/// 疎行列操作のユーティリティ関数群
/// 内部実装専用（SWIG非公開）
/// </summary>
class SparseMatrixUtils {
public:
    // 行列分割関数（2x2分割）
    static void splitMatrixWithResize(
        const Eigen::SparseMatrix<double>& A,
        const std::vector<int>& fixed_indices,
        Eigen::SparseMatrix<double>& free_matrix);

    static void splitMatrixWithResize(
        const Eigen::SparseMatrix<double>& A,
        const std::vector<int>& fixed_indices,
        Eigen::SparseMatrix<double>& free_matrix,
        Eigen::SparseMatrix<double>& free_fixed_matrix,
        Eigen::SparseMatrix<double>& fixed_matrix);

    // 行列分割関数（3x3分割）
    static void splitMatrix3x3(
        const Eigen::SparseMatrix<double>& A,
        const std::vector<int>& indices_group1,
        const std::vector<int>& indices_group2,
        Eigen::SparseMatrix<double>& mat_11,
        Eigen::SparseMatrix<double>& mat_12,
        Eigen::SparseMatrix<double>& mat_13,
        Eigen::SparseMatrix<double>& mat_22,
        Eigen::SparseMatrix<double>& mat_23,
        Eigen::SparseMatrix<double>& mat_33);

    // 行列結合関数（2x2結合）
    static void mergeMatrixWithResize(
        const Eigen::SparseMatrix<double>& free_free,
        const Eigen::SparseMatrix<double>& free_fixed,
        const Eigen::SparseMatrix<double>& fixed_fixed,
        Eigen::SparseMatrix<double>& A);

    // 行列連結関数（縦方向）
    static Eigen::SparseMatrix<double> vstack(
        const Eigen::SparseMatrix<double>& A,
        const Eigen::SparseMatrix<double>& B);

    // 行列連結関数（横方向）
    static Eigen::SparseMatrix<double> hstack(
        const Eigen::SparseMatrix<double>& A,
        const Eigen::SparseMatrix<double>& B);

private:
    SparseMatrixUtils() = delete;  // インスタンス化禁止
};

#endif
