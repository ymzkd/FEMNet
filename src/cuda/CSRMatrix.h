#ifndef _CSR_MATRIX_H_
#define _CSR_MATRIX_H_

#include <Eigen/Sparse>
#include <vector>

// CSR matrix on host (no CUDA dependency)
struct CSRMatrix {
    int n = 0;
    int nnz = 0;
    std::vector<int> rowPtr;
    std::vector<int> colIdx;
    std::vector<double> values;
};

// Convert Eigen CSC upper-triangular to full symmetric CSR
CSRMatrix eigenUpperToFullCSR(const Eigen::SparseMatrix<double>& A);

// Convert Eigen CSC upper-triangular to lower-triangular CSR
CSRMatrix eigenUpperToLowerCSR(const Eigen::SparseMatrix<double>& A);

#endif
