#include "CudaSolverBase.h"
#include <map>
#include <algorithm>

CSRMatrix eigenUpperToFullCSR(const Eigen::SparseMatrix<double>& A) {
    int n = A.rows();

    // Collect all entries from upper triangle, then mirror to lower
    // Use a map to accumulate: (row, col) -> value
    std::vector<std::vector<std::pair<int, double>>> rowEntries(n);

    // Eigen SparseMatrix is CSC: iterate by columns
    for (int col = 0; col < A.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            int row = it.row();
            double val = it.value();

            // Add the entry itself
            rowEntries[row].emplace_back(col, val);

            // Add the mirror entry (if off-diagonal)
            if (row != col) {
                rowEntries[col].emplace_back(row, val);
            }
        }
    }

    // Sort each row by column index
    for (int i = 0; i < n; ++i) {
        std::sort(rowEntries[i].begin(), rowEntries[i].end());
    }

    // Build CSR
    CSRMatrix csr;
    csr.n = n;
    csr.rowPtr.resize(n + 1, 0);

    for (int i = 0; i < n; ++i) {
        csr.rowPtr[i + 1] = csr.rowPtr[i] + static_cast<int>(rowEntries[i].size());
    }
    csr.nnz = csr.rowPtr[n];
    csr.colIdx.resize(csr.nnz);
    csr.values.resize(csr.nnz);

    for (int i = 0; i < n; ++i) {
        int offset = csr.rowPtr[i];
        for (size_t j = 0; j < rowEntries[i].size(); ++j) {
            csr.colIdx[offset + j] = rowEntries[i][j].first;
            csr.values[offset + j] = rowEntries[i][j].second;
        }
    }

    return csr;
}

CSRMatrix eigenUpperToLowerCSR(const Eigen::SparseMatrix<double>& A) {
    // For a CSC matrix stored as upper triangular:
    // CSC stores: outerIndexPtr (column pointers), innerIndexPtr (row indices), values
    // Reinterpreting these as CSR gives us:
    //   rowPtr = outerIndexPtr (was column pointers, now row pointers)
    //   colIdx = innerIndexPtr (was row indices, now column indices)
    //   values = values
    // This effectively transposes the matrix: CSR(A^T) = CSC(A)
    // Since A is upper triangular and symmetric, A^T is lower triangular
    // So this gives us CSR of the lower triangular part

    int n = A.cols();
    int nnz = A.nonZeros();

    CSRMatrix csr;
    csr.n = n;
    csr.nnz = nnz;
    csr.rowPtr.resize(n + 1);
    csr.colIdx.resize(nnz);
    csr.values.resize(nnz);

    const int* outerPtr = A.outerIndexPtr();
    const int* innerPtr = A.innerIndexPtr();
    const double* valPtr = A.valuePtr();

    // CSC outer = column pointers -> CSR row pointers (transposed)
    for (int i = 0; i <= n; ++i) {
        csr.rowPtr[i] = outerPtr[i];
    }
    for (int i = 0; i < nnz; ++i) {
        csr.colIdx[i] = innerPtr[i];
        csr.values[i] = valPtr[i];
    }

    return csr;
}
