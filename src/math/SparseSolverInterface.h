#ifndef _SPARSE_SOLVER_INTERFACE_H_
#define _SPARSE_SOLVER_INTERFACE_H_

#include <Eigen/Sparse>
#include <memory>
#include <string>

enum class SolverBackend {
    Default,     // Eigen SimplicialLLT (or PardisoLLT if MKL)
    CudaPCG      // Preconditioned CG (IC0) on GPU
};

class ISparseSolver {
public:
    virtual ~ISparseSolver() = default;
    virtual bool compute(const Eigen::SparseMatrix<double>& A) = 0;
    virtual Eigen::VectorXd solve(const Eigen::VectorXd& b) = 0;
    virtual Eigen::MatrixXd solveMulti(const Eigen::MatrixXd& B) = 0;
    virtual bool success() const = 0;
    virtual std::string name() const = 0;
};

#endif // _SPARSE_SOLVER_INTERFACE_H_
