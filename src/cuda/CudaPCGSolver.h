#ifndef _CUDA_PCG_SOLVER_H_
#define _CUDA_PCG_SOLVER_H_

#include "SparseSolverInterface.h"

struct CudaPCGSolverImpl;

class CudaPCGSolver : public ISparseSolver {
private:
    std::unique_ptr<CudaPCGSolverImpl> impl_;
    bool success_ = false;
    double tol_;
    int maxIter_;

public:
    CudaPCGSolver(double tol = 1e-10, int maxIter = 10000);
    ~CudaPCGSolver() override;

    bool compute(const Eigen::SparseMatrix<double>& A) override;
    Eigen::VectorXd solve(const Eigen::VectorXd& b) override;
    Eigen::MatrixXd solveMulti(const Eigen::MatrixXd& B) override;
    bool success() const override { return success_; }
    std::string name() const override { return "CUDA PCG (IC0)"; }
};

#endif
