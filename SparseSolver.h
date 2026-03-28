#ifndef _SPARSE_SOLVER_H_
#define _SPARSE_SOLVER_H_

#include "SparseSolverInterface.h"
#include <iostream>

#ifdef EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif

class EigenSolver : public ISparseSolver {
private:
#ifdef EIGEN_USE_MKL_ALL
    Eigen::PardisoLLT<Eigen::SparseMatrix<double>> solver_;
#else
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper> solver_;
#endif
    bool success_ = false;

public:
    bool compute(const Eigen::SparseMatrix<double>& A) override {
        solver_.compute(A);
        success_ = (solver_.info() == Eigen::Success);
        if (!success_) {
            std::cerr << "EigenSolver::compute() failed" << std::endl;
        }
        return success_;
    }

    Eigen::VectorXd solve(const Eigen::VectorXd& b) override {
        return solver_.solve(b);
    }

    Eigen::MatrixXd solveMulti(const Eigen::MatrixXd& B) override {
        return solver_.solve(B);
    }

    bool success() const override { return success_; }

    std::string name() const override {
#ifdef EIGEN_USE_MKL_ALL
        return "MKL PardisoLLT";
#else
        return "Eigen SimplicialLLT";
#endif
    }
};

#ifdef USE_CUDA
#include "cuda/CudaPCGSolver.h"
#endif

inline std::unique_ptr<ISparseSolver> createSolver(SolverBackend backend =
#ifdef USE_CUDA
    SolverBackend::CudaPCG
#else
    SolverBackend::Default
#endif
) {
#ifdef USE_CUDA
    if (backend == SolverBackend::CudaPCG) {
        return std::make_unique<CudaPCGSolver>();
    }
#endif
    (void)backend;
    return std::make_unique<EigenSolver>();
}

#endif // _SPARSE_SOLVER_H_
