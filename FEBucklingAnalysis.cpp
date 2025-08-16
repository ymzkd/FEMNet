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

#include "FEBucklingAnalysis.h"

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

        // using OpType = Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse, Eigen::Upper>;
        // using BOpType = Spectra::SparseSymMatProd<double, Eigen::Upper>;
        // OpType A_op(ka, kg);
        // BOpType B_op(kg);
        // Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::Buckling>
        //    geigs(A_op, B_op, computed_num, ncv, 0.1);

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
        // int nconv = geigs.compute(Spectra::SortRule::LargestMagn);
        // int nconv = geigs.compute(Spectra::SortRule::SmallestMagn);
        //  int nconv = geigs.compute(Spectra::SortRule::SmallestAlge);
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
    else if (solver_selection == 2)
    {
        // Dense Solver
        Eigen::MatrixXd kg_dense = Eigen::MatrixXd(kg);
        Eigen::MatrixXd k_dense = Eigen::MatrixXd(ka);

        // 行列演算子を作成
        Spectra::DenseSymMatProd<double, Eigen::Upper> kg_dense_op(-kg_dense);
        Spectra::DenseCholesky<double, Eigen::Upper> k_dense_op(k_dense);

        Spectra::SymGEigsSolver<Spectra::DenseSymMatProd<double, Eigen::Upper>,
                                Spectra::DenseCholesky<double, Eigen::Upper>,
                                Spectra::GEigsMode::Cholesky>
            solver(kg_dense_op, k_dense_op, computed_num, 2 * computed_num + 1);

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
    else
    {
        // 密行列への変換
        Eigen::MatrixXd kg_dense = Eigen::MatrixXd(-kg);
        Eigen::MatrixXd k_dense = Eigen::MatrixXd(ka);

        // 対称行列の場合、下三角部分を補完
        kg_dense.triangularView<Eigen::Lower>() = kg_dense.triangularView<Eigen::Upper>().transpose();
        k_dense.triangularView<Eigen::Lower>() = k_dense.triangularView<Eigen::Upper>().transpose();

        // 一般固有値問題を解く
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(kg_dense, k_dense);

        if (solver.info() == Eigen::Success)
        {
            Eigen::VectorXd eigs_arr = solver.eigenvalues();
            for (size_t i = 0; i < mode_num; i++)
            {
                double eig = eigs_arr[eigs_arr.size() - i - 1];
                eigs.push_back(1.0 / eig);
            }

            // for (size_t i = 0; i < mode_num; i++)
            //{
            //     double eig = eigs_arr[eigs_arr.size() - i - 1];
            //     eigs.push_back(1.0 / eig);
            // }

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
        else
        {
            return -1;
        }
    }

    return mode_num;
}
