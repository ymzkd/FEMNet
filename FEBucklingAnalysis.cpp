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

    // インデックスの取得（RigidLinkを考慮）
    std::vector<int> slave_indices = model->RigidLinkData->SlaveDOFIndices();
    std::vector<int> free_indices = model->FreeIndices(true);  // rigid_link=true
    std::vector<int> fixed_indices = model->FixIndices();

    // 剛性行列の組み立て
    Eigen::SparseMatrix<double> k_full = model->AssembleStiffnessMatrix();
    if (InitailDeformOp != nullptr)
        k_full += model->AssembleGeometricStiffnessMatrix(InitailDeformOp->GetDisplacements());

    // 幾何剛性行列の組み立て
    Eigen::SparseMatrix<double> kg_full =
        model->AssembleGeometricStiffnessMatrix(deform_case->GetDisplacements());

    // 変換行列の取得
    Eigen::SparseMatrix<double> linkTransMat =
        model->RigidLinkData->TransformationMatrix().sparseView(1e-10);
    int master_dof_num = linkTransMat.cols();

    Eigen::SparseMatrix<double> ka, kg;

    if (master_dof_num > 0) {

        // RigidLinkがある場合: 3x3ブロックに分割して縮小
        Eigen::SparseMatrix<double> k11, k12, k13, k22, k23, k33;
        SparseMatrixUtils::splitMatrix3x3(k_full, slave_indices, free_indices,
            k11, k12, k13, k22, k23, k33);

        Eigen::SparseMatrix<double> g11, g12, g13, g22, g23, g33;
        SparseMatrixUtils::splitMatrix3x3(kg_full, slave_indices, free_indices,
            g11, g12, g13, g22, g23, g33);

        // 剛性行列の縮小
        Eigen::SparseMatrix<double> kaa, kab;
        kaa = (linkTransMat.transpose() * k11.selfadjointView<Eigen::Upper>() * linkTransMat)
              .triangularView<Eigen::Upper>();
        kab = (linkTransMat.transpose() * k12);
        SparseMatrixUtils::mergeMatrixWithResize(kaa, kab, k22, ka);

        // 幾何剛性行列の縮小
        Eigen::SparseMatrix<double> gaa, gab;
        gaa = (linkTransMat.transpose() * g11.selfadjointView<Eigen::Upper>() * linkTransMat)
              .triangularView<Eigen::Upper>();
        gab = (linkTransMat.transpose() * g12);
        SparseMatrixUtils::mergeMatrixWithResize(gaa, gab, g22, kg);
    }
    else {
        // RigidLinkがない場合: 従来通り2x2分割
        SparseMatrixUtils::splitMatrixWithResize(k_full, fixed_indices, ka);
        SparseMatrixUtils::splitMatrixWithResize(kg_full, fixed_indices, kg);
    }

    int ncv = 2 * computed_num + 1; // Recommended value

    using OpType = Spectra::SparseSymMatProd<double, Eigen::Upper>;
    using BOpType = Spectra::SparseCholesky<double, Eigen::Upper>;
    OpType A_op(-kg); // Invert
    BOpType B_op(ka); // Invert
    Spectra::SymGEigsSolver<OpType, BOpType, Spectra::GEigsMode::Cholesky>
        geigs(A_op, B_op, computed_num, ncv);

    geigs.init();
    int nconv = geigs.compute(Spectra::SortRule::LargestAlge);

    if (geigs.info() == Spectra::CompInfo::Successful)
    {
        Eigen::MatrixXd part_eigen_vectors = geigs.eigenvectors();
        Eigen::MatrixXd eigs_vector = Eigen::MatrixXd::Zero(model->DOFNum(), computed_num);

        if (master_dof_num > 0) {
            // RigidLinkがある場合: master DOFをslave DOFに展開
            for (size_t i = 0; i < nconv; i++) {
                Eigen::VectorXd part_vec = part_eigen_vectors.col(i);
                Eigen::VectorXd d_master = part_vec.head(master_dof_num);
                Eigen::VectorXd d_free = part_vec.tail(free_indices.size());
                Eigen::VectorXd d_slave = linkTransMat * d_master;

                for (size_t j = 0; j < slave_indices.size(); j++)
                    eigs_vector(slave_indices[j], i) = d_slave(j);
                for (size_t j = 0; j < free_indices.size(); j++)
                    eigs_vector(free_indices[j], i) = d_free(j);
            }
        }
        else {
            // RigidLinkがない場合: 従来通り
            for (size_t i = 0; i < free_indices.size(); i++)
                eigs_vector.row(free_indices[i]) = part_eigen_vectors.row(i);
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
    return mode_num;
}
