#include <Eigen/Sparse>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/SymGEigsSolver.h>

#include "FEBucklingAnalysis.h"
#include "PardisoTriOp.h"

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

    // 剛性が全く付かない自由度（トラス節点の回転等）で K が特異になるのを防ぐ。
    // PSD の組立行列では対角ゼロ⇔行・列全体ゼロ（完全非連成）なので、幾何剛性も
    // 持たない死自由度なら対角に正値を置いても他自由度の解・固有値は変わらない。
    {
        Eigen::VectorXd kdiag = ka.diagonal();
        std::vector<bool> kg_active(kg.rows(), false);
        for (int c = 0; c < kg.outerSize(); c++)
            for (Eigen::SparseMatrix<double>::InnerIterator it(kg, c); it; ++it)
                if (it.value() != 0.0) {
                    kg_active[it.row()] = true;
                    kg_active[it.col()] = true;
                }
        std::vector<Eigen::Triplet<double>> reg;
        for (int i = 0; i < kdiag.size(); i++)
        {
            if (kdiag(i) > 0.0)
                continue;
            if (kg_active[i])
                return -1; // 幾何剛性があるのに弾性剛性ゼロの自由度は解けない
            reg.emplace_back(i, i, 1.0);
        }
        if (!reg.empty()) {
            Eigen::SparseMatrix<double> kreg(ka.rows(), ka.cols());
            kreg.setFromTriplets(reg.begin(), reg.end());
            ka += kreg;
        }
    }

    int mat_size = (int)ka.rows();
    if (computed_num > mat_size - 1)
        computed_num = mat_size - 1;
    if (computed_num < 1 || mat_size - 2 < computed_num)
        return -1;

    int ncv = 2 * computed_num + 1; // Recommended value
    if (ncv > mat_size) ncv = mat_size;

    Eigen::SparseMatrix<double> neg_kg = -kg;
    using OpType = Spectra::SparseSymMatProd<double, Eigen::Upper>;
    OpType A_op(neg_kg);

    int nconv = 0;
    Eigen::MatrixXd part_eigen_vectors;
    std::vector<double> part_eigen_values;

    try {
#ifdef EIGEN_USE_MKL_ALL
        // Cholesky モード + Pardiso 部分求解 (phase 331/333)。
        // K の分解は並列で高速、反復は三角求解2回のみで、RegularInverse モードで
        // 必要だった直交化の B 内積（K·v の SpMV）も不要。
        PardisoTriOp B_op(ka);
        if (B_op.info() != Spectra::CompInfo::Successful)
            return -1;
        Spectra::SymGEigsSolver<OpType, PardisoTriOp, Spectra::GEigsMode::Cholesky>
            geigs(A_op, B_op, computed_num, ncv);
#else
        // MKL なし: 従来通り SimplicialLLT ベースの Cholesky モード
        Spectra::SparseCholesky<double, Eigen::Upper> B_op(ka);
        Spectra::SymGEigsSolver<OpType, Spectra::SparseCholesky<double, Eigen::Upper>,
            Spectra::GEigsMode::Cholesky> geigs(A_op, B_op, computed_num, ncv);
#endif

        geigs.init();
        nconv = geigs.compute(Spectra::SortRule::LargestAlge);

        if (geigs.info() != Spectra::CompInfo::Successful)
            return -1;

        part_eigen_vectors = geigs.eigenvectors();
        for (double v : geigs.eigenvalues())
            part_eigen_values.push_back(v);
    }
    catch (const std::exception&) {
        // 分解・求解の失敗（特異行列等）
        return -1;
    }

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
    for (double v : part_eigen_values)
        eigs.push_back(1.0 / v);

    mode_num = nconv;
    return mode_num;
}
