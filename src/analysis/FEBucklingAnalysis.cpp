#include <Eigen/Sparse>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/SymGEigsSolver.h>

#include "FEBucklingAnalysis.h"
#include "PardisoTriOp.h"
#include "ReducedSystem.h"

int FEBucklingAnalysis::SolveBuckling()
{
    int computed_num = mode_num;

    // 縮約系の構築（RigidLinkを考慮）
    ReducedSystem rs(*model);

    // 剛性行列の組み立て
    Eigen::SparseMatrix<double> k_full = model->AssembleStiffnessMatrix();
    if (InitailDeformOp != nullptr)
        k_full += model->AssembleGeometricStiffnessMatrix(InitailDeformOp->GetDisplacements());

    // 幾何剛性行列の組み立て
    Eigen::SparseMatrix<double> kg_full =
        model->AssembleGeometricStiffnessMatrix(deform_case->GetDisplacements());

    Eigen::SparseMatrix<double> ka, kg;
    rs.Reduce(k_full, ka);
    rs.Reduce(kg_full, kg);

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

    // 全体DOFへの固有ベクトルを構築(master DOFはslave DOFへ展開)
    Eigen::MatrixXd eigs_vector = Eigen::MatrixXd::Zero(model->DOFNum(), computed_num);
    for (int i = 0; i < nconv; i++)
        eigs_vector.col(i) = rs.ExpandVector(part_eigen_vectors.col(i), model->DOFNum());

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
