#include "FEVibrationAnalysis.h"
#include "ReducedSystem.h"

#include <Spectra/SymEigsSolver.h>

std::vector<double> FEVibrationAnalysis::ParticipationFactors()
{
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();

    std::vector<double> participation_factors;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        double beta_i = Eigen::VectorXd::Ones(v.size()).dot(vM) / v.dot(vM);

        participation_factors.push_back(beta_i);
    }

    return participation_factors;
}

std::vector<double> FEVibrationAnalysis::ParticipationFactors(Vector direction)
{
    // MassMatrixの組み立て
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();

    std::vector<double> participation_factors;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();

            f(p) = direction.x;
            f(p + 1) = direction.y;
            f(p + 2) = direction.z;
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        double beta_i = f.dot(vM) / v.dot(vM);

        participation_factors.push_back(beta_i);
    }

    return participation_factors;
}

std::vector<Displacement> FEVibrationAnalysis::ParticipationDirectedFactors()
{
    // MassMatrixの組み立て
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();
    std::vector<Displacement> participation_factors;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        vM /= v.dot(vM);
        Displacement facs;
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            // 方向別足し合わせ
            Displacement dj(
                vM(p), vM(p + 1), vM(p + 2),
                vM(p + 3), vM(p + 4), vM(p + 5));
            facs += dj;
        }
        participation_factors.push_back(facs);
    }

    return participation_factors;
}

std::vector<double> FEVibrationAnalysis::EffectiveMassRates()
{
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();

    double total_mass = mass_mat.diagonal().sum();

    std::vector<double> mass_rates;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        double mi = v.dot(vM);
        double beta_i = Eigen::VectorXd::Ones(v.size()).dot(vM) / mi;

        mass_rates.push_back(beta_i * beta_i * mi / total_mass);
    }

    return mass_rates;
}

std::vector<Displacement> FEVibrationAnalysis::EffectiveDirectedMassRates()
{
    // MassMatrixの組み立て
    Eigen::SparseMatrix<double> mass_mat = model->AssembleMassMatrix();
    double total_mass = mass_mat.diagonal().sum();
    std::vector<Displacement> mass_rates;

    for (size_t i = 0; i < ModeNum(); i++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(model->NodeNum() * 6);
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            v(p) = mode_vectors[i][j].Dx();
            v(p + 1) = mode_vectors[i][j].Dy();
            v(p + 2) = mode_vectors[i][j].Dz();
            v(p + 3) = mode_vectors[i][j].Rx();
            v(p + 4) = mode_vectors[i][j].Ry();
            v(p + 5) = mode_vectors[i][j].Rz();
        }

        Eigen::VectorXd vM = mass_mat.selfadjointView<Eigen::Upper>() * v;
        double vMv = v.dot(vM);
        Displacement facs;
        for (size_t j = 0; j < model->NodeNum(); j++)
        {
            int p = j * 6;
            // 方向別足し合わせ
            Displacement dj(
                vM(p), vM(p + 1), vM(p + 2),
                vM(p + 3), vM(p + 4), vM(p + 5));
            facs += dj;
        }

        facs = Displacement(
            facs.Dx() * facs.Dx() / vMv,
            facs.Dy() * facs.Dy() / vMv,
            facs.Dz() * facs.Dz() / vMv,
            facs.Rx() * facs.Rx() / vMv,
            facs.Ry() * facs.Ry() / vMv,
            facs.Rz() * facs.Rz() / vMv);

        facs = facs / (total_mass / 3.0);
        mass_rates.push_back(facs);
    }

    return mass_rates;
}

namespace {

// 振動固有値問題 K φ = ω^2 M φ を、質量あり自由度上の標準固有値問題
// A y = μ y, A = L^T P^T K^{-1} P L, μ = 1/ω^2 として解くための Spectra 用作用素。
// P は質量あり自由度の選択、L は質量ブロック M_mm = L L^T のコレスキー因子
// （元の質量行列は対角だが、剛体リンク縮約後のマスタ自由度には連成が入るため
// 対角とは限らない）。シューア補元（静的縮約）を陽に作らないため K の疎性が
// 保たれ、K の分解は1回で済む。
class VibrationInverseOp {
public:
    using Scalar = double;

    VibrationInverseOp(ISparseSolver& solver, const std::vector<int>& massed_indices,
        const Eigen::SparseMatrix<double>& mass_cholL, int full_size)
        : solver_(solver), massed_indices_(massed_indices),
        mass_cholL_(mass_cholL), full_size_(full_size) {}

    Eigen::Index rows() const { return (Eigen::Index)massed_indices_.size(); }
    Eigen::Index cols() const { return (Eigen::Index)massed_indices_.size(); }

    // y_out = L^T P^T K^{-1} P L x_in
    void perform_op(const double* x_in, double* y_out) const
    {
        int n = (int)massed_indices_.size();
        Eigen::VectorXd v = mass_cholL_ * Eigen::Map<const Eigen::VectorXd>(x_in, n);

        Eigen::VectorXd full_rhs = Eigen::VectorXd::Zero(full_size_);
        for (int i = 0; i < n; i++)
            full_rhs(massed_indices_[i]) = v(i);

        Eigen::VectorXd sol = solver_.solve(full_rhs);

        Eigen::VectorXd g(n);
        for (int i = 0; i < n; i++)
            g(i) = sol(massed_indices_[i]);

        Eigen::Map<Eigen::VectorXd>(y_out, n) = mass_cholL_.transpose() * g;
    }

private:
    ISparseSolver& solver_;
    const std::vector<int>& massed_indices_;
    const Eigen::SparseMatrix<double>& mass_cholL_;
    int full_size_;
};

} // namespace


// 固有値解析の実装。旧 FEModel::SolveVibration をOperator側に移管したもの。
int FEVibrationAnalysis::Compute(const int nev)
{
    int computed_num = nev;

    eigs.clear();
    mode_vectors.clear();
    m_computed = false;

    // 縮約系の構築（RigidLinkを考慮）
    ReducedSystem rs(*model);

    Eigen::SparseMatrix<double> ka, ma;
    rs.Reduce(model->AssembleStiffnessMatrix(), ka);
    rs.Reduce(model->AssembleMassMatrix(), ma);

    // 質量あり自由度の抽出
    // 質量ゼロ自由度は行列を縮約せず、VibrationInverseOp が暗黙に静的縮約と
    // 同じ固有値問題を解く（対角質量では対角ゼロ⇔行・列全体ゼロなので厳密）
    std::vector<int> massed_indices, massless_indices;
    Eigen::VectorXd mdiag = ma.diagonal();
    for (int i = 0; i < mdiag.size(); i++)
    {
        if (mdiag(i) >= 0.0000001)
            massed_indices.push_back(i);
        else
            massless_indices.push_back(i);
    }

    int mat_size = (int)massed_indices.size();
    if (computed_num > mat_size - 1)
        computed_num = mat_size - 1;
    if (computed_num < 1 || mat_size - 2 < computed_num)
        return -1;

    // 質量あり自由度ブロック M_mm のコレスキー分解 M_mm = L L^T
    // （剛体リンクのマスタ自由度は質量が連成するため対角とは限らない。
    //   並べ替えを行わない NaturalOrdering で置換の扱いを不要にする）
    Eigen::SparseMatrix<double> m_sh;
    SparseMatrixUtils::splitMatrixWithResize(ma, massless_indices, m_sh);
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper,
        Eigen::NaturalOrdering<int>> mass_llt(m_sh);
    if (mass_llt.info() != Eigen::Success)
        return -1;
    Eigen::SparseMatrix<double> mass_cholL = mass_llt.matrixL();

    // 剛性が全く付かない自由度（トラス節点の回転、板のドリリング等）で
    // K が特異になるのを防ぐ。PSD の組立行列では対角ゼロ⇔行・列全体ゼロ
    // （完全非連成）なので、質量ゼロの死自由度の対角に正値を置いても
    // 他自由度の解は変わらず、当該モード成分は 0 になる。
    {
        Eigen::VectorXd kdiag = ka.diagonal();
        std::vector<Eigen::Triplet<double>> reg;
        for (int i = 0; i < kdiag.size(); i++)
        {
            if (kdiag(i) > 0.0)
                continue;
            if (mdiag(i) >= 0.0000001)
                return -1; // 質量があるのに剛性ゼロの自由度は解けない
            reg.emplace_back(i, i, 1.0);
        }
        if (!reg.empty()) {
            Eigen::SparseMatrix<double> kreg(ka.rows(), ka.cols());
            kreg.setFromTriplets(reg.begin(), reg.end());
            ka += kreg;
        }
    }

    // K の分解は全体でこの1回のみ
    auto solver_vib = createSolver();
    if (!solver_vib->compute(ka))
        return -1;

    int ncv = 2 * computed_num + 1; // Recommended value
    if (ncv > mat_size) ncv = mat_size;

    // A = L^T P^T K^{-1} P L の最大固有値 μ = 1/ω^2 を求める（ω 昇順で得られる）
    VibrationInverseOp op(*solver_vib, massed_indices, mass_cholL, (int)ka.rows());
    Spectra::SymEigsSolver<VibrationInverseOp> geigs(op, computed_num, ncv);

    geigs.init();
    int nconv = geigs.compute();

    if (geigs.info() != Spectra::CompInfo::Successful)
        return -1;

    Eigen::VectorXd mu = geigs.eigenvalues();   // μ 降順 = ω 昇順
    Eigen::MatrixXd u1s = geigs.eigenvectors(); // 正規直交 (y^T y = 1)

    // 固有値を元の固有値問題に戻す（ω = 角振動数）
    for (int i = 0; i < nconv; i++)
        eigs.push_back(1.0 / sqrt(mu(i)));

    // 縮小空間（master + free）の固有ベクトルを復元: φ = K^{-1} P L y / μ
    // 質量ゼロ自由度の成分も静的縮約関係を満たす形で同時に得られる。
    // φ^T M φ = y^T y = 1 となり質量正規化が構成上厳密に成立する。
    Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(ka.rows(), nconv);
    for (int j = 0; j < nconv; j++)
    {
        Eigen::VectorXd v = mass_cholL * u1s.col(j);
        for (int i = 0; i < mat_size; i++)
            rhs(massed_indices[i], j) = v(i);
    }

    Eigen::MatrixXd reduced_vectors = solver_vib->solveMulti(rhs);
    for (int j = 0; j < nconv; j++)
        reduced_vectors.col(j) /= mu(j);

    // 全体DOFへの固有ベクトルを構築(master DOFはslave DOFへ展開)
    Eigen::MatrixXd eigs_vector = Eigen::MatrixXd::Zero(model->DOFNum(), nconv);
    for (int i = 0; i < nconv; i++)
        eigs_vector.col(i) = rs.ExpandVector(reduced_vectors.col(i), model->DOFNum());

    // 固有ベクトルを格納（φ^T M φ = 1 の質量正規化済み）
    for (int i = 0; i < nconv; i++)
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

    m_computed = true;
    return nconv;
}

