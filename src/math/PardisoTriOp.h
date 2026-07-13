#ifndef _PARDISO_TRI_OP_H_
#define _PARDISO_TRI_OP_H_

#ifdef EIGEN_USE_MKL_ALL

#include <stdexcept>
#include <string>
#include <Eigen/Sparse>
#include <Spectra/Util/CompInfo.h>
#include <mkl.h>

/// <summary>
/// Spectra の Cholesky モード（SymGEigsSolver&lt;..., GEigsMode::Cholesky&gt;）互換の
/// B 作用素を MKL pardiso の部分求解 (phase 331/333) で実装するクラス。
///
/// mtype=2 (SPD) の分解 B = P^T L L^T P に対し、各 phase は元の並びのベクトルへ
/// 置換込みで作用するため、331/333 のペアは B = R^T R の R^{-T}・R^{-1} として
/// 一貫に振る舞う（R は三角行列である必要はなく、B の平方根因子であれば
/// Lanczos の標準形変換 C = R^{-T} A R^{-1} は成立する）。
///
/// Spectra 内蔵の SparseCholesky（SimplicialLLT・シングルスレッド）に対して
/// 分解が並列で大幅に速く、RegularInverse モードと異なり直交化の B 内積
/// （K·v の SpMV）も不要なため、座屈固有値解析の反復コストを最小にできる。
/// </summary>
class PardisoTriOp {
public:
    using Scalar = double;

    /// B は上三角格納（または全体格納）の対称正定値疎行列
    explicit PardisoTriOp(const Eigen::SparseMatrix<double>& B) : a_(B) {
        a_.makeCompressed();
        n_ = (MKL_INT)B.rows();
        MKL_INT mtype = 2;
        pardisoinit(pt_, &mtype, iparm_);
        iparm_[34] = 1;  // 0-based indexing
        iparm_[7] = 0;   // 反復改良なし（部分求解では元々行われない）
        MKL_INT maxfct = 1, mnum = 1, phase = 12, nrhs = 1, msglvl = 0, error = 0;
        MKL_INT idum = 0;
        double ddum = 0;
        pardiso(pt_, &maxfct, &mnum, &mtype, &phase, &n_,
                (void*)a_.valuePtr(), (MKL_INT*)a_.outerIndexPtr(), (MKL_INT*)a_.innerIndexPtr(),
                &idum, &nrhs, iparm_, &msglvl, &ddum, &ddum, &error);
        info_ = (error == 0) ? Spectra::CompInfo::Successful
                             : Spectra::CompInfo::NumericalIssue;
    }

    ~PardisoTriOp() {
        MKL_INT mtype = 2, maxfct = 1, mnum = 1, phase = -1, nrhs = 1, msglvl = 0, error = 0;
        MKL_INT idum = 0;
        double ddum = 0;
        pardiso(pt_, &maxfct, &mnum, &mtype, &phase, &n_, &ddum, &idum, &idum,
                &idum, &nrhs, iparm_, &msglvl, &ddum, &ddum, &error);
    }

    PardisoTriOp(const PardisoTriOp&) = delete;
    PardisoTriOp& operator=(const PardisoTriOp&) = delete;

    Eigen::Index rows() const { return n_; }
    Eigen::Index cols() const { return n_; }
    Spectra::CompInfo info() const { return info_; }

    /// y = R^{-T} x（前進代入, phase 331）
    void lower_triangular_solve(const double* x_in, double* y_out) const {
        run_phase(331, x_in, y_out);
    }
    /// y = R^{-1} x（後退代入, phase 333）
    void upper_triangular_solve(const double* x_in, double* y_out) const {
        run_phase(333, x_in, y_out);
    }

private:
    void run_phase(int ph, const double* x_in, double* y_out) const {
        MKL_INT mtype = 2, maxfct = 1, mnum = 1, phase = (MKL_INT)ph, nrhs = 1, msglvl = 0, error = 0;
        MKL_INT idum = 0;
        pardiso(pt_, &maxfct, &mnum, &mtype, &phase, &n_,
                (void*)a_.valuePtr(), (MKL_INT*)a_.outerIndexPtr(), (MKL_INT*)a_.innerIndexPtr(),
                &idum, &nrhs, iparm_, &msglvl,
                const_cast<double*>(x_in), y_out, &error);
        if (error != 0)
            throw std::runtime_error("PardisoTriOp: partial solve failed, error=" +
                                     std::to_string((long long)error));
    }

    Eigen::SparseMatrix<double, Eigen::RowMajor> a_;  // 上三角 CSR
    MKL_INT n_ = 0;
    Spectra::CompInfo info_ = Spectra::CompInfo::NotComputed;
    mutable void* pt_[64] = {};
    mutable MKL_INT iparm_[64] = {};
};

#endif // EIGEN_USE_MKL_ALL

#endif // _PARDISO_TRI_OP_H_
