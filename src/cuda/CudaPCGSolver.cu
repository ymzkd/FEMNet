#include "CudaPCGSolver.h"
#include "CudaSolverBase.h"
#include <cmath>

struct CudaPCGSolverImpl {
    cusparseHandle_t cusparseHandle = nullptr;
    cublasHandle_t cublasHandle = nullptr;

    // Matrix (full symmetric CSR)
    DeviceBuffer<double> d_values;
    DeviceBuffer<int> d_rowPtr;
    DeviceBuffer<int> d_colIdx;
    cusparseSpMatDescr_t matDescr = nullptr;

    // CG vectors
    DeviceBuffer<double> d_r;
    DeviceBuffer<double> d_p;
    DeviceBuffer<double> d_Ap;
    DeviceBuffer<double> d_x;
    DeviceBuffer<double> d_b;
    DeviceBuffer<char> d_spmvBuffer;

    // IC(0) preconditioner
    DeviceBuffer<double> d_icValues;
    DeviceBuffer<int> d_icRowPtr;
    DeviceBuffer<int> d_icColIdx;
    DeviceBuffer<double> d_z;
    DeviceBuffer<double> d_tmp;

    // IC(0) factorization (legacy API - still available in CUDA 12.8)
    csric02Info_t ic02Info = nullptr;
    cusparseMatDescr_t icDescr = nullptr;
    DeviceBuffer<char> d_ic02Buffer;

    // Generic SpSV for triangular solves (new API)
    cusparseSpMatDescr_t spMatL = nullptr;   // lower triangular L
    cusparseSpSVDescr_t spsvDescrL = nullptr;
    cusparseSpSVDescr_t spsvDescrLT = nullptr;
    DeviceBuffer<char> d_svBufferL;
    DeviceBuffer<char> d_svBufferLT;

    // Dense vector descriptors for SpSV (preconditioner)
    cusparseDnVecDescr_t dnVecTmp = nullptr;
    cusparseDnVecDescr_t dnVecZ = nullptr;
    cusparseDnVecDescr_t dnVecR_sv = nullptr;

    // Dense vector descriptors for SpMV (CG loop) - created once, reused
    cusparseDnVecDescr_t dnVecP = nullptr;
    cusparseDnVecDescr_t dnVecAp = nullptr;

    int n = 0;
    int nnz = 0;
    bool computed = false;
    bool precondReady = false;

    CudaPCGSolverImpl() {
        cusparseCreate(&cusparseHandle);
        cublasCreate(&cublasHandle);
    }

    ~CudaPCGSolverImpl() {
        cleanup();
        if (cusparseHandle) cusparseDestroy(cusparseHandle);
        if (cublasHandle) cublasDestroy(cublasHandle);
    }

    void cleanup() {
        if (matDescr) { cusparseDestroySpMat(matDescr); matDescr = nullptr; }
        if (ic02Info) { cusparseDestroyCsric02Info(ic02Info); ic02Info = nullptr; }
        if (icDescr) { cusparseDestroyMatDescr(icDescr); icDescr = nullptr; }
        if (spMatL) { cusparseDestroySpMat(spMatL); spMatL = nullptr; }
        if (spsvDescrL) { cusparseSpSV_destroyDescr(spsvDescrL); spsvDescrL = nullptr; }
        if (spsvDescrLT) { cusparseSpSV_destroyDescr(spsvDescrLT); spsvDescrLT = nullptr; }
        if (dnVecP) { cusparseDestroyDnVec(dnVecP); dnVecP = nullptr; }
        if (dnVecAp) { cusparseDestroyDnVec(dnVecAp); dnVecAp = nullptr; }
        if (dnVecTmp) { cusparseDestroyDnVec(dnVecTmp); dnVecTmp = nullptr; }
        if (dnVecZ) { cusparseDestroyDnVec(dnVecZ); dnVecZ = nullptr; }
        if (dnVecR_sv) { cusparseDestroyDnVec(dnVecR_sv); dnVecR_sv = nullptr; }

        d_values.free(); d_rowPtr.free(); d_colIdx.free();
        d_r.free(); d_p.free(); d_Ap.free(); d_x.free(); d_b.free();
        d_spmvBuffer.free();
        d_icValues.free(); d_icRowPtr.free(); d_icColIdx.free();
        d_z.free(); d_tmp.free();
        d_ic02Buffer.free(); d_svBufferL.free(); d_svBufferLT.free();

        computed = false;
        precondReady = false;
    }

    bool setupSpMV() {
        cusparseStatus_t st = cusparseCreateCsr(
            &matDescr, n, n, nnz,
            d_rowPtr.get(), d_colIdx.get(), d_values.get(),
            CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);
        if (st != CUSPARSE_STATUS_SUCCESS) return false;

        cusparseCreateDnVec(&dnVecP, n, d_p.get(), CUDA_R_64F);
        cusparseCreateDnVec(&dnVecAp, n, d_Ap.get(), CUDA_R_64F);

        double alpha = 1.0, beta = 0.0;
        size_t bufferSize = 0;
        st = cusparseSpMV_bufferSize(
            cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            &alpha, matDescr, dnVecP, &beta, dnVecAp,
            CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize);
        if (st != CUSPARSE_STATUS_SUCCESS) return false;

        if (bufferSize > 0) {
            if (!d_spmvBuffer.alloc(bufferSize)) return false;
        }
        return true;
    }

    bool setupPreconditioner() {
        // Copy matrix for IC(0)
        if (!d_icValues.alloc(nnz) || !d_icRowPtr.alloc(n + 1) || !d_icColIdx.alloc(nnz)) return false;

        cudaMemcpy(d_icRowPtr.get(), d_rowPtr.get(), (n + 1) * sizeof(int), cudaMemcpyDeviceToDevice);
        cudaMemcpy(d_icColIdx.get(), d_colIdx.get(), nnz * sizeof(int), cudaMemcpyDeviceToDevice);
        cudaMemcpy(d_icValues.get(), d_values.get(), nnz * sizeof(double), cudaMemcpyDeviceToDevice);

        if (!d_z.alloc(n) || !d_tmp.alloc(n)) return false;

        // General descriptor for IC(0) (legacy API)
        cusparseCreateMatDescr(&icDescr);
        cusparseSetMatType(icDescr, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(icDescr, CUSPARSE_INDEX_BASE_ZERO);

        cusparseCreateCsric02Info(&ic02Info);

        // IC(0) buffer
        int ic02BufferSize = 0;
        cusparseStatus_t st = cusparseDcsric02_bufferSize(
            cusparseHandle, n, nnz, icDescr,
            d_icValues.get(), d_icRowPtr.get(), d_icColIdx.get(),
            ic02Info, &ic02BufferSize);
        if (st != CUSPARSE_STATUS_SUCCESS) return false;
        if (!d_ic02Buffer.alloc(ic02BufferSize)) return false;

        // IC(0) analysis
        st = cusparseDcsric02_analysis(
            cusparseHandle, n, nnz, icDescr,
            d_icValues.get(), d_icRowPtr.get(), d_icColIdx.get(),
            ic02Info, CUSPARSE_SOLVE_POLICY_USE_LEVEL,
            d_ic02Buffer.get());
        if (st != CUSPARSE_STATUS_SUCCESS) return false;

        // IC(0) factorization
        st = cusparseDcsric02(
            cusparseHandle, n, nnz, icDescr,
            d_icValues.get(), d_icRowPtr.get(), d_icColIdx.get(),
            ic02Info, CUSPARSE_SOLVE_POLICY_USE_LEVEL,
            d_ic02Buffer.get());
        if (st != CUSPARSE_STATUS_SUCCESS) return false;

        // Setup generic SpSV for L and L^T triangular solves
        // Create sparse matrix descriptor for L (lower triangular, using IC values)
        st = cusparseCreateCsr(
            &spMatL, n, n, nnz,
            d_icRowPtr.get(), d_icColIdx.get(), d_icValues.get(),
            CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);
        if (st != CUSPARSE_STATUS_SUCCESS) return false;

        // Set lower triangular fill mode
        cusparseFillMode_t fillMode = CUSPARSE_FILL_MODE_LOWER;
        cusparseDiagType_t diagType = CUSPARSE_DIAG_TYPE_NON_UNIT;
        cusparseSpMatSetAttribute(spMatL, CUSPARSE_SPMAT_FILL_MODE,
            &fillMode, sizeof(fillMode));
        cusparseSpMatSetAttribute(spMatL, CUSPARSE_SPMAT_DIAG_TYPE,
            &diagType, sizeof(diagType));

        // Create dense vector descriptors
        cusparseCreateDnVec(&dnVecR_sv, n, d_r.get(), CUDA_R_64F);
        cusparseCreateDnVec(&dnVecTmp, n, d_tmp.get(), CUDA_R_64F);
        cusparseCreateDnVec(&dnVecZ, n, d_z.get(), CUDA_R_64F);

        // Create SpSV descriptors and analyze
        double one = 1.0;

        // L solve: L * tmp = r
        cusparseSpSV_createDescr(&spsvDescrL);
        size_t bufSizeL = 0;
        st = cusparseSpSV_bufferSize(
            cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            &one, spMatL, dnVecR_sv, dnVecTmp,
            CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL, &bufSizeL);
        if (st != CUSPARSE_STATUS_SUCCESS) return false;
        if (!d_svBufferL.alloc(bufSizeL > 0 ? bufSizeL : 1)) return false;

        st = cusparseSpSV_analysis(
            cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            &one, spMatL, dnVecR_sv, dnVecTmp,
            CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL, d_svBufferL.get());
        if (st != CUSPARSE_STATUS_SUCCESS) return false;

        // L^T solve: L^T * z = tmp
        cusparseSpSV_createDescr(&spsvDescrLT);
        size_t bufSizeLT = 0;
        st = cusparseSpSV_bufferSize(
            cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE,
            &one, spMatL, dnVecTmp, dnVecZ,
            CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrLT, &bufSizeLT);
        if (st != CUSPARSE_STATUS_SUCCESS) return false;
        if (!d_svBufferLT.alloc(bufSizeLT > 0 ? bufSizeLT : 1)) return false;

        st = cusparseSpSV_analysis(
            cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE,
            &one, spMatL, dnVecTmp, dnVecZ,
            CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrLT, d_svBufferLT.get());
        if (st != CUSPARSE_STATUS_SUCCESS) return false;

        precondReady = true;
        return true;
    }

    bool applyPreconditioner() {
        // Solve L * tmp = r, then L^T * z = tmp
        double one = 1.0;

        cusparseStatus_t st = cusparseSpSV_solve(
            cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            &one, spMatL, dnVecR_sv, dnVecTmp,
            CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL);
        if (st != CUSPARSE_STATUS_SUCCESS) return false;

        st = cusparseSpSV_solve(
            cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE,
            &one, spMatL, dnVecTmp, dnVecZ,
            CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrLT);
        if (st != CUSPARSE_STATUS_SUCCESS) return false;

        return true;
    }
};

CudaPCGSolver::CudaPCGSolver(double tol, int maxIter)
    : impl_(std::make_unique<CudaPCGSolverImpl>()), tol_(tol), maxIter_(maxIter) {}

CudaPCGSolver::~CudaPCGSolver() = default;

bool CudaPCGSolver::compute(const Eigen::SparseMatrix<double>& A) {
    auto& d = *impl_;
    d.cleanup();
    success_ = false;

    CSRMatrix csr = eigenUpperToFullCSR(A);
    d.n = csr.n;
    d.nnz = csr.nnz;

    if (!d.d_rowPtr.alloc(d.n + 1) || !d.d_colIdx.alloc(d.nnz) || !d.d_values.alloc(d.nnz)) return false;
    if (cudaMemcpy(d.d_rowPtr.get(), csr.rowPtr.data(), (d.n + 1) * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemcpy(d.d_colIdx.get(), csr.colIdx.data(), d.nnz * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemcpy(d.d_values.get(), csr.values.data(), d.nnz * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;

    if (!d.d_r.alloc(d.n) || !d.d_p.alloc(d.n) || !d.d_Ap.alloc(d.n) ||
        !d.d_x.alloc(d.n) || !d.d_b.alloc(d.n)) return false;

    if (!d.setupSpMV()) {
        std::cerr << "CudaPCGSolver: SpMV setup failed" << std::endl;
        return false;
    }
    d.computed = true;

    if (!d.setupPreconditioner()) {
        std::cerr << "CudaPCGSolver: preconditioner setup failed" << std::endl;
        return false;
    }

    success_ = true;
    return true;
}

Eigen::VectorXd CudaPCGSolver::solve(const Eigen::VectorXd& b) {
    auto& d = *impl_;
    Eigen::VectorXd x(d.n);
    x.setZero();

    if (!d.computed || !d.precondReady) {
        std::cerr << "CudaPCGSolver::solve() called without compute()" << std::endl;
        return x;
    }

    cudaMemcpy(d.d_b.get(), b.data(), d.n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(d.d_x.get(), 0, d.n * sizeof(double));
    cudaMemcpy(d.d_r.get(), d.d_b.get(), d.n * sizeof(double), cudaMemcpyDeviceToDevice);

    double bnorm = 0.0;
    cublasDdot(d.cublasHandle, d.n, d.d_b.get(), 1, d.d_b.get(), 1, &bnorm);
    bnorm = std::sqrt(bnorm);
    if (bnorm < 1e-30) return x;
    double tol2 = tol_ * tol_ * bnorm * bnorm;

    // z = M^{-1} r
    d.applyPreconditioner();

    // p = z
    cudaMemcpy(d.d_p.get(), d.d_z.get(), d.n * sizeof(double), cudaMemcpyDeviceToDevice);

    double rz = 0.0;
    cublasDdot(d.cublasHandle, d.n, d.d_r.get(), 1, d.d_z.get(), 1, &rz);

    for (int iter = 0; iter < maxIter_; ++iter) {
        // Ap = A * p
        double alpha_spmv = 1.0, beta_spmv = 0.0;
        cusparseSpMV(d.cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            &alpha_spmv, d.matDescr, d.dnVecP, &beta_spmv, d.dnVecAp,
            CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, d.d_spmvBuffer.get());

        double pAp = 0.0;
        cublasDdot(d.cublasHandle, d.n, d.d_p.get(), 1, d.d_Ap.get(), 1, &pAp);

        double alpha = rz / pAp;
        cublasDaxpy(d.cublasHandle, d.n, &alpha, d.d_p.get(), 1, d.d_x.get(), 1);

        double neg_alpha = -alpha;
        cublasDaxpy(d.cublasHandle, d.n, &neg_alpha, d.d_Ap.get(), 1, d.d_r.get(), 1);

        double rr_new = 0.0;
        cublasDdot(d.cublasHandle, d.n, d.d_r.get(), 1, d.d_r.get(), 1, &rr_new);
        if (rr_new < tol2) break;

        // z = M^{-1} r
        d.applyPreconditioner();

        double rz_new = 0.0;
        cublasDdot(d.cublasHandle, d.n, d.d_r.get(), 1, d.d_z.get(), 1, &rz_new);

        double beta_cg = rz_new / rz;
        cublasDscal(d.cublasHandle, d.n, &beta_cg, d.d_p.get(), 1);
        double one = 1.0;
        cublasDaxpy(d.cublasHandle, d.n, &one, d.d_z.get(), 1, d.d_p.get(), 1);

        rz = rz_new;
    }

    cudaMemcpy(x.data(), d.d_x.get(), d.n * sizeof(double), cudaMemcpyDeviceToHost);
    return x;
}

Eigen::MatrixXd CudaPCGSolver::solveMulti(const Eigen::MatrixXd& B) {
    Eigen::MatrixXd X(impl_->n, B.cols());
    for (int i = 0; i < B.cols(); ++i) {
        X.col(i) = solve(B.col(i));
    }
    return X;
}
