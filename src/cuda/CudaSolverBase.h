#ifndef _CUDA_SOLVER_BASE_H_
#define _CUDA_SOLVER_BASE_H_

#include "CSRMatrix.h"
#include <iostream>

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cusolverSp.h>
#include <cublas_v2.h>

// CUDA error checking macros
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = (call); \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ \
                      << " - " << cudaGetErrorString(err) << std::endl; \
            return false; \
        } \
    } while(0)

#define CUSPARSE_CHECK(call) \
    do { \
        cusparseStatus_t status = (call); \
        if (status != CUSPARSE_STATUS_SUCCESS) { \
            std::cerr << "cuSPARSE error at " << __FILE__ << ":" << __LINE__ \
                      << " - code " << status << std::endl; \
            return false; \
        } \
    } while(0)

#define CUSOLVER_CHECK(call) \
    do { \
        cusolverStatus_t status = (call); \
        if (status != CUSOLVER_STATUS_SUCCESS) { \
            std::cerr << "cuSOLVER error at " << __FILE__ << ":" << __LINE__ \
                      << " - code " << status << std::endl; \
            return false; \
        } \
    } while(0)

#define CUBLAS_CHECK(call) \
    do { \
        cublasStatus_t status = (call); \
        if (status != CUBLAS_STATUS_SUCCESS) { \
            std::cerr << "cuBLAS error at " << __FILE__ << ":" << __LINE__ \
                      << " - code " << status << std::endl; \
            return false; \
        } \
    } while(0)

// RAII wrapper for device memory
template<typename T>
class DeviceBuffer {
private:
    T* ptr_ = nullptr;
    size_t size_ = 0;

public:
    DeviceBuffer() = default;
    ~DeviceBuffer() { free(); }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    DeviceBuffer(DeviceBuffer&& other) noexcept : ptr_(other.ptr_), size_(other.size_) {
        other.ptr_ = nullptr;
        other.size_ = 0;
    }

    DeviceBuffer& operator=(DeviceBuffer&& other) noexcept {
        if (this != &other) {
            free();
            ptr_ = other.ptr_;
            size_ = other.size_;
            other.ptr_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }

    bool alloc(size_t count) {
        free();
        size_ = count;
        cudaError_t err = cudaMalloc(&ptr_, count * sizeof(T));
        if (err != cudaSuccess) {
            ptr_ = nullptr;
            size_ = 0;
            return false;
        }
        return true;
    }

    void free() {
        if (ptr_) {
            cudaFree(ptr_);
            ptr_ = nullptr;
            size_ = 0;
        }
    }

    T* get() { return ptr_; }
    const T* get() const { return ptr_; }
    size_t size() const { return size_; }
};

// CSRMatrix and conversion functions declared in CSRMatrix.h

#endif // _CUDA_SOLVER_BASE_H_
