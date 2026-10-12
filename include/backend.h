#ifndef BACKEND_H
#define BACKEND_H

/*
    Execution backend: CUDA (GPU) or CPU.

    Every loop over the lattice sites of a rank is written once, as a lambda, and executed with
        for_each_site( [=] HOST_DEVICE (int n){ ... } );                 //see utils.h
        double s = lattice_sum<double>( [=] HOST_DEVICE (int n){ return ...; } );
    where n = x*(width_t+2)+t is the local site index (halo included). When the code is compiled with
    -DUSE_CUDA=ON these loops are CUDA kernels executed by the GPU assigned to the rank, with one thread
    per lattice site. Otherwise they are ordinary sequential loops on the CPU.

    Rules for the lambdas (with CUDA they are executed on the GPU):
      - Capture by value ([=]) only raw pointers to the data of the fields (field.val) and plain values
        (int, double, c_double, LocalLattice). Fields (spinor_template), global variables (mpi::, sim_params::)
        and class members (this->...) are not accessible on the GPU: copy them to local variables first.
      - Only HOST_DEVICE functions can be called inside (c_double arithmetic, get_neighbors, sin, cos, ...).

    With CUDA the fields live in GPU memory, so field.val can only be dereferenced inside these lambdas.
    Use field.toHost() and field.fromHost() to read or write a field on the CPU (I/O, random numbers, tests).
*/

#include <cstddef>
#include <cstdio>
#include <cstring>
#include <string>
#include "mpi.h"

#ifdef USE_CUDA
    #ifndef __CUDACC__
        #error "USE_CUDA requires all the source files to be compiled with nvcc (LANGUAGE CUDA in CMakeLists.txt)"
    #endif
    #include <cuda_runtime.h>
    #define HOST_DEVICE __host__ __device__
    //Abort all the ranks if a CUDA call fails
    #define CUDA_CHECK(call)                                                                    \
        do {                                                                                    \
            cudaError_t cuda_check_err = (call);                                                \
            if (cuda_check_err != cudaSuccess) {                                                \
                std::fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,           \
                             cudaGetErrorString(cuda_check_err));                               \
                MPI_Abort(MPI_COMM_WORLD, 1);                                                   \
            }                                                                                   \
        } while (0)
#else
    #define HOST_DEVICE
#endif

namespace backend {

constexpr int block_size = 128;         //threads per block of the GPU kernels
constexpr int max_reduce_blocks = 1024; //maximum number of blocks of a reduction kernel

/*
    Selects the GPU of this rank (node-local rank modulo the number of visible GPUs) and allocates the
    buffers for the reductions. initializeMPI() calls it, so it runs before any field is allocated.
    It does nothing in the CPU version.
*/
void init();

//Description of the backend for the output, e.g. "CUDA, NVIDIA A100-SXM4-40GB" or "CPU"
std::string description();

//Waits until the GPU has finished all the work submitted so far (nothing to do on the CPU)
void synchronize();

//--------Memory in the execution space (GPU memory with CUDA, host memory otherwise)--------//
void* allocate_bytes(std::size_t bytes);
void free_bytes(void* ptr);
void copy_bytes(void* dst, const void* src, std::size_t bytes); //both pointers in the execution space
void zero_bytes(void* ptr, std::size_t bytes);
void copy_bytes_to_host(void* host_dst, const void* src, std::size_t bytes);
void copy_bytes_from_host(void* dst, const void* host_src, std::size_t bytes);

//Buffers that live during the whole run, e.g. the ones handed to MPI. With CUDA they are allocated outside
//of the memory pool used for the fields, since some CUDA-aware MPI libraries cannot use memory from the pool.
void* allocate_comm_bytes(std::size_t bytes);
void free_comm_bytes(void* ptr);
//Page-locked host memory (fast copies between GPU and CPU). Ordinary host memory in the CPU version.
void* allocate_pinned_bytes(std::size_t bytes);
void free_pinned_bytes(void* ptr);

//Partial sums of the reductions: max_reduce_blocks values of at most 16 bytes
void* reduction_buffer();      //execution space
void* reduction_buffer_host(); //host

#ifdef USE_CUDA
//----------------------------------------GPU kernels----------------------------------------//

//f(n) for every site of a width_x x width_t local lattice (halo excluded), one thread per site
template <typename F>
__global__ void sites_kernel(const int width_x, const int width_t, F f){
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < width_x * width_t){
        const int x = i / width_t + 1;
        const int t = i % width_t + 1;
        f(x * (width_t + 2) + t);
    }
}

//f(i) for i = 0, ..., N-1
template <typename F>
__global__ void range_kernel(const int N, F f){
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N)
        f(i);
}

//Warp shuffle for types made of doubles (double and c_double)
template <typename T>
__device__ T shfl_down(T v, const int offset){
    constexpr int ndoubles = sizeof(T) / sizeof(double);
    double d[ndoubles];
    memcpy(d, &v, sizeof(T));
    for (int k = 0; k < ndoubles; k++)
        d[k] = __shfl_down_sync(0xffffffffu, d[k], offset);
    memcpy(&v, d, sizeof(T));
    return v;
}

//Sum of v over the threads of a block. The result is only valid in thread 0.
template <typename T>
__device__ T block_sum(T v){
    __shared__ __align__(16) unsigned char buffer[32 * sizeof(T)];
    T* warp_sums = reinterpret_cast<T*>(buffer);
    const int lane = threadIdx.x % 32;
    const int warp = threadIdx.x / 32;
    for (int offset = 16; offset > 0; offset /= 2)
        v += shfl_down(v, offset);
    if (lane == 0)
        warp_sums[warp] = v;
    __syncthreads();
    if (warp == 0){
        v = (lane < blockDim.x / 32) ? warp_sums[lane] : T();
        for (int offset = 16; offset > 0; offset /= 2)
            v += shfl_down(v, offset);
    }
    return v;
}

//partial[block] = sum of f(n) over the sites handled by the block
template <typename T, typename F>
__global__ void sum_sites_kernel(const int width_x, const int width_t, F f, T* partial){
    T sum = T();
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < width_x * width_t; i += gridDim.x * blockDim.x){
        const int x = i / width_t + 1;
        const int t = i % width_t + 1;
        sum += f(x * (width_t + 2) + t);
    }
    sum = block_sum(sum);
    if (threadIdx.x == 0)
        partial[blockIdx.x] = sum;
}
#endif

//----------------------------------------Parallel loops----------------------------------------//

//f(n) for every site n = x*(width_t+2)+t, x = 1, ..., width_x, t = 1, ..., width_t
template <typename F>
void for_sites(const int width_x, const int width_t, F f){
#ifdef USE_CUDA
    const int N = width_x * width_t;
    if (N == 0) return;
    sites_kernel<<<(N + block_size - 1) / block_size, block_size>>>(width_x, width_t, f);
    CUDA_CHECK(cudaGetLastError());
#else
    for (int x = 1; x <= width_x; x++)
        for (int t = 1; t <= width_t; t++)
            f(x * (width_t + 2) + t);
#endif
}

//f(i) for i = 0, ..., N-1
template <typename F>
void for_range(const int N, F f){
#ifdef USE_CUDA
    if (N == 0) return;
    range_kernel<<<(N + block_size - 1) / block_size, block_size>>>(N, f);
    CUDA_CHECK(cudaGetLastError());
#else
    for (int i = 0; i < N; i++)
        f(i);
#endif
}

/*
    Sum of f(n) over the sites of the local lattice (T = double or c_double).
    On the GPU every block computes a partial sum, and the partial sums are added up on the CPU in a fixed
    order, so the result does not change from run to run.
*/
template <typename T, typename F>
T sum_sites(const int width_x, const int width_t, F f){
    static_assert(sizeof(T) % sizeof(double) == 0 && sizeof(T) <= 2 * sizeof(double),
                  "sum_sites supports double and c_double");
    T sum = T();
#ifdef USE_CUDA
    const int N = width_x * width_t;
    if (N == 0) return sum;
    int blocks = (N + block_size - 1) / block_size;
    if (blocks > max_reduce_blocks) blocks = max_reduce_blocks;
    T* partial = static_cast<T*>(reduction_buffer());
    sum_sites_kernel<T><<<blocks, block_size>>>(width_x, width_t, f, partial);
    CUDA_CHECK(cudaGetLastError());
    T* partial_host = static_cast<T*>(reduction_buffer_host());
    CUDA_CHECK(cudaMemcpy(partial_host, partial, blocks * sizeof(T), cudaMemcpyDeviceToHost));
    for (int b = 0; b < blocks; b++)
        sum += partial_host[b];
#else
    for (int x = 1; x <= width_x; x++)
        for (int t = 1; t <= width_t; t++)
            sum += f(x * (width_t + 2) + t);
#endif
    return sum;
}

} //namespace backend

#endif
