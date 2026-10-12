#include "backend.h"
#include <cstdint>
#include <iostream>
#include <new>

namespace backend {

namespace {
    bool initialized = false;
    bool use_pool = false;            //GPU memory comes from the stream-ordered memory pool
    int device = 0;                   //GPU used by this rank
#ifndef USE_CUDA
    constexpr std::align_val_t host_alignment{64};
#endif

    //Partial sums of the reductions. Created by init() as a function-local static, so it is destroyed
    //at the end of the program before the CUDA runtime shuts down.
    struct ReductionBuffers {
        static constexpr std::size_t bytes = max_reduce_blocks * 2 * sizeof(double);
        void* partial = allocate_comm_bytes(bytes);
        void* partial_host = allocate_pinned_bytes(bytes);
        ~ReductionBuffers(){
            free_comm_bytes(partial);
            free_pinned_bytes(partial_host);
        }
    };

    ReductionBuffers& reduction_buffers(){
        static ReductionBuffers buffers;
        return buffers;
    }
}

void init(){
    if (initialized) return;
    initialized = true;
#ifdef USE_CUDA
    //Ranks running on the same node use different GPUs (if there are fewer GPUs than ranks, they are shared)
    MPI_Comm node_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
    int local_rank;
    MPI_Comm_rank(node_comm, &local_rank);
    MPI_Comm_free(&node_comm);

    int num_devices = 0;
    CUDA_CHECK(cudaGetDeviceCount(&num_devices));
    device = local_rank % num_devices;
    CUDA_CHECK(cudaSetDevice(device));

    //Solvers allocate temporary fields at every call. With the stream-ordered memory pool, and a release
    //threshold that keeps the freed memory in the pool, these allocations are almost free.
    int pools_supported = 0;
    CUDA_CHECK(cudaDeviceGetAttribute(&pools_supported, cudaDevAttrMemoryPoolsSupported, device));
    if (pools_supported){
        cudaMemPool_t pool;
        CUDA_CHECK(cudaDeviceGetDefaultMemPool(&pool, device));
        std::uint64_t threshold = UINT64_MAX;
        CUDA_CHECK(cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &threshold));
        use_pool = true;
    }
#endif
    reduction_buffers();
}

std::string description(){
#ifdef USE_CUDA
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, device));
    std::string text = "CUDA, " + std::string(prop.name) + " (compute capability "
                     + std::to_string(prop.major) + "." + std::to_string(prop.minor) + ")";
    #ifdef CUDA_AWARE_MPI
    text += ", CUDA-aware MPI";
    #endif
    return text;
#else
    return "CPU";
#endif
}

void synchronize(){
#ifdef USE_CUDA
    CUDA_CHECK(cudaDeviceSynchronize());
#endif
}

void* allocate_bytes(std::size_t bytes){
    if (bytes == 0) return nullptr;
#ifdef USE_CUDA
    if (!initialized){
        //Fields cannot be created before initializeMPI(), otherwise every rank would use the first GPU
        std::cerr << "GPU memory requested before backend::init() (initializeMPI)" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    void* ptr = nullptr;
    if (use_pool)
        CUDA_CHECK(cudaMallocAsync(&ptr, bytes, 0));
    else
        CUDA_CHECK(cudaMalloc(&ptr, bytes));
    return ptr;
#else
    return ::operator new(bytes, host_alignment);
#endif
}

void free_bytes(void* ptr){
    if (ptr == nullptr) return;
#ifdef USE_CUDA
    //Errors are ignored on purpose: this is called by destructors, possibly when the program is exiting
    if (use_pool)
        cudaFreeAsync(ptr, 0);
    else
        cudaFree(ptr);
#else
    ::operator delete(ptr, host_alignment);
#endif
}

void copy_bytes(void* dst, const void* src, std::size_t bytes){
    if (bytes == 0) return;
#ifdef USE_CUDA
    CUDA_CHECK(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToDevice, 0));
#else
    std::memcpy(dst, src, bytes);
#endif
}

void zero_bytes(void* ptr, std::size_t bytes){
    if (bytes == 0) return;
#ifdef USE_CUDA
    CUDA_CHECK(cudaMemsetAsync(ptr, 0, bytes, 0));
#else
    std::memset(ptr, 0, bytes);
#endif
}

void copy_bytes_to_host(void* host_dst, const void* src, std::size_t bytes){
    if (bytes == 0) return;
#ifdef USE_CUDA
    CUDA_CHECK(cudaMemcpy(host_dst, src, bytes, cudaMemcpyDeviceToHost));
#else
    std::memcpy(host_dst, src, bytes);
#endif
}

void copy_bytes_from_host(void* dst, const void* host_src, std::size_t bytes){
    if (bytes == 0) return;
#ifdef USE_CUDA
    CUDA_CHECK(cudaMemcpy(dst, host_src, bytes, cudaMemcpyHostToDevice));
#else
    std::memcpy(dst, host_src, bytes);
#endif
}

void* allocate_comm_bytes(std::size_t bytes){
#ifdef USE_CUDA
    void* ptr = nullptr;
    CUDA_CHECK(cudaMalloc(&ptr, bytes));
    return ptr;
#else
    return ::operator new(bytes, host_alignment);
#endif
}

void* allocate_pinned_bytes(std::size_t bytes){
#ifdef USE_CUDA
    void* ptr = nullptr;
    CUDA_CHECK(cudaMallocHost(&ptr, bytes));
    return ptr;
#else
    return ::operator new(bytes, host_alignment);
#endif
}

void free_comm_bytes(void* ptr){
    if (ptr == nullptr) return;
#ifdef USE_CUDA
    cudaFree(ptr); //errors ignored, see free_bytes
#else
    ::operator delete(ptr, host_alignment);
#endif
}

void free_pinned_bytes(void* ptr){
    if (ptr == nullptr) return;
#ifdef USE_CUDA
    cudaFreeHost(ptr); //errors ignored, see free_bytes
#else
    ::operator delete(ptr, host_alignment);
#endif
}

void* reduction_buffer(){ return reduction_buffers().partial; }
void* reduction_buffer_host(){ return reduction_buffers().partial_host; }

} //namespace backend
