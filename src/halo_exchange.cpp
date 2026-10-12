#include "halo_exchange.h"

/*
    The halo is filled in segments: the four edges and, if needed, the four corners. Segment s contains sites
    on the border of the local lattice that are sent to the rank dest[s]; at the same time we receive segment s
    from the rank on the opposite side, src[s], and store it in our halo (check boundary.h for the conventions):

      segment   inner sites sent       to          received from    stored in the halo
         0      row x = 1              top         bot              row x = width_x+1
         1      row x = width_x        bot         top              row x = 0
         2      column t = 1           left        right            column t = width_t+1
         3      column t = width_t     right       left             column t = 0
         4      (width_x, 1)           bot_left    top_right        (0, width_t+1)
         5      (1, width_t)           top_right   bot_left         (width_x+1, 0)
         6      (width_x, width_t)     bot_right   top_left         (0, 0)
         7      (1, 1)                 top_left    bot_right        (width_x+1, width_t+1)

    All the segments are packed into one contiguous buffer, so MPI only handles contiguous memory and, with CUDA,
    there is a single copy GPU -> CPU and a single copy CPU -> GPU per exchange (none with CUDA-aware MPI).
    With a single rank every neighbor is the rank itself, and the halo is filled directly (periodic copy).
*/

namespace {

constexpr int max_segments = 8;
constexpr int max_dof = 2; //maximum number of complex numbers per site

//Size and position of the segments in the buffers
struct HaloLayout {
    int width_x, width_t;
    int dof;                  //complex numbers per site
    int segments;             //4 (edges) or 8 (edges and corners)
    int count[max_segments];  //sites in each segment
    int offset[max_segments]; //position (in sites) of each segment in the buffers
    int sites;                //total number of sites
};

HaloLayout make_layout(const int dof, const bool corners){
    HaloLayout h;
    h.width_x = mpi::width_x;
    h.width_t = mpi::width_t;
    h.dof = dof;
    h.segments = corners ? 8 : 4;
    const int count[max_segments] = {mpi::width_t, mpi::width_t, mpi::width_x, mpi::width_x, 1, 1, 1, 1};
    h.sites = 0;
    for(int s = 0; s < max_segments; s++){
        h.count[s] = (s < h.segments) ? count[s] : 0;
        h.offset[s] = h.sites;
        h.sites += h.count[s];
    }
    return h;
}

//inner: site in position k of segment s that is sent; halo: site of the halo where it is stored by the receiver
HOST_DEVICE inline void segment_sites(const HaloLayout& h, const int s, const int k, int& inner, int& halo){
    const int row = h.width_t+2;
    const int wx = h.width_x, wt = h.width_t;
    switch(s){
        case 0:  inner = 1*row + k+1;     halo = (wx+1)*row + k+1;     break;
        case 1:  inner = wx*row + k+1;    halo = k+1;                  break;
        case 2:  inner = (k+1)*row + 1;   halo = (k+1)*row + wt+1;     break;
        case 3:  inner = (k+1)*row + wt;  halo = (k+1)*row;            break;
        case 4:  inner = wx*row + 1;      halo = wt+1;                 break;
        case 5:  inner = 1*row + wt;      halo = (wx+1)*row;           break;
        case 6:  inner = wx*row + wt;     halo = 0;                    break;
        default: inner = 1*row + 1;       halo = (wx+1)*row + wt+1;    break;
    }
}

//Position i of the buffers (in complex numbers) -> segment s, site k inside the segment, component c
HOST_DEVICE inline void locate(const HaloLayout& h, const int i, int& s, int& k, int& c){
    const int site = i / h.dof;
    c = i % h.dof;
    s = 0;
    while (s < h.segments-1 && site >= h.offset[s+1])
        s++;
    k = site - h.offset[s];
}

//Communication buffers for the largest exchange (max_dof complex numbers per site, corners included)
struct HaloBuffers {
    const int capacity = max_dof * make_layout(max_dof, true).sites; //in complex numbers
    c_double* send = static_cast<c_double*>(backend::allocate_comm_bytes(capacity * sizeof(c_double)));
    c_double* recv = static_cast<c_double*>(backend::allocate_comm_bytes(capacity * sizeof(c_double)));
    #if defined(USE_CUDA) && !defined(CUDA_AWARE_MPI)
    //Host copies of the buffers for a MPI library that is not CUDA-aware
    c_double* send_host = static_cast<c_double*>(backend::allocate_pinned_bytes(capacity * sizeof(c_double)));
    c_double* recv_host = static_cast<c_double*>(backend::allocate_pinned_bytes(capacity * sizeof(c_double)));
    #endif

    ~HaloBuffers(){
        backend::free_comm_bytes(send);
        backend::free_comm_bytes(recv);
        #if defined(USE_CUDA) && !defined(CUDA_AWARE_MPI)
        backend::free_pinned_bytes(send_host);
        backend::free_pinned_bytes(recv_host);
        #endif
    }
};

//Allocated at the first exchange between ranks, freed at the end of the program
HaloBuffers& halo_buffers(){
    static HaloBuffers buffers;
    return buffers;
}

//Halo exchange of a field with dof complex numbers per site
void exchange(c_double* field, const int dof, const bool corners){
    if (dof > max_dof){
        std::cerr << "exchange_halo: at most " << max_dof << " complex numbers per site" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    const HaloLayout h = make_layout(dof, corners);
    const int n = h.sites * dof; //complex numbers to exchange

    if (mpi::size == 1){
        //All the neighbors are this rank: copy the border of the lattice to the opposite halo
        backend::for_range(n, [=] HOST_DEVICE (int i){
            int s, k, c, inner, halo;
            locate(h, i, s, k, c);
            segment_sites(h, s, k, inner, halo);
            field[halo*dof + c] = field[inner*dof + c];
        });
        return;
    }

    HaloBuffers& buffers = halo_buffers();

    //Pack the border of the local lattice
    c_double* send = buffers.send;
    c_double* recv = buffers.recv;
    backend::for_range(n, [=] HOST_DEVICE (int i){
        int s, k, c, inner, halo;
        locate(h, i, s, k, c);
        segment_sites(h, s, k, inner, halo);
        send[i] = field[inner*dof + c];
    });

    #if defined(USE_CUDA) && !defined(CUDA_AWARE_MPI)
    //MPI works on host copies of the buffers
    backend::copy_bytes_to_host(buffers.send_host, send, n * sizeof(c_double));
    c_double* send_mpi = buffers.send_host;
    c_double* recv_mpi = buffers.recv_host;
    #else
    //MPI reads and writes the buffers directly (GPU memory with CUDA-aware MPI)
    backend::synchronize(); //the packing has to be finished
    c_double* send_mpi = send;
    c_double* recv_mpi = recv;
    #endif

    const int dest[max_segments] = {mpi::top, mpi::bot, mpi::left, mpi::right,
                                    mpi::bot_left, mpi::top_right, mpi::bot_right, mpi::top_left};
    const int src[max_segments]  = {mpi::bot, mpi::top, mpi::right, mpi::left,
                                    mpi::top_right, mpi::bot_left, mpi::top_left, mpi::bot_right};
    MPI_Request requests[2*max_segments];
    for(int s = 0; s < h.segments; s++)
        MPI_Irecv(recv_mpi + h.offset[s]*dof, h.count[s]*dof, MPI_DOUBLE_COMPLEX, src[s], s,
                  mpi::cart_comm, &requests[s]);
    for(int s = 0; s < h.segments; s++)
        MPI_Isend(send_mpi + h.offset[s]*dof, h.count[s]*dof, MPI_DOUBLE_COMPLEX, dest[s], s,
                  mpi::cart_comm, &requests[h.segments + s]);
    MPI_Waitall(2*h.segments, requests, MPI_STATUSES_IGNORE);

    #if defined(USE_CUDA) && !defined(CUDA_AWARE_MPI)
    backend::copy_bytes_from_host(recv, buffers.recv_host, n * sizeof(c_double));
    #endif

    //Unpack into the halo
    backend::for_range(n, [=] HOST_DEVICE (int i){
        int s, k, c, inner, halo;
        locate(h, i, s, k, c);
        segment_sites(h, s, k, inner, halo);
        field[halo*dof + c] = recv[i];
    });
}

} //namespace


void exchange_halo(c_double* phi){
    exchange(phi, 2, false);
}

void exchange_halo_vec(c_double* v){
    exchange(v, 1, true);
}

void exchange_halo_gauge(c_double* U){
    exchange(U, 2, true);
}
