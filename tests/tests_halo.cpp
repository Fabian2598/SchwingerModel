#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include "halo_exchange.h"

/*
    The interior sites of a field store their global coordinates, (gx + 1000*c, gt) for component c.
    After the halo exchange, every halo site has to contain the value of its periodic image on the global lattice.
    This can be run with 1 rank (periodic copies inside the rank) or 4 ranks (MPI communication).
*/

c_double label(const int gx, const int gt, const int c){
    return c_double(gx + 1000.0*c, gt);
}

//Field with dof complex numbers per site, labelled with the global coordinates of the sites
spinor labelled_field(const int dof){
    spinor field(dof*mpi::sitesH);
    std::vector<c_double> field_host(field.size);
    for(int x = 1; x<=mpi::width_x; x++){
        for(int t = 1; t<=mpi::width_t; t++){
            int n = x*(mpi::width_t+2)+t;
            int gx = mpi::coords[0]*mpi::width_x + x-1;
            int gt = mpi::coords[1]*mpi::width_t + t-1;
            for(int c = 0; c < dof; c++)
                field_host[dof*n+c] = label(gx, gt, c);
        }
    }
    field.fromHost(field_host);
    return field;
}

//Checks the halo (corners only if requested) on all ranks
bool halo_is_correct(const spinor& field, const int dof, const bool corners){
    std::vector<c_double> field_host = field.toHost();
    int correct = 1;
    for(int x = 0; x<=mpi::width_x+1; x++){
        for(int t = 0; t<=mpi::width_t+1; t++){
            bool x_halo = (x == 0 || x == mpi::width_x+1);
            bool t_halo = (t == 0 || t == mpi::width_t+1);
            if (!x_halo && !t_halo) continue;           //interior site
            if (x_halo && t_halo && !corners) continue; //corner
            int n = x*(mpi::width_t+2)+t;
            int gx = mod(mpi::coords[0]*mpi::width_x + x-1, LV::Nx);
            int gt = mod(mpi::coords[1]*mpi::width_t + t-1, LV::Nt);
            for(int c = 0; c < dof; c++)
                if (field_host[dof*n+c] != label(gx, gt, c))
                    correct = 0;
        }
    }
    int global_correct;
    MPI_Allreduce(&correct, &global_correct, 1, MPI_INT, MPI_LAND, mpi::cart_comm);
    return global_correct == 1;
}

TEST_CASE("Halo exchange of a spinor", "[halo]") {
    spinor phi = labelled_field(2);
    exchange_halo(phi.val);
    REQUIRE(halo_is_correct(phi, 2, false));
}

TEST_CASE("Halo exchange of a vector, corners included", "[halo]") {
    spinor v = labelled_field(1);
    exchange_halo_vec(v.val);
    REQUIRE(halo_is_correct(v, 1, true));
}

TEST_CASE("Halo exchange of the gauge links, corners included", "[halo]") {
    spinor U = labelled_field(2);
    exchange_halo_gauge(U.val);
    REQUIRE(halo_is_correct(U, 2, true));
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::size);

    if (mpi::size != 1 && mpi::size != 4){
        if (mpi::rank == 0)
            std::cerr << "ERROR: This test is meant to be run with 1 or 4 ranks" << std::endl;
        MPI_Finalize();
        return 0;
    }
    mpi::ranks_x = (mpi::size == 1) ? 1 : 2;
    mpi::ranks_t = (mpi::size == 1) ? 1 : 2;
    initializeMPI();

    int result = Catch::Session().run(argc, argv);
    MPI_Finalize();

    return result;
}
