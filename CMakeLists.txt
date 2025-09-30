cmake_minimum_required(VERSION 3.26.5)
set(CMAKE_CXX_STANDARD 20) 
set(CMAKE_CXX_STANDARD_REQUIRED True)
project(Schwinger_model VERSION 1.0)


set(SOURCES
    src/statistics.cpp          #Functions to perform statistical analysis
    src/variables.cpp           #Relevant parameters for the lattice
    src/main.cpp                #Main program
    src/gauge_conf.cpp          #Gauge configuration (computes staples, plaquettes and gauge action)
    src/dirac_operator.cpp      #Implementation of D, D^+ and 2* Re ( left^+ d D / d omega(z) right )
    src/conjugate_gradient.cpp  #Conjugate gradient and Bi-cgstab
    src/hmc.cpp                 #Hybrid Monte Carlo
)

set(NS "64")   #Lattice dimension x
set(NT "64")   #Lattice dimension t

configure_file(${CMAKE_SOURCE_DIR}/include/config.h.in ${CMAKE_SOURCE_DIR}/include/config.h)


add_executable("SM_${NS}x${NT}" ${SOURCES})


target_include_directories("SM_${NS}x${NT}"
        PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)
target_compile_options("SM_${NS}x${NT}" PRIVATE
			-O3
            #-g #For profiling
			-fopenmp
)

find_package(OpenMP REQUIRED)
target_link_libraries("SM_${NS}x${NT}" PUBLIC OpenMP::OpenMP_CXX)

find_package(MPI REQUIRED)
target_link_libraries("SM_${NS}x${NT}" PUBLIC MPI::MPI_CXX)



#---------------------------------------#

#PROJECT_SOURCE_DIR: it indicates where src is located
#cmake -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER="C:/msys64/ucrt64/bin/g++.exe" -DCMAKE_C_COMPILER="C:/msys64/ucrt64/bin/gcc.exe" ../
#dir/s/b *.txt > confs.txt

#My workstation has 12 physical cores. MPI ranks * OMP threads 
#should not exceed that number
#export OMP_NUM_THREADS=2
#export OMP_PROC_BIND=close
#export OMP_PLACES=cores

#mpirun --bind-to core --map-by core -np 4 ./SM_64x64