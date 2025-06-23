cmake_minimum_required(VERSION 3.26.5)
set(CMAKE_CXX_STANDARD 11) 
set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(CMAKE_C_COMPILER "C:/msys64/ucrt64/bin/gcc.exe")
#set(CMAKE_CXX_COMPILER "C:/msys64/ucrt64/bin/g++.exe")
project(Schwinger_model VERSION 1.0)

# Create a sources variable with a link to all cpp files to compile
set(SRC_FILES 
    src/variables.cpp
    src/gauge_conf.cpp
    src/dirac_operator.cpp
    src/bi_cgstab.cpp
    src/sap.cpp
    src/amg.cpp
    src/main.cpp  
    src/statistics.cpp
    src/gmres.cpp
    src/fgmres.cpp
)
#file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR} src/*.cpp)

#Lattice dimensions
set(NS "16") 
set(NT "16") 

#Lattice blocks (used for aggregation)
#N_i / BLOCK_i has to be an integer number
set(BLOCK_X "2")
set(BLOCK_T "2")
#The number of aggregates is 2 * BlockX * BlockT 

#Lattice blocks for Schwarz alternating method. They don't need to be the same as BLOCK_X and BLOCK_T
set(SAP_BLOCK_X "4")
set(SAP_BLOCK_T "4")
#N_i / SAP_BLOCK_i has to be an integer number
#One block contains 2 * N_x / SAP_BLOCK_x * N_t / SAP_BLOCK_t variables 
#For parallelization we fix one block per core (i.e. SAP_BLOCK_X * SAP_BLOCK_T / 2 cores)

#Number of test vectors
set(NTEST "20") 
configure_file(${CMAKE_SOURCE_DIR}/include/config.h.in ${CMAKE_SOURCE_DIR}/include/config.h)


add_executable("SM_${NS}x${NT}" ${SRC_FILES})

target_include_directories("SM_${NS}x${NT}"
        PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)

target_compile_options("SM_${NS}x${NT}" PRIVATE
            #-Wall
			-O3
			#-fopenmp
)

find_package(MPI REQUIRED)

#find_package(BLAS)
#include_directories(" /usr/include/cblas/")
#target_link_libraries("SM_${NS}x${NT}" PUBLIC MPI::MPI_CXX ${BLAS_LIBRARIES})

target_link_libraries("SM_${NS}x${NT}" PUBLIC MPI::MPI_CXX)

#find_package(OpenMP REQUIRED)
#target_link_libraries("SM_${NS}x${NT}" PUBLIC OpenMP::OpenMP_CXX)


#PROJECT_SOURCE_DIR: it indicates where src is located
#cmake -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER="C:/msys64/ucrt64/bin/g++.exe" -DCMAKE_C_COMPILER="C:/msys64/ucrt64/bin/gcc.exe" ../


