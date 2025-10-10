cmake_minimum_required(VERSION 3.26.5)
set(CMAKE_CXX_STANDARD 20) 
set(CMAKE_CXX_STANDARD_REQUIRED True)
project(CG VERSION 1.0)


set(SOURCES
    src/variables.cpp           #Relevant parameters for the lattice
    src/main.cpp                #Main program
    src/gauge_conf.cpp          #Gauge configuration (computes staples, plaquettes and gauge action)
    src/dirac_operator.cpp      #Implementation of D, D^+ and 2* Re ( left^+ d D / d omega(z) right )
    src/conjugate_gradient.cpp  #Conjugate gradient and Bi-cgstab
)

set(NS "512")   #Lattice dimension x
set(NT "512")   #Lattice dimension t

configure_file(${CMAKE_SOURCE_DIR}/include/config.h.in ${CMAKE_SOURCE_DIR}/include/config.h)


add_executable("CG_${NS}x${NT}" ${SOURCES})


target_include_directories("CG_${NS}x${NT}"
        PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)
target_compile_options("CG_${NS}x${NT}" PRIVATE
			-O3
            #-g #For profiling
)

find_package(MPI REQUIRED)
target_link_libraries("CG_${NS}x${NT}" PUBLIC MPI::MPI_CXX)

