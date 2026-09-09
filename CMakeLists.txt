cmake_minimum_required(VERSION 3.26.5)
set(CMAKE_CXX_STANDARD 20) 
set(CMAKE_CXX_STANDARD_REQUIRED True)
project(Schwinger_model VERSION 1.0)

option(TWISTED_MASS "Use twisted mass fermions" OFF)
option(CLOVER "Use clover term" OFF)
option(BUILD_TESTS "Build unit tests" OFF)


set(COMMON_SOURCES
    src/utils.cpp               #Utility functions
    src/variables.cpp           #Relevant parameters for the lattice
    src/gauge_conf.cpp          #Gauge configuration (computes staples, plaquettes and gauge action)
    src/dirac_operator.cpp      #Implementation of D, D^+ and 2* Re ( left^+ d D / d omega(z) right )
    src/halo_exchange.cpp       #Halo exchange
    src/conjugate_gradient.cpp  #Conjugate gradient
)

set(SIMULATION_SOURCES
    ${COMMON_SOURCES}
    src/hmc.cpp                 #Hybrid Monte Carlo
    src/main.cpp
)

set(NS "64")   #Lattice dimension x
set(NT "64")   #Lattice dimension t

configure_file(${CMAKE_SOURCE_DIR}/include/config.h.in ${CMAKE_SOURCE_DIR}/include/config.h)


add_executable("SM_${NS}x${NT}" ${SIMULATION_SOURCES} )

if(TWISTED_MASS)
    target_compile_definitions("SM_${NS}x${NT}" PRIVATE TWISTED_MASS)
endif()

if(CLOVER)
    target_compile_definitions("SM_${NS}x${NT}" PRIVATE CLOVER)
endif()

target_include_directories("SM_${NS}x${NT}"
        PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)
target_compile_options("SM_${NS}x${NT}" PRIVATE
			-O3
            #-g #For profiling
)


find_package(MPI REQUIRED)
target_link_libraries("SM_${NS}x${NT}" PUBLIC MPI::MPI_CXX)


#PROJECT_SOURCE_DIR: it indicates where src is located
#cmake -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER="C:/msys64/ucrt64/bin/g++.exe" -DCMAKE_C_COMPILER="C:/msys64/ucrt64/bin/gcc.exe" ../
#dir/s/b *.txt > confs.txt


#-----------------For measuring the pion correlator-----------------#
set(SOURCES_MASS
    ${COMMON_SOURCES}
    src/mass.cpp                #Main program for computing pion and pcac correlators
    src/bi_cgstab.cpp           #Bi-cgstab
)

add_executable("mass_${NS}x${NT}" ${SOURCES_MASS})

target_include_directories("mass_${NS}x${NT}"
        PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)
target_compile_options("mass_${NS}x${NT}" PRIVATE
			-O3
)

find_package(MPI REQUIRED)
target_link_libraries("mass_${NS}x${NT}" PUBLIC MPI::MPI_CXX)


#-----------------Testing with Catch2-----------------#

if(BUILD_TESTS)

    include(FetchContent)

    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.8.1
    )

    FetchContent_MakeAvailable(Catch2)

    add_executable(tests_D
     ${COMMON_SOURCES}
     tests/test_dirac_operator.cpp
     )

    target_include_directories(tests_D
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
    )

    target_compile_options(tests_D PRIVATE
        -O3
    )

    add_executable(tests_CG ${TEST_SOURCES_CG})

    target_include_directories(tests_CG
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
    )

    target_compile_options(tests_CG PRIVATE
        -O3
    )

    if(TWISTED_MASS)
        target_compile_definitions(tests_D PRIVATE TWISTED_MASS)
        target_compile_definitions(tests_CG PRIVATE TWISTED_MASS)
    endif()

    if(CLOVER)
        target_compile_definitions(tests_D PRIVATE CLOVER)
        target_compile_definitions(tests_CG PRIVATE CLOVER)
    endif()

    target_link_libraries(tests_D
        PRIVATE
        MPI::MPI_CXX
        Catch2::Catch2WithMain
    )

    target_link_libraries(tests_CG
        PRIVATE
        MPI::MPI_CXX
        Catch2::Catch2WithMain
    )

    

    include(CTest)
    include(Catch)

endif()