cmake_minimum_required(VERSION 3.26.5)
set(CMAKE_CXX_STANDARD 20) 
set(CMAKE_CXX_STANDARD_REQUIRED True)
project(Schwinger_model VERSION 1.0)

option(TWISTED_MASS "Use twisted mass fermions" OFF)
option(CLOVER "Use clover term" OFF)
option(BUILD_TESTS "Build unit tests" OFF)

set(NS "64")   #Lattice dimension x
set(NT "64")   #Lattice dimension t
configure_file(${CMAKE_SOURCE_DIR}/include/config.h.in ${CMAKE_SOURCE_DIR}/include/config.h)

find_package(MPI REQUIRED)

#------We create a common library with src files shared by different executables
add_library(common_lib
    src/utils.cpp               #Utility functions
    src/variables.cpp           #Relevant parameters for the lattice
    src/gauge_conf.cpp          #Gauge configuration (computes staples, plaquettes and gauge action)
    src/dirac_operator.cpp      #Implementation of D, D^+ and 2* Re ( left^+ d D / d omega(z) right )
    src/halo_exchange.cpp       #Halo exchange
    src/conjugate_gradient.cpp  #Conjugate gradient
    src/bi_cgstab.cpp           #BiCGstab
    src/hmc.cpp                 #Hybrid Monte Carlo
)

target_include_directories(common_lib
    PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(common_lib
    PUBLIC
    MPI::MPI_CXX
)

target_compile_options(common_lib PRIVATE
    -O3
    #-g #for profiling
)

if(TWISTED_MASS)
    target_compile_definitions(common_lib PUBLIC TWISTED_MASS)
endif()

if(CLOVER)
    target_compile_definitions(common_lib PUBLIC CLOVER)
endif()


#HMC simulation
add_executable("SM_${NS}x${NT}"
    src/main.cpp
)

target_link_libraries("SM_${NS}x${NT}"
    PRIVATE
    common_lib
)

#Pion and PCAC correlators
add_executable("mass_${NS}x${NT}" 
    src/mass.cpp
)

target_link_libraries("mass_${NS}x${NT}"
    PRIVATE
    common_lib
)

#-----------------Testing with Catch2-----------------#

if(BUILD_TESTS)

    include(FetchContent)

    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.8.1
    )

    FetchContent_MakeAvailable(Catch2)

    #Dirac operator
    add_executable(tests_D
        tests/test_dirac_operator.cpp
    )

    target_link_libraries(tests_D
    PRIVATE
        common_lib
        Catch2::Catch2WithMain
    )

    #Test CG
    add_executable(tests_CG
        tests/tests_CG.cpp
    )

    target_link_libraries(tests_CG
    PRIVATE
        common_lib
        Catch2::Catch2WithMain
    )

    #Test gauge conf
    add_executable(tests_gconf
        tests/tests_gconf.cpp
    )

    target_link_libraries(tests_gconf
    PRIVATE
        common_lib
        Catch2::Catch2WithMain
    )

    #Test HMC
    add_executable(tests_hmc
        tests/tests_hmc.cpp
    )

    target_link_libraries(tests_hmc
    PRIVATE
        common_lib
        Catch2::Catch2WithMain
    )

endif()