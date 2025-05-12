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
    src/matrix_operations.cpp
    src/conjugate_gradient.cpp
    src/sap.cpp
    src/amg.cpp
    src/main.cpp  
    src/statistics.cpp
    src/gmres.cpp
)
#file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR} src/*.cpp)

#Lattice dimensions
set(NS "64") 
set(NT "64") 

#Number of lattice blocks (used for aggregation)
#N_i / BLOCK_i has to be an integer number
set(BLOCK_X "2")
set(BLOCK_T "2")

#Lattice blocks for Schwarz alternating method. They don't need to be the same as BLOCK_X and BLOCK_T
set(SAP_BLOCK_X "4")
set(SAP_BLOCK_T "4")
#N_i / SAP_BLOCK_i has to be an integer number
#One block contains 2 * N_x / SAP_BLOCK_x * N_t / SAP_BLOCK_t variables 

#Number of test vectors
set(NTEST "10") 
set(NAGG "8") #BlockX * BlockT for aggregation scheme 1 and 2 * BlockX * BlockT for scheme 2
configure_file(${CMAKE_SOURCE_DIR}/include/config.h.in ${CMAKE_SOURCE_DIR}/include/config.h)



# Add an executable with the above sources
add_executable("SM_${NS}x${NT}" ${SRC_FILES})

# Set the directories that should be included in the build command for this target
# when running g++ these will be included as -I/directory/path/
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
target_link_libraries("SM_${NS}x${NT}" PUBLIC MPI::MPI_CXX)
#find_package(OpenMP REQUIRED)
#target_link_libraries("SM_${NS}x${NT}" PUBLIC OpenMP::OpenMP_CXX)


#PROJECT_SOURCE_DIR: it indicates where src is located
#cmake -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER="C:/msys64/ucrt64/bin/g++.exe" -DCMAKE_C_COMPILER="C:/msys64/ucrt64/bin/gcc.exe" ../


