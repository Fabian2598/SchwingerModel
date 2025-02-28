cmake_minimum_required(VERSION 3.26.5)
set(CMAKE_CXX_STANDARD 11) 
set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(CMAKE_C_COMPILER "C:/msys64/ucrt64/bin/gcc.exe")
#set(CMAKE_CXX_COMPILER "C:/msys64/ucrt64/bin/g++.exe")
project(Schwinger_model VERSION 1.0)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    src/variables.cpp
    src/gauge_conf.cpp
    src/matrix_operations.cpp
    src/conjugate_gradient.cpp
    src/amg.cpp
    src/main.cpp  
)

#Lattice dimensions
set(NS "32") 
set(NT "32") 
#Number of lattice blocks (used for aggregation)
#N_i / BLOCK_i has to be an integer number
set(BLOCK_X "8")
set(BLOCK_T "8")
#Number of test vectors
set(NTEST "15") 
configure_file(${CMAKE_SOURCE_DIR}/include/config.h.in ${CMAKE_SOURCE_DIR}/include/config.h)



# Add an executable with the above sources
add_executable("SM_${NS}x${NT}" ${SOURCES} "include/matrix_operations.h" 
 "include/variables.h" 
 "include/gauge_conf.h" 
 "include/conjugate_gradient.h"
 "include/amg.h"
 "include/operator_overloads.h" 
)

# Set the directories that should be included in the build command for this target
# when running g++ these will be included as -I/directory/path/
target_include_directories("SM_${NS}x${NT}"
        PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)
target_compile_options("SM_${NS}x${NT}" PRIVATE
			#-O3
			#-fopenmp
)

#find_package(OpenMP REQUIRED)
#target_link_libraries("SM_${NS}x${NT}" PUBLIC OpenMP::OpenMP_CXX)


#PROJECT_SOURCE_DIR: it indicates where src is located
#cmake -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER="C:/msys64/ucrt64/bin/g++.exe" -DCMAKE_C_COMPILER="C:/msys64/ucrt64/bin/gcc.exe" ../


