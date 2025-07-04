cmake_minimum_required(VERSION 3.26.5)
set(CMAKE_CXX_STANDARD 20) 
set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(CMAKE_C_COMPILER "C:/msys64/ucrt64/bin/gcc.exe")
#set(CMAKE_CXX_COMPILER "C:/msys64/ucrt64/bin/g++.exe")
project(Schwinger_model VERSION 1.0)


set(SOURCES
    src/statistics.cpp
    src/variables.cpp
    src/main.cpp   
    src/gauge_conf.cpp
    src/dirac_operator.cpp
    src/conjugate_gradient.cpp
    src/hmc.cpp
)

set(NS "16") 
set(NT "16") 

configure_file(${CMAKE_SOURCE_DIR}/include/config.h.in ${CMAKE_SOURCE_DIR}/include/config.h)


add_executable("SM_${NS}x${NT}" ${SOURCES})

# This is where the header files are located
target_include_directories("SM_${NS}x${NT}"
        PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)
target_compile_options("SM_${NS}x${NT}" PRIVATE
			-O3
			#-fopenmp
)

#-----------------mass-----------------#
set(SOURCES_MASS
    src/statistics.cpp
    src/variables.cpp
    src/mass.cpp   
    src/gauge_conf.cpp
    src/dirac_operator.cpp
    src/conjugate_gradient.cpp
)

add_executable("mass_${NS}x${NT}" ${SOURCES_MASS})

target_include_directories("mass_${NS}x${NT}"
        PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)
target_compile_options("mass_${NS}x${NT}" PRIVATE
			-O3
)
#---------------------------------------#


#find_package(OpenMP REQUIRED)
#target_link_libraries("SM_${NS}x${NT}" PUBLIC OpenMP::OpenMP_CXX)


#PROJECT_SOURCE_DIR: it indicates where src is located
#cmake -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER="C:/msys64/ucrt64/bin/g++.exe" -DCMAKE_C_COMPILER="C:/msys64/ucrt64/bin/gcc.exe" ../
#dir/s/b *.txt > confs.txt

