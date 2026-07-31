#!/bin/bash

#Script has to be run in the same directory where CMakeLists.txt is located
NX=64 #lattice dimensions
NT=64 
SOURCEFILE="readBinConf.cpp"
CONF_PATH="2D_U1_Ns64_Nt64_b40000_m02000_0.ctxt" #Change this accordingly
NAME="human_readable_conf.txt"
sed -i "s/constexpr int Nx=.*/constexpr int Nx= ${NX};/" "$SOURCEFILE"
sed -i "s/constexpr int Nt=.*/constexpr int Nt= ${NT};/" "$SOURCEFILE"
g++ ${SOURCEFILE} -O3 -o readBinConf
printf "${CONF_PATH}\n${NAME}" >> filenames
./readBinConf < filenames
rm filenames
rm readBinConf