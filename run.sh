#!/bin/bash

#Script has to be run in the same directory where CMakeLists.txt is located
NX=32 #lattice dimensions
NT=32 
RANKS_X=4 #This has to exactly divide NX
RANKS_T=2 #The same
RANKS=$((RANKS_X*$RANKS_T)) #Total number of cores
CMAKELISTS="CMakeLists.txt"
M0=0 #bare mass
BETA=2 #beta
MD_STEPS=10 #Moleculear dynamics steps
TAU=1.0 #Trajectory length
NTHERM=1000 #Thermalization
NMEAS=1000 #Measurements
NSTEPS=10 #Decorrelation steps between measurements
SAVE=0


sed -i "17s/set(NS \".*\")/set(NS \"${NX}\")/" "$CMAKELISTS"
sed -i "18s/set(NT \".*\")/set(NT \"${NT}\")/" "$CMAKELISTS"

DIR="build"
if [ ! -d "$DIR" ]; then
    mkdir build
    cd build 
    cmake ../
    cd ..
fi

echo "--------Compiling for Nx=${NX}, Nt=${NT}--------"
cd build
cmake --build .
mv SM_${NX}x${NT} ../
cd ../
printf "${RANKS_X}\n${RANKS_T}\n${M0}\n${MD_STEPS}\n${TAU}\n${BETA}\n${NTHERM}\n${NMEAS}\n${NSTEPS}\n${SAVE}" >> parameters
mpirun -n ${RANKS} SM_${NX}x${NT} < parameters
rm parameters