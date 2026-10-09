#!/bin/bash
#
# Exact test of the four-flavor HMC on a 4x4 lattice. Since det(D)^4 = det(D)^2 * det(D)^2,
#     <O>_{Nf=4} = <O det(D)^2>_{Nf=2} / <det(D)^2>_{Nf=2},
# with det(D) computed from the dense Dirac matrix. The script builds the two-flavor and the four-flavor
# Wilson versions in a temporary copy of the sources (so the lattice size of the main build is not touched),
# runs independent HMC chains with both, and compares the direct four-flavor averages with the reweighted
# two-flavor averages (rw_analysis.py, requires Python 3 with NumPy).
#
# Usage (from any directory):
#     tests/reweighting/run_reweighting_test.sh
# Optional environment variables:
#     NCHAINS   independent chains per build and parameter set (default 4)
#     NTRAJ     trajectories per chain after thermalization (default 50000)
#     NPROC     chains running at the same time (default: number of cores)
#     WORKDIR   directory for the build and the data (default: a new temporary directory)
#     LAUNCHER  command used to start each chain, e.g. "mpirun -n 1" or "srun -n 1" (default: none)
#
set -euo pipefail

NCHAINS=${NCHAINS:-4}
NTRAJ=${NTRAJ:-50000}
NPROC=${NPROC:-$(nproc)}
LAUNCHER=${LAUNCHER:-}
NTHERM=1000
MD_STEPS=10
TAU=1.0
L=4                              # lattice size (L x L)
PARAMS=("2.0 0.0" "1.0 0.1")     # (beta, m0) pairs

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
WORKDIR=${WORKDIR:-$(mktemp -d)}
mkdir -p "${WORKDIR}/data"

echo "Copying the sources to ${WORKDIR}"
cp -r "${REPO}/src" "${REPO}/include" "${REPO}/tests" "${REPO}/CMakeLists.txt" "${WORKDIR}/"
sed -i "s/^set(NS \"[0-9]*\")/set(NS \"${L}\")/" "${WORKDIR}/CMakeLists.txt"
sed -i "s/^set(NT \"[0-9]*\")/set(NT \"${L}\")/" "${WORKDIR}/CMakeLists.txt"
cat >> "${WORKDIR}/CMakeLists.txt" <<'EOF'

#Reweighting test (added by tests/reweighting/run_reweighting_test.sh)
add_executable(rw_chain tests/reweighting/rw_chain.cpp)
target_link_libraries(rw_chain PRIVATE common_lib)
EOF

for NF in 2 4; do
    FLAG=OFF
    [ "${NF}" = 4 ] && FLAG=ON
    echo "Building Nf=${NF} (${L}x${L} lattice)"
    LOG="${WORKDIR}/cmake_nf${NF}.log"
    { cmake -S "${WORKDIR}" -B "${WORKDIR}/build_nf${NF}" -DFOUR_FLAVORS=${FLAG} -DCLOVER=OFF -DTWISTED_MASS=OFF \
            -DBUILD_TESTS=OFF && cmake --build "${WORKDIR}/build_nf${NF}" --target rw_chain -j "${NPROC}"; } > "${LOG}" 2>&1 \
        || { echo "Build failed, see ${LOG}"; exit 1; }
done

JOBS="${WORKDIR}/jobs.txt"
: > "${JOBS}"
for P in "${PARAMS[@]}"; do
    read -r BETA M0 <<< "${P}"
    for NF in 2 4; do
        for C in $(seq 1 "${NCHAINS}"); do
            echo "${LAUNCHER} ${WORKDIR}/build_nf${NF}/rw_chain ${BETA} ${M0} ${NTHERM} ${NTRAJ} ${MD_STEPS} ${TAU} ${WORKDIR}/data/nf${NF}_b${BETA}_m${M0}_c${C}.dat" >> "${JOBS}"
        done
    done
done
echo "Running $(wc -l < "${JOBS}") chains of ${NTRAJ} trajectories, ${NPROC} at a time"
xargs -P "${NPROC}" -I{} sh -c '{}' < "${JOBS}"

python3 "${REPO}/tests/reweighting/rw_analysis.py" "${WORKDIR}/data"
