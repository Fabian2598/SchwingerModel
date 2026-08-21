# HMC for the Schwinger model

This project implements an MPI-parallel Monte Carlo simulation of the two-flavor Schwinger model using a Hybrid Monte Carlo (HMC) algorithm. The simulation uses Wilson fermions, pseudofermions, and iterative solvers for the fermion operator. For an OpenMP implementation, check the OpenMP branch.

Gauge configurations are written in binary format with records of the form:

$$
(x, t, \mu, \mathrm{Re}(U_\mu(t,x)), \mathrm{Im}(U_\mu(t,x)))
$$

Here $\mu=0$ denotes the time direction and $\mu=1$ denotes the spatial direction. See [HMC_doc.pdf](HMC_doc.pdf) for details of the HMC formulation.

The `mass_analysis` directory contains Python helpers and a notebook for computing $m_\pi$ and $m_{\mathrm{PCAC}}$ from the correlators. Results are in the same directory.

## Requirements

- CMake
- A C++20 compiler
- An MPI implementation such as OpenMPI, MPICH, or Microsoft MPI
- Python 3 with Jupyter, NumPy, Matplotlib, and SciPy for the mass analysis
- Linux, Windows with a compatible toolchain, or another Unix-like environment

The shell scripts require Bash and common Unix utilities such as `sed` and `mv`. On Windows, run them from an MSYS2/MinGW environment or follow the manual CMake commands below.

## Build the project

From the repository root, configure and build the project:

```bash
cmake -S . -B build
cmake --build build
```

The lattice dimensions are set in `CMakeLists.txt`:

```cmake
set(NS "64")
set(NT "64")
```

Change `NS` and `NT` before configuring if a different lattice is required. The build produces `SM_${NS}x${NT}` for the HMC simulation and `mass_${NS}x${NT}` for correlator computation. On Windows, the executables have an `.exe` suffix.

## Run the simulation

For the default lattice, run the HMC executable with MPI:

```bash
mpirun -n <number-of-ranks> ./build/SM_64x64
```

The program will prompt for the simulation parameters. A typical example is:

```text
----------------------------
|  Two-flavor Schwinger model   |
| Hybrid Monte Carlo simulation |
----------------------------
Nx NS Nt NT
ranks_x: number of processes on the x direction
ranks_t: number of processes on the t direction
m0: 0
Molecular dynamics steps: 10
Trajectory length: 1
beta: 2
Thermalization: 1000
Measurements: 1000
Step (sweeps between measurements): 10
Save configurations yes/no (1 or 0): 1
```

### Parameter descriptions

- `ranks_x` and `ranks_t`: number of MPI ranks in the $x$ and $t$ directions. The total number of processes is `ranks_x * ranks_t`.
- The lattice dimensions must be divisible by the corresponding rank count, so the workload is balanced across processes.
- `m0`: bare mass parameter.
- `Molecular dynamics steps`: number of leapfrog integration steps.
- `Trajectory length`: integration length in lattice units.
- `beta`: inverse gauge coupling.
- `Thermalization`: number of configurations discarded before measurements begin.
- `Measurements`: number of configurations used for measurements.
- `Step`: number of sweeps discarded between saved measurements.
- `Save configurations`: set to `1` to write configurations to disk, or `0` to skip writing them.

The script `run.sh` edits the lattice dimensions, configures the build if needed, builds the HMC executable, supplies a sample parameter set, and runs it. Review its variables before use. It currently moves the executable into the repository root and writes output there.

## Critical mass values

The bare mass parameter must remain above the critical mass to avoid unphysical configurations. The following values are useful guides:

| $\beta$ | $-m_{\mathrm{crit}}$ |
| :-----: | :------------------: |
| 1 | 0.3204(7) |
| 2 | 0.1968(9) |
| 3 | 0.1351(2) |
| 4 | 0.1033(1) |
| 5 | 0.0840(1) |
| 6 | 0.0719(1) |

These values are from N. Christian, K. Jansen, K. Nagai and B. Pollakowski, “Scaling test of the fermion actions in the Schwinger Model”, *Nucl. Phys. B* 739 (2006).

## Output files

When `Save configurations` is `1`, the simulation writes binary `.ctxt` files and a summary file such as `_SimData.txt`.

## Convert binary configurations

`readBinConf.cpp` converts one binary configuration to text. Update `NX`, `NT`, and `CONF_PATH` in `readBin.sh`, then run it from the repository root:

```bash
./readBin.sh
```

The converter writes records in the form:

```text
x, t, mu, Re(U_mu), Im(U_mu)
```

## Windows (MSYS2/MinGW)

From an MSYS2 shell with MinGW and an MPI implementation available on `PATH`:

```bash
cmake -S . -B build -G "MinGW Makefiles" \
  -DCMAKE_CXX_COMPILER=C:/msys64/ucrt64/bin/g++.exe \
  -DCMAKE_C_COMPILER=C:/msys64/ucrt64/bin/gcc.exe
cmake --build build --config Release
mpiexec -n 4 ./build/SM_64x64.exe
```

The exact compiler, MPI implementation, and generator flags may vary. The Bash helper scripts are not native PowerShell scripts.

## Compute pion and PCAC masses

Once the simluation is completed and the configurations written to disk, the `mass_NSxNT` program determines the necessary correlators, which are later analyzed with Python (check mass_analysis.ipynb) to compute $m_\pi$ and $m_{\mathrm{PCAC}}$. When executed, the program will ask for some parameters 

```text
----------------------------
|  Pion correlator computation   |
----------------------------
Nx NS Nt NT
ranks_x: number of processes on the x direction
ranks_t: number of processes on the t direction
m0: same as the simulation
File with list of confs (ls -1 *.ctxt > confFiles.txt): confFiles.txt
```

The last one corresponds to a list with the name of the configurations that will be considered for the analysis. This can be easily created in Linux by typing 

```bash
ls -1 -v *.ctxt > confFiles.txt
```

in the directory with the configurations. The output of the program are two .txt files with the correlators for the pion and PCAC. The corresponding names are `2D_U1_${NX}x${NT}"_b${BETA}_m${M0}_corr.txt` and `2D_U1_${NX}x${NT}"_b${BETA}_m${M0}_corrPCAC.txt`. This files are used in the Python program to calculate $m_\pi$ and $m_{\mathrm{PCAC}}$.

## Troubleshooting

### MPI errors

If CMake or the runtime reports missing MPI libraries, verify that MPI is installed and that the compiler and MPI launcher are on `PATH`.

### Wrong lattice dimensions

Use the same `NS` and `NT` values when building the simulation, building the mass executable, and converting configurations. Also keep `readBinConf.cpp` consistent with the configuration dimensions.

### Low acceptance rate

Reduce the trajectory length or increase the number of molecular-dynamics steps.

### Missing output

Set `Save configurations` to `1` and verify that the process can write to the working directory.

## References

For background on the physics and algorithm, see [HMC_doc.pdf](HMC_doc.pdf) and the scaling-test reference cited above.