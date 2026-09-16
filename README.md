# HMC for the Schwinger model

This project implements an MPI-parallel Monte Carlo simulation of the two-flavor Schwinger model using a Hybrid Monte Carlo (HMC) algorithm. The simulation uses Wilson fermions, pseudofermions, and iterative solvers for the fermion operator. For an OpenMP implementation, check the OpenMP branch.

The number of leapfrog steps is tuned automatically during thermalization to reach a target acceptance rate, at fixed trajectory length. See the "Automatic step-size tuning" section below.

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

### Twisted mass

Twisted mass is also implemented in the code. To compile the twisted mass version build:

```bash
cmake -S . -B build -DTWISTED_MASS=ON
cmake --build build
```

### Clover term

The clover (Sheikholeslami-Wohlert) term is also implemented. To compile with the clover term enabled:

```bash
cmake -S . -B build -DCLOVER=ON
cmake --build build
```

Compiling with `-DCLOVER=ON` adds an extra runtime prompt for `csw`, the Sheikholeslami-Wohlert constant (see the parameter list below). `-DCLOVER=ON` and `-DTWISTED_MASS=ON` can be combined:

```bash
cmake -S . -B build -DCLOVER=ON -DTWISTED_MASS=ON
cmake --build build
```

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

In case the program was compiled for the twisted mass operator (`-DTWISTED_MASS=ON`), an extra parameter (`mu0`, prompted as "mu (twisted mass)") will be requested. In case it was compiled with the clover term (`-DCLOVER=ON`), an extra parameter (`csw`) will be requested.

### Parameter descriptions

- `ranks_x` and `ranks_t`: number of MPI ranks in the $x$ and $t$ directions. The total number of processes is `ranks_x * ranks_t`.
- The lattice dimensions must be divisible by the corresponding rank count, so the workload is balanced across processes.
- `m0`: bare mass parameter.
- `mu0` (only if compiled with `-DTWISTED_MASS=ON`): twisted mass parameter. Prompted as "mu (twisted mass)".
- `csw` (only if compiled with `-DCLOVER=ON`): Sheikholeslami-Wohlert constant for the clover term.
- `Molecular dynamics steps`: initial number of leapfrog integration steps. This value is only a starting guess: it is adjusted automatically during thermalization (see below), and the tuned value is the one actually used for the measurement phase.
- `Trajectory length`: integration length in lattice units. Unlike the number of molecular dynamics steps, this value is kept fixed throughout the run.
- `beta`: inverse gauge coupling.
- `Thermalization`: number of configurations discarded before measurements begin. This also sets how many trajectories are available for tuning the number of molecular dynamics steps, so it should not be too small (a few hundred at least) for the tuning to converge.
- `Measurements`: number of configurations used for measurements.
- `Step`: number of sweeps discarded between saved measurements.
- `Save configurations`: set to `1` to write configurations to disk, or `0` to skip writing them.

### Automatic step-size tuning

At fixed trajectory length, the number of molecular dynamics steps controls the acceptance rate. Rather than choosing it by hand for every set of parameters, the code tunes it automatically during thermalization, targeting an acceptance rate of 0.78. The tuning uses dual averaging and runs in three phases within the thermalization loop: an initial burn-in with adaptation enabled, so that the chain can move away from the random initial configuration; a reset of the averaging statistics once the chain is no longer at its initial hot start; and a final adaptation phase before the number of molecular dynamics steps is frozen for the remaining thermalization and measurement trajectories. The implementation is in `include/hmc_tuner.h`.

Because the tuner needs a reasonable number of trajectories to converge, `Thermalization` should be at least a few hundred for the tuned acceptance rate to be close to the target. With a very small `Thermalization` (for example, in a quick test run) the reported number of molecular dynamics steps may not have converged.

The script `run.sh` edits the lattice dimensions, configures the build if needed, builds the HMC executable, supplies a sample parameter set, and runs it. Review its variables before use. It currently moves the executable into the repository root and writes output there.

## Critical mass values

The bare mass parameter must remain above the critical mass to avoid unphysical configurations. The following values are useful guides. They are valid for `csw = 0` (i.e. without the clover term, or with `-DCLOVER=ON` and `csw` set to 0); the critical mass shifts when the clover term is switched on with `csw != 0`.

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

Once the simulation is completed and the configurations are written to disk, the `mass_NSxNT` program determines the necessary correlators, which are later analyzed with Python (check mass_analysis.ipynb) to compute $m_\pi$ and $m_{\mathrm{PCAC}}$. When executed, the program will ask for some parameters 

```text
--------------------------------------------
|  Pion and PCAC correlators computation   |
--------------------------------------------
Nx NS Nt NT
ranks_x: number of processes on the x direction
ranks_t: number of processes on the t direction
m0: same as the simulation
beta: same as the simulation (using for file naming)
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

The number of molecular dynamics steps is tuned automatically during thermalization (see "Automatic step-size tuning" above), so this should not normally require manual intervention. If the acceptance rate at the end of thermalization is still far from the target, increase `Thermalization` so that the tuner has more trajectories to converge, or check that the trajectory length is reasonable (values much shorter than 1 lattice unit lead to poor decorrelation between measurements even once the acceptance rate is on target).

### Missing output

Set `Save configurations` to `1` and verify that the process can write to the working directory.

## References

For background on the physics and algorithm, see [HMC_doc.pdf](HMC_doc.pdf) and the scaling-test reference cited above. Check [results.pdf](mass_analysis/results.pdf) for some results.