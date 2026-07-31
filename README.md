# HMC for the Schwinger model

This project implements an MPI-parallel Monte Carlo simulation of the two-flavor Schwinger model using a Hybrid Monte Carlo (HMC) algorithm. The simulation uses Wilson fermions, pseudofermions, and a conjugate-gradient solver to invert $(DD^\dagger)^{-1}$. For an OpenMP implementation check the OpenMP branch.

The code generates gauge configurations in a binary format, with the following structure:

$$
(x, t, \mu, \mathrm{Re}(U_\mu(t,x)), \mathrm{Im}(U_\mu(t,x)))
$$

The action implemented is

$$
S[\Phi,\Phi^\dagger,U_\mu]
= \beta \sum_{x\in V}\sum_{\mu<\nu}\mathrm{Re}\left(1-U_{\mu\nu}(x)\right)
+ \Phi_\alpha^\dagger(\mathbf{n}') (DD^\dagger)^{-1\,\alpha\beta}_{\mathbf{n}',\mathbf{n}} \Phi_\beta(\mathbf{n})
$$

where $U_{\mu\nu}$ is the plaquette and $\Phi$ is the pseudofermion field. The Dirac operator has the form

$$
D[\mathbf{n}',\mathbf{n}]^{\alpha\beta}
= (m_0 + 2)\delta^{\alpha\beta}\delta_{\mathbf{n}',\mathbf{n}}
- \frac{1}{2}\sum_{\mu=\{0,1\}}
\left[
\left(1-\sigma_\mu\right)^{\alpha\beta}U_\mu(\mathbf{n}')\delta_{\mathbf{n}' + \hat{\mu},\mathbf{n}}
+\left(1+\sigma_\mu\right)^{\alpha\beta}U_\mu^\dagger(\mathbf{n}'-\hat{\mu})\delta_{\mathbf{n}' - \hat{\mu},\mathbf{n}}
\right].
$$

with

$$
\sigma_0 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix},
\qquad
\sigma_1 = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}.
$$

Here $\mu=0$ denotes the time direction and $\mu=1$ denotes the spatial direction.

---

## Requirements

- CMake
- C++ compiler
- MPI implementation (OpenMPI, MPICH, etc.)
- Linux or Windows with a compatible toolchain

---

## Build the project

Create a build directory:

```bash
mkdir build
cd build
```

Then configure and build:

```bash
cmake ..
cmake --build .
```

This produces an executable named `SM_NSxNT` (or a name defined in `CMakeLists.txt`). The lattice dimensions are fixed in `CMakeLists.txt`, so you should edit them there before building.

---

## Run the simulation

Execute the program with MPI:

```bash
mpirun -n <number_of_cores> ./SM_NSxNT
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

### Notes on tuning

The molecular dynamics steps and trajectory length must be tuned to achieve a good acceptance rate, typically between $0.6$ and $0.8$. This indicates a reasonable level of decorrelation between configurations.

For many simulations, a reasonable starting point is:

- `MD steps = 10`
- `trajectory length = 1.0`

Close to the critical mass, longer runs and careful tuning are often required. Increasing the MD steps typically increases the cost of each trajectory.

---

## Critical mass values

The bare mass parameter must remain above the critical mass to avoid unphysical configurations. The critical values for various $\beta$ are:

| $\beta$ | $-m_{\mathrm{crit}}$ |
| :-----: | :------------------: |
|    1    |       0.3204(7)      |
|    2    |       0.1968(9)      |
|    3    |       0.1351(2)      |
|    4    |       0.1033(1)      |
|    5    |       0.0840(1)      |
|    6    |       0.0719(1)      |

These values are from the literature and are useful as a guide when choosing a physical mass parameter. See N. Christian, K. Jansen, K. Nagai and B. Pollakowski. “Scaling test of the fermion actions in the Schwinger Model”, Nucl. Phys. B, 739, (2006).

---

## Output files

Configurations are stored in binary format. After a successful run, the simulation produces gauge configurations and a summary file such as `_SimData.txt`.

To convert the binary gauge configurations to a human-readable text format, compile and run `readBinConf.cpp`. You may need to update the lattice dimensions inside that source file before conversion.

A simple convenience script is provided:

```bash
./run.sh
```

This script shows how the program can be compiled and executed in practice.

---

## Converting binary configurations to text

The repository includes a helper script and a C++ converter:

- `readBinConf.cpp`
- `readBin.sh`

The script compiles the converter and reads a binary configuration file. The output is written in the form:

```text
x, t, mu, Re(U_mu), Im(U_mu)
```

This is useful for inspecting or post-processing saved gauge configurations.

---

## Windows

The build steps are similar on Windows. The main difference is that the CMake configuration must point to the correct C and C++ compiler paths.

For example, with MinGW:

```bash
cmake -G "MinGW Makefiles" \
  -DCMAKE_CXX_COMPILER=C:\msys64\ucrt64\bin\g++ \
  -DCMAKE_C_COMPILER=C:\msys64\ucrt64\bin\gcc \
  ..
```

Then build the project and run the generated executable:

```bash
SM_NSxNT.exe
```

The exact compiler and generator flags may vary depending on your setup.

---

## Repository structure

```text
.
├── CMakeLists.txt
├── README.md
├── run.sh
├── readBin.sh
├── readBinConf.cpp
├── include/
│   ├── config.h
│   ├── conjugate_gradient.h
│   ├── dirac_operator.h
│   ├── gauge_conf.h
│   ├── hmc.h
│   ├── mpi_setup.h
│   ├── statistics.h
│   └── variables.h
├── src/
│   ├── conjugate_gradient.cpp
│   ├── dirac_operator.cpp
│   ├── gauge_conf.cpp
│   ├── hmc.cpp
│   ├── main.cpp
│   ├── statistics.cpp
│   └── variables.cpp
├── build/
└── data files and output logs
```

---

## Troubleshooting

### MPI errors

If CMake or the runtime reports missing MPI libraries, check that MPI is installed and that your compiler environment is configured properly.

### Wrong lattice dimensions

If the simulation or conversion code does not match the configuration files, verify the lattice dimensions in `CMakeLists.txt` and in `readBinConf.cpp`.

### Acceptance rate is too low

If the acceptance rate is poor, reduce the trajectory length or increase the number of molecular dynamics integration steps.

### Missing or corrupted output

Check the save flag in the simulation input and confirm that the process has write permissions in the working directory.

---

## References

For more background on the physics and the algorithm, see the project documentation and the cited work on scaling tests of fermion actions in the Schwinger model.
