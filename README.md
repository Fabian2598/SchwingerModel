# HMC for the Schwinger model

**Note:** The CUDA implementation was written with AI and has not been tested by a human being yet. 

This project implements an MPI-parallel Monte Carlo simulation of the two-flavor Schwinger model using a Hybrid Monte Carlo (HMC) algorithm. The simulation uses Wilson fermions, pseudofermions, and iterative solvers for the fermion operator. For an OpenMP implementation, check the OpenMP branch.

The loops over the lattice sites of every MPI rank can run on NVIDIA GPUs with CUDA (one GPU per rank) or on the CPU. The backend is chosen at compile time, see "GPU (CUDA)" below.

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
- For the GPU version only: the CUDA Toolkit 12 or newer (the code uses C++20) and an NVIDIA GPU supported by it
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

### GPU (CUDA)

To run the lattice loops on NVIDIA GPUs, compile with `-DUSE_CUDA=ON` (it can be combined with `-DCLOVER` and `-DTWISTED_MASS`):

```bash
cmake -S . -B build -DUSE_CUDA=ON
cmake --build build
```

By default the code is compiled for the GPU of the machine where `cmake` runs. On login nodes without a GPU, give the architecture explicitly, for example `-DCMAKE_CUDA_ARCHITECTURES=80` for A100 or `-DCMAKE_CUDA_ARCHITECTURES=90` for H100/GH200.

All the lattice fields live in GPU memory during the whole simulation. Only the halos are communicated between ranks: they are packed into a contiguous buffer on the GPU, copied to the CPU and sent with MPI. If the MPI library is CUDA-aware, add `-DCUDA_AWARE_MPI=ON` and the GPU buffers are given directly to MPI, which skips the copies. With a single rank the halos are filled on the GPU and MPI is not involved. Random numbers (momenta, pseudofermions, initial configuration) are still generated on the CPU and copied to the GPU, so the GPU and CPU versions produce the same Markov chain for the same seeds, up to rounding.

