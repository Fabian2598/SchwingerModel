# SchwingerModel

This branch only implements parallelization in 1D. It was useful for learning purposes, I don't recommend using it.

To compile create a new folder 

```
mkdir build
```

## Linux
In the `build` folder run the following commands:

```
cmake ../
cmake --build .
```

This will create an executable for you to run. The **lattice dimensions** are fixed in the **CMakeLists.txt**.
You can change the dimensions there as well as the executable name.

A running example with HMC is shown below

```
./SM_NSxNT.exe
----------------------------
|  Two-flavor Schwinger model   |
| Hybrid Monte Carlo simulation |
----------------------------
Ns NS Nt NT
m0 min: 0
m0 max: 0
Number of masses in [m0_min, m0_max] 1
Molecular dynamics steps: 8
Trajectory length: 1
beta: 2
Thermalization: 500
Measurements: 1000
Step (sweeps between measurements): 10
Save configurations yes/no (1 or 0): 1
```

