# SchwingerModel

MPI parallelized simulation of the two-flavor Schwinger Model with degenerate fermions. The simulation is performed using pseudofermions, HMC and conjugate gradient to invert $(DD^\dagger)^{-1}$. Check HMC_doc.pdf for a detailed explanation of the Hybrid Monte Carlo in this context. 

The code implements Wilson fermions 

$$S[\Phi,\Phi^\dagger,U_\mu]=\beta\sum_{x\in V}\sum_{\mu<\nu}\textrm{Re}\left(1-U_{\mu\nu}(x)\right)+\Phi_\alpha^\dagger(\textbf{n}') (DD^\dagger)_{\textbf{n'},\textbf{n}}^{-1\,\alpha\beta}\Phi_{\beta} (\textbf{n})$$

where $U_{\mu\nu}$ is the plaquette and $\Phi$ represents the pseudofermion field. The Dirac operator has the following structure 

$$D\left[\textbf{n}',\textbf{n}\right]^{\alpha\beta} = \left(m_0 + 2\right)\delta^{\alpha\beta}\delta_{\textbf{n}',\textbf{n}} - \frac{1}{2} \sum_{\mu=\{0,1\}}
	\left[
		\left(1-\sigma_{\mu}\right)^{\alpha\beta} U_{\mu}(\textbf{n}') \delta_{\textbf{n}' + \hat{\mu},\textbf{n} }
		+\left(1+\sigma_{\mu}\right)^{\alpha\beta} U_{\mu}^\dagger(\textbf{n}'-\hat{\mu}) \delta_{\textbf{n}' - \hat{\mu},\textbf{n} }
	\right].$$

$$
\sigma_0 = \begin{pmatrix} 0 & 1 \\ 
1 & 0 \end{pmatrix}, \quad \sigma_1 = \begin{pmatrix}0 & -i \\ 
i & 0\end{pmatrix}.
$$ 

$\mu=0$ refers to time and $\mu=1$ to space. 


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

This will create an executable `SM_NSxNT`. The **lattice dimensions** are fixed in the **CMakeLists.txt**.
You can change the dimensions there as well as the executable name.

The simulation is executed as follows:

```
mpirun -n #cores SM_NSxNT
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

* `ranks_x` (`ranks_t`) is the number of cores assigned to the $x$-direction ($t$-direction). The total number of cores (#cores) is `ranks_x` x `ranks_t`. The program requires that $N_x$ ($N_t$) is exactly divisible by `ranks_x` (`ranks_t`). This is done to balance the workload, since each core deals with exactly $N_x\cdot N_t/$(ranks_x ranks_t) lattice sites. 
* $m_0$ is the bare mass parameter. The critical values were determined in [N. Christian, K. Jansen, K. Nagai and B. Pollakowski. “Scaling test of the fermion actions in the Schwinger Model”, Nucl. Phys. B, 739, (2006)] for different values of $\beta$ (see the table below). For $m_0<-m_{crit}$ the simulations are non-physical. Close to $m_{crit}$ the Dirac operator becomes highly ill-conditioned. 


| $\beta$ | $-m_{crit}$ |
| :-----: | :---------: |
|    1    |  0.3204(7)  |
|    2    |  0.1968(9)  |
|    3    |  0.1351(2)  |
|    4    |  0.1033(1)  |
|    5    |  0.0840(1)  |
|    6    |  0.0719(1)  |

* Molecular dynamics steps: Number of steps for the leapgrog integrator (integer number).
* Trajectory length: length of the integration trajectory (double number).

Both the MD steps and the trajectory length have to be manually tuned to obtain an acceptance rate between $0.6 - 0.8$, which indicates a good decorrelation between configurations. For most cases, $\tau=1.0$ and MD steps = 10 works fine. Close to the critical mass, these numbers have to be carefully tuned by increasing the number of steps and decreasing the trajectory length. Increasing the steps will also increase the simulation time. 

* $\beta$: inverse gauge coupling. The continuum limit is obtained when $\beta\rightarrow \infty$.
* Thermalization: number of configurations discarded for *thermalizing* the system. 
* Measurements: number of configurations used for taking measurements. If save_conf = 1, this is the number of confs that are written to disk after thermalization. For save_conf = 0 no confs are stored, which is useful for testing purposes. 
* Steps (sweep between measurements): number of configurations discarded between each measurement. This is done to decorrelate the configurations.

Configurations are stored in a binary format. To convert them to human readable text compile and run `readBinConf.cpp`



We provide a simple bash script (`run.sh`) to portray the compilation and execution of the simulation. 

## Windows

The instructions are essentially the same. The CMakeLists.txt only needs the address of your C++ and C compiler on lines 4 and 5. 
Then, in the `build` folder, run the following commands:
```
cmake -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER=C:\msys64\ucrt64\bin\g++ -DCMAKE_C_COMPILER=C:\msys64\ucrt64\bin\gcc ../
```

 This command depends on the compiler you are using. In this case, we are using MinGW. If you are using another compiler 
 you have to change the `-G` flag. The `-DCMAKE_CXX_COMPILER` and `-DCMAKE_C_COMPILER` flags are the address of the compiler.

Then you can run the executable

```
SM_NSxNT.exe
```