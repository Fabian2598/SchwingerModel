# SchwingerModel

Simulation of the two-flavor Schwinger Model with degenerate fermions. The simulation is performed using pseudofermions, HMC and (for the moment) conjugate gradient to invert $(DD^\dagger)^{-1}$. Later multigrid will replace conjugate gradient, so temporarily ignore the AMG directory. Check HMC_doc.pdf for a detailed explanation of the Hybrid Monte Carlo in this context.

The code implements Wilson fermions 

$$S[\Phi,\Phi^\dagger,U_\mu]=\beta\sum_{x\in V}\sum_{\mu<\nu}\textrm{Re}\left(1-U_{\mu\nu}(x)\right)+\Phi_\alpha^\dagger(\textbf{n}') (DD^\dagger)_{\textbf{n'},\textbf{n}}^{-1\,\alpha\beta}\Phi{\beta} (\textbf{n})$$

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

The other executable

 ```
mass_NSxNT.exe
```

measures the pion mass correlator, given a set of gauge configurations. The latter are generated during the simulation.
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

**Only the average plaquette value is measured, one can implement other observables.**
