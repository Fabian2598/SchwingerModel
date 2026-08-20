import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


def Datos(path):
    """
    Load a data file containing three columns: t, C(t), and dC(t).
    Returns an array with shape (3, N).
    """
    return np.loadtxt(path).T


def EffectiveMass(Corr, dCorr, Nt):
    """
    Calculate the effective mass and its error using linear
    error propagation from the correlator uncertainties.
    
    Note: The error propagation is not entirely correct, since the measurements 
    from the correlator could be autocorrelated. However, writing a routine to
    perform the fit and get the effective mass for each jackknife block 
    would be very time consuming.

    Parameters
    ----------
    Corr : array
        Correlator C(t).
    dCorr : array
        Error of the correlator dC(t), obtained from jackknife.
    Nt : int
        Temporal lattice extent.

    Returns
    -------
    meff : array
        Effective mass m_eff(t).
    dmeff : array
        Error of the effective mass.
    """

    def f(x, Corr, t):
        #We need the root of this function, which corresponds to meff
        return (Corr[t] / Corr[(t + 1) % Nt]- np.cosh(x * (t - Nt / 2))/ np.cosh(x * ((t + 1) - Nt / 2)))

    meff = np.full(Nt, np.nan)
    dmeff = np.full(Nt, np.nan)

    for t in range(Nt):
        tp1 = (t + 1) % Nt
        # Avoid division by zero
        if Corr[tp1] == 0:
            continue
        # Ratio C(t)/C(t+1)
        ratio = Corr[t] / Corr[tp1]
        # Solve for m_eff
        try:
            sol = optimize.root_scalar(f,args=(Corr, t),method='secant',x0=0.5,x1=1.0,xtol=1e-12,rtol=1e-12,maxiter=1000)
        except (ValueError, RuntimeError):
            print(f"Solution not found for t={t}")
            continue
        if not sol.converged:
            print(f"Solution not found for t={t}")
            continue

        m = sol.root
        meff[t] = m

        # ---------------------------------------------#
        # Error propagation
        # R = C(t)/C(t+1)
        # R(m) = cosh[m(t-Nt/2)] / cosh[m(t+1-Nt/2)]
        # dm/dR = 1 / (dR/dm)
        # ---------------------------------------------#
        a = t - Nt / 2
        b = (t + 1) - Nt / 2
        numerator = np.cosh(m * a)
        denominator = np.cosh(m * b)
        R = numerator / denominator
        # dR/dm
        dRdm = R * (a * np.tanh(m * a) - b * np.tanh(m * b))
        # If derivative is essentially zero, error propagation
        # becomes numerically unstable.
        if np.abs(dRdm) < 1e-14:
            dmeff[t] = np.nan
            continue
        # dR/dC(t)
        dR_dCt = 1.0 / Corr[tp1]
        # dR/dC(t+1)
        dR_dCtp1 = -Corr[t] / Corr[tp1]**2
        # Propagate C(t) and C(t+1) errors to R
        dR = np.sqrt((dR_dCt * dCorr[t])**2+(dR_dCtp1 * dCorr[tp1])**2)
        # Propagate R error to m
        dmeff[t] = dR / np.abs(dRdm)

    return meff, dmeff

def kappa2mass(kappa):
    #Convert hopping parameter to bare mass
    return 1 / (2 * kappa) - 2

def mass2kappa(m0):
    #Convert bare mass to hoping parameter
    return 1 / (2 * (m0 + 2))

def format_number(number):
    #Formatting needed for the corr.txt files
    s = f"{number:.4f}"
    s = s.replace('.', '')
    return s

def plot_correlators(correlator_dir,obs, masses, beta, Nx, Nt, save=False):
    """
    Parameters
    ----------
    correlator_dir: string 
                    directory with the .txt files with the correlators
    obs: string
            observable (pion or pcac)
    masses: list 
            bare masses to be analyzed
    """
    fig = plt.figure(dpi=100)
    txt_string = ""
    sizeFont = 15
    if obs=="pion":
        txt_string = ""
        plt.title(r'Correlation function for pion, $N_x$={0}, $N_t={1}$, $\beta$={2}'.format(Nx, Nt, beta),size=sizeFont)
        plt.ylabel(r'$|c(t)|$',size=sizeFont)
        plt.yscale('log')
    elif obs == "pcac":
        txt_string = "PCAC"
        plt.title(r"$m_{\mathrm{PCAC}}(t)$"+ r", $N_x$={0}, $N_t={1}$, $\beta$={2}".format(Nx, Nt, beta),size=sizeFont)
        plt.ylabel(r'$|m_{\mathrm{PCAC}}(t)|$',size=sizeFont)
    else: 
        print("Provide a valid observable, pion or pcac")
        return
    
    plt.xlabel(r"$t$",size=sizeFont)
    
    for m0 in masses:
        path = (correlator_dir+ "2D_U1_{0}x{1}_b{2}_m{3}_corr{4}.txt" .format(Nx, Nt, beta, format_number(m0),txt_string))
        t, Corr, dCorr = Datos(path)
        Corr = np.abs(Corr)
        dCorr = np.abs(dCorr)
        plt.errorbar(t,Corr,yerr=dCorr,fmt='*',markersize=5,elinewidth=0.5,solid_capstyle='projecting',
            label='$m_0=${0}'.format(np.round(m0, 4)),capsize=1.5)
    
    plt.legend()
    plt.show()

    if save:
        fig.savefig("correlators_{0}_b{1}_{2}x{3}.pdf".format(txt_string,beta, Nx, Nt))


def print_correlator(correlator_dir, obs,m0, beta, Nx, Nt):
    txt_string = ""
    if obs=="pion":
        txt_string = ""
    elif obs == "pcac":
        txt_string = "PCAC"
    else: 
        print("Provide a valid observable, pion or pcac")
        return
    
    path = (correlator_dir+ "2D_U1_{0}x{1}_b{2}_m{3}_corr{4}.txt" .format(Nx, Nt, beta, format_number(m0),txt_string))
    t, Corr, dCorr = Datos(path)
    for i in range(Nt):
        print('C({0}) = {1} +- {2}'.format(i, Corr[i], dCorr[i]))

def plot_effective_mass(correlator_dir,masses,beta,Nx,Nt,mean_ranges,save=False):
    sizeFont = 15
    fig = plt.figure(dpi=100)
    plt.title(r'Effective mass, $N_x$={0}, $N_t={1}$, $\beta$={2}'.format(Nx, Nt, beta),size=sizeFont)
    plt.ylabel(r'$m_{\mathrm{eff}}(t)$',size=sizeFont)
    plt.xlabel(r"$t$",size=sizeFont)
    t0, t1 = mean_ranges
    for m0 in masses:
        path = (correlator_dir+ "2D_U1_{0}x{1}_b{2}_m{3}_corr.txt".format(Nx, Nt, beta, format_number(m0)))
        t, Corr, dCorr = Datos(path)
        Corr = np.abs(Corr)
        dCorr = np.abs(dCorr)
        # Calculate m_eff and its propagated error
        meff, dmeff = EffectiveMass(Corr,dCorr,Nt)

        # Plot effective mass
        plt.errorbar(t,meff,yerr=dmeff,fmt='*',markersize=5,elinewidth=0.5,solid_capstyle='projecting',
            label='$m_0=${0}'.format(np.round(m0, 4)),capsize=1.5)
        
        # Plateau mass
        plateau_meff = meff[t0:t1]
        plateau_dmeff = dmeff[t0:t1]

        # Remove invalid values
        valid = (np.isfinite(plateau_meff) & np.isfinite(plateau_dmeff))

        plateau_meff = plateau_meff[valid]
        plateau_dmeff = plateau_dmeff[valid]

        if len(plateau_meff) == 0:
            print("No valid effective masses in plateau range "f"for m0={m0}")
            continue

        # Plateau mass
        mpi = np.mean(plateau_meff)

        # Error of the mean, assuming the individual m_eff(t)
        # values are statistically independent.
        dmpi = (np.sqrt(np.sum(plateau_dmeff**2))/ len(plateau_dmeff))

        print("m0 = {0}, kappa = {1}, mpi = {2} +- {3}".format(
            np.round(m0, 5),mass2kappa(m0),np.round(mpi, 4),np.round(dmpi, 4))
             )

    plt.legend()
    plt.show()
    if save:
        fig.savefig("meff_b{0}_{1}x{2}.pdf".format(beta, Nx, Nt))


def jackknife_mean(data):
    """
    Jackknife estimate and error of the mean of a 1D array.
    This version considers means of size len(data)-1 ...
    """
    data = np.asarray(data, dtype=float)
    if data.ndim != 1:
        raise ValueError("data must be a one-dimensional array.")
    N = len(data)
    if N < 2:
        raise ValueError("At least two measurements are required for jackknife.")
    total = np.sum(data)
    samples = (total - data) / (N - 1)
    mean = np.mean(samples)

    error = np.sqrt((N - 1) / N* np.sum((samples - mean)**2))
    return mean, error, samples

def compute_pcac(correlator_dir, masses, beta, Nx, Nt,mean_ranges):
    """
    Computes the PCAC mass from the data in the .txt file
    """
    t0, t1 = mean_ranges
    m, dm = [], []
    for m0 in masses:
        path = (correlator_dir+ "2D_U1_{0}x{1}_b{2}_m{3}_corr{4}.txt" .format(Nx, Nt, beta, format_number(m0),"PCAC"))
        t, mpcac, dpcac = Datos(path)
        mpcac = np.abs(mpcac)
        _,error,_ = jackknife_mean(mpcac[t0:t1])
        print("mpcac = {0} +- {1}".format(np.round(np.mean(mpcac[t0:t1]),4),np.round(error,4)))
        m.append(np.mean(mpcac[t0:t1]))
        dm.append(error)
    m, dm = np.array(m), np.array(dm)
    return m, dm
    

def compute_mpi(correlator_dir,masses,beta,Nx,Nt,mean_ranges):
    """
    Computes the pion mass for different bare masses
    """
    t0, t1 = mean_ranges
    Mpi, dMpi = [], []
    for m0 in masses:
        path = (correlator_dir+ "2D_U1_{0}x{1}_b{2}_m{3}_corr.txt".format(Nx, Nt, beta, format_number(m0)))
        t, Corr, dCorr = Datos(path)
        Corr = np.abs(Corr)
        dCorr = np.abs(dCorr)
        # Calculate m_eff and its propagated error
        meff, dmeff = EffectiveMass(Corr,dCorr,Nt)
        # Plateau mass
        plateau_meff = meff[t0:t1]
        plateau_dmeff = dmeff[t0:t1]
        # Remove invalid values
        valid = (np.isfinite(plateau_meff) & np.isfinite(plateau_dmeff))
        plateau_meff = plateau_meff[valid]
        plateau_dmeff = plateau_dmeff[valid]
        if len(plateau_meff) == 0:
            print("No valid effective masses in plateau range "f"for m0={m0}")
            continue
        # Plateau mass
        mpi = np.mean(plateau_meff)
        dmpi = (np.sqrt(np.sum(plateau_dmeff**2))/ len(plateau_dmeff))
        Mpi.append(mpi)
        dMpi.append(dmpi)
    Mpi = np.array(Mpi)
    dMpi = np.array(dMpi)
    return Mpi, dMpi

def mpi_vs_mpcac(correlator_dir,masses,beta,Nx,Nt,mean_ranges):
    t0, t1 = mean_ranges
    Mpi, dMpi = compute_mpi(correlator_dir,masses,beta,Nx,Nt,mean_ranges)
    Mpcac, dMpcac = compute_pcac(correlator_dir, masses, beta, Nx, Nt,mean_ranges)
   
    fig = plt.figure(dpi=100)
    plt.title(r'$N_x$={0}, $N_t={1}$, $\beta$={2}'.format(Nx, Nt, beta),size=15)
    plt.ylabel(r'$m_{\pi}$',size=15)
    plt.xlabel(r"$m_{\mathrm{PCAC}}$",size=15)
    x0, x1 = 0, 0.16
    plt.xlim([x0,x1])
    plt.ylim([0,0.5])
    plt.errorbar(Mpcac,Mpi,xerr=dMpcac,yerr=dMpi,fmt='*',markersize=5,elinewidth=0.5,solid_capstyle='projecting',capsize=1.5,label='HMC simulation')

    #Prediction by Smilga
    x = np.linspace(x0,x1,500)
    g = 1/np.sqrt(beta)
    y = 2.008*(x**2*g)**(1/3)
    plt.plot(x,y,label=r"Smilga prediction $m_\pi=2.008 \left( m^2 g\right)^{1/3}$")
    plt.legend()