from settings import ERG,C,COSMO,MAXZ,KBOLTZMANN,NINTERP_Z,LITTLEH,HPLANCK
from settings import KPC,JY,PI
import numpy as np
#EMISSAGNFUNC=[None]
#EMISSSFGFUNC=[None]
#SPLINE_DICT={}
from settings import SPLINE_DICT
def sfr_rate(sfr,z,luminosity=False):
    '''
    star formation rate density (dN/dlogSFR dV) in Mpc^{-3} log SFR^{-1}
    detailed in equation 2 of Mancuso 2015 and Aversa 2015 eqs 9/10 and
    table 1.
    Args:
        sfr, float, (msolar/year)
        z, float, redshift
        luminosity,bool,just keep it false.
    Returns:
        sfr rate, float, dN/dlogSFR dV in h^3 Mpc^{-3} log SFR{-1}
    '''
    xi=np.log10(1+z)
    xi2=xi*xi
    xi3=xi*xi2
    phiz=10.**(-2.4 - 2.3*xi + 6.2*xi2 - 4.9*xi3)
    sfrz=10.**(1.09+3.2*xi-1.4*xi2-2.1*xi3)
    if luminosity:
        sfrz*=10.**9.8
    alphaz=1.2+0.5*xi-0.5*xi2+0.2*xi3
    omegaz=0.7-0.15*xi+0.16*xi2+0.01*xi3
    phisfr=phiz*(sfr/sfrz)**(1-alphaz)*np.exp(-(sfr/sfrz)**omegaz)
    return phisfr/LITTLEH**3.

def sfr_tot(z,luminosity=False):
    '''
    Total SFR rate density (msolar/year/Mpc^3)
    Args:
        z, float, redshift
        luminosity, bool, returns luminosity density instead of SFR
    Returns:
        float, total SFR rate density (msolar/year/Mpc^3)
    '''
    g=lambda x: sfr_rate(np.exp(x),z,luminosity)*np.exp(x)
    return integrate.quad(g,-100.,100.)[0]
def gaunt(nu,t=1e4,zi=1):
    '''
    Gaunt factor (equation 2 from Bonato 2017)
    Args: nu, frequency Hz
          t, temperature (kelvin)
          zi, proton number
    Returns gaunt factor (unitless)
    '''
    return np.log( np.exp(5.960-np.sqrt(3)/np.pi*np.log(zi*nu/1e9*(t/1e4)**-1.5))+np.exp(1.))
def lff(sfr,nu,t=1e4,zi=1):
    '''
    free-free radio luminosity as a function of star formation
    Args:
         sfr, star formation rate (msolar/yr)
         t, temperature (kelvin)
         zi, proton count
    Returns specific luminosity Watts/Hz
    '''
    return 3.75e26*(sfr)*(t/1e4)**0.3*gaunt(nu,t,zi)\
    *np.exp(-HPLANCK*nu/(KBOLTZMANN*t))*ERG
def lstarsync(sfr,nu):
    '''
    linear synchrotron function
    Args:
        sfr (msolar/year)
        t, temperature (kelvin)
        zi, proton count
    Returns specific luminosity Watts/Hz
    '''
    return 1.9e28*sfr*(nu/1e9)**-0.85*(1+(nu/20e9)**0.5)**-1*ERG
def lsync(sfr,nu,alpha=.866,beta=3):
    '''
    non-linear sychrotron function with power laws.
    Args:
        sfr, (msolar/year)
        t, temperature (kelvin)
        zi (proton count)
        alpha (coefficient)
        beta (power law in non-linear regime)
    Returns specific luminosity of synchrotron Watts/Hz.
    '''
    return lstarsync(1.,nu)/((alpha/sfr)**beta+alpha/sfr)
def lsyncpw(sfr,nu,alpha=.866,beta=3):
    '''
    piecewise parameterization of synchrotron emission
    sfr, (msolar/year)
    t, temperature (kelvin)
    zi, proton count
    alpha, coefficient
    beta, power law
    Returns specific luminosity Watts/Hz
    '''
    if type(sfr)==np.ndarray:
        selection=sfr<=alpha
        output=np.zeros_like(sfr)
        output[selection]=lstarsync(1,nu)*sfr[selection]**beta/alpha**(beta-1)
        output[np.invert(selection)]=lstarsync(sfr[np.invert(selection)],nu)
        return output
    else:
        if sfr<=alpha:
            return lstarsync(1,nu)*sfr**beta/alpha**(beta-1)
        else:
            return lstarsync(sfr,nu)
def sfr(l,nu,alpha=.866,beta=3.):
    '''
    Solve for SFR given luminosity and frequency
    Args:
        l, float, luminosity (Watts/Hz)
        nu, float, frequency (Hz)
    Returns:
        Solves for SFR given luminosity at frequency nu.
    '''
    lf=lff(1.,nu)
    ls=lstarsync(1.,nu)
    x0=l/(lf+ls)
    a1=np.sqrt(3.)*np.sqrt(27.*l*l*ls**4.*alpha**4.+4.*lf**3.*ls**3.*alpha**6.)
    b1=(9*l*ls*ls+a1)**(1./3)
    x1=-(2./3.)**(1./3.)*lf*alpha**2./b1+b1/(2.**(1./3.)*3.**(2./3.)*ls)
    if type(x0)==np.ndarray:
        return np.max([x0,x1],axis=1)
    else:
        return max(x0,x1)

def lfunc_sfg(l,nu,z):
    '''
    dN/dV/dlogL(freq)
    Args:
        l, luminosity (ergs/Hz)
        nu, observed frequency (Hz)
        z, redshift
    Returns:
        float luminosity function of star forming galaxies
        at redshift z and luminosity l in
        number Mpc^-3 h^3
    '''
    return sfr_rate(sfr(l,nu),z)
def dvc_dzdomega(z):
    '''
    comoving volume element
    dV_c/(d\Omega dz) at redshift z (Mpc^3/h^3)
    Args:
        z, redshift, float.
    Returns:
        float, comoving volume element (Mpc^3/h^3)
    '''
    return C*1e-3/COSMO.Hz(z)*LITTLEH\
    *COSMO.luminosityDistance(z)**2./(1.+z)**2.

def ztop(l,ztop0,dztop0,l0):
    '''
    helper funciton from massardi 2010
    Args:
        l, float, luminosity (1.4 GHz)
        ztop0, float, redshift of luminosity peak at l0
        dztop0, float, evolving luminosity peak redshift (with luminosity)
        l0, float, reference luminosity
    Returns:
        float, z of luminosity maximum.
    '''
    return ztop0+dztop0/(1+l0/l)
def lstar_agn(l,z,l0,kevo,ztop0,dztop0,mev):
    '''
    turnover luminosity from massardi 2010
    Args:
        l, float, luminosity density (1.4 GHz)
        z, float, redshift
        l0, float, reference luminosity
        dztop0, float, evolution of redshift with luminosity.
        mev, float, something
    Returns:
        turnover luminosity, lstar from Massardi 2010
    '''
    zt=ztop(l,ztop0,dztop0,l0)
    return l0*10**(2.*kevo*(z*zt-z**(1+mev)*zt**(1-mev)/(1.+mev)))
def lfunc_mass_agn(l,z,a,b,n0,l0,kevo,ztop0,dztop0,mev):
    '''
    1.4 GHz luminosity function of AGN from Massardi et al. 2010.
    in h^3 MPc^-3
    Args:
        l,float,luminosity (1.4 GHz)
        z,float,redshift
        a,.... floats, parameters in model
    Returns:
        luminosity function dN/dV/dlogL1.4GHz h^3 Mpc^-3
    '''
    ls=lstar_agn(l,z,l0,kevo,ztop0,dztop0,mev)
    return n0/((l/ls)**a+(l/ls)**b)/LITTLEH**3.#*dl0dlz(l,z,l0,kevo,ztop0,dztop0,mev)
def lfunc_mass_FSRQ(l,z):
    '''
    rest frame 1.4 GHz luminosity function for
    FSRQ AGN from Massardi+2010 (see table 1.)
    Args:
        l, float, luminosity (1.4GHz)
        z, float, redshift
    Returns:
        float, dN/dlogL1.4GHz/dV(h^3/Mpc^3)
    '''
    return lfunc_mass_agn(l,
                          z,
                          a=0.760,
                          b=2.508,
                          n0=10.**-10.382,
                          l0=10.**34.323,
                          kevo=-0.996,
                          ztop0=1.882,
                          dztop0=0.018,
                          mev=-0.166)
def lfunc_mass_bllac(l,z):
    '''
    rest frame 1.4 GHz luminosity function for BLLac
    AGN from Massardi+2010 (see table 1.)
    Args:
        l,float, luminosity (1.4GHz)
        z, float, redshift
    Returns:
        float, dN/dlogL1.4GHz/dV(h^3/Mpc^3)
    '''
    return lfunc_mass_agn(l,
                          z,
                          a=0.723,
                          b=1.618,
                          n0=10.**-6.879,
                          l0=10.**32.638,
                          kevo=0.208,
                          ztop0=1.282,
                          dztop0=0.,
                          mev=1.)
def lfunc_mass_ssagn(l,z):
    '''
    rest frame 1.4 GHz luminosity function for
    steep spectrum AGN in Massardi+2010. Table 1. params
    do not reproduce sensible counts and brightness temperature
    so I had to refit these parameters.
    Args:
        l,float,luminosity(1.4GHz)
        z, float, redshift
    Returns:
        float, dN/dlogL1.4GHz/dV(h^3/Mpc^3)
    '''
    return lfunc_mass_agn(l,
                          z,
                          a=6.09218673e-01,
                          b=1.94805350,
                          n0=8.30084893e-07,
                          l0=2.89354987e32,
                          kevo=3.20206060,
                          ztop0=1.00006968,
                          dztop0=6.82286615e-01,
                          mev=8.15187123e-02
                         )


#******************************************************
#Functions to compute the number of sources in a Comoving
#interval with observed fluxes Sobs.
#******************************************************
def dldsobs_agn(freq_obs,zem,alpha,F0=1.4e9):
    '''
    Compute the derivive of luminosity with respect to observed flux
    for agn.
    Args:
        freq_obs observed frequency
        zem, emission redshift.
        alpha, power law for emission: L=L0(f/F0)^alpha
        F0, reference frequency
    Returns: dL/dS_obs (W/Hz/Jy)
    '''
    output = 4.*PI*COSMO.luminosityDistance(zem)**2.\
    *(F0/(1.+zem)/freq_obs)**alpha/(1.+zem)/LITTLEH**2.#This is output in Mpc^2
    output=output*JY*(1e3*KPC)**2.
    return output


def l_agn(sobs,freq_obs,zem,alpha,F0=1.4e9):
    '''
    luminosity of AGN at redshift zem
    observed at frequency freq_obs
    Args:
        observed flux (Jy)
        freq_obs, observed frequency
        zem, agn redshift
        alpha, power law of emission L=L0(f/F0)^alpha
        F0, reference frequency (Hz)
    Returns:
        luminosity (W/Hz)
    '''
    return sobs*dldsobs_agn(freq_obs,zem,alpha,F0)
