'''
Functions for calculating power_spectra.
'''

from utils import fftlogbessel
from utils import massfunc
from utils import bias
from settings import SPLINE_DICT,NINTERP_KPERP,\
KPERP_INTERP_MIN,KPERP_INTERP_MAX,NINTERP_KPERP,COSMO,LITTLEH
from utils import dict2tuple
from pspec_helpers import *
from _radio_background import RADIO_BACKGROUND_MODELS
import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interp


#************************************
#The Radio Background power spectrum
#************************************
def ps_radio(kperp,z21,params):
    '''
    The power spectrum of known radio point sources
    Args:
        kperp, wavenumber perpindicular to line of sight (h/Mpc)
        z21, redshift where wavelength is 21cm.
        params, dictionary of parameters.
    Returns:
        power spectrum (K^2 h^3/Mpc^3)
    '''
    splkey=('ps_radio')+dict2tuple(params)

    if not mkey in SPLINE_DICT.keys():
        if params['INCLUDE_ARCADE']:
            arccoeff=0.
        else:
            arccoeff=1.
        kperpvals=np.logspace(K_PERP_INTERP_MIN,
        KPERP_INTERP_MAX,NINTERP_KPERP)
        psqvals=np.zeros_like(kperpvals)
        prefactor=(C/F21*(1+z21)/(1e3*KPC))**4/(64.*PI*PI*KBOLTZMANN**2.)\
        *LITTLEH*1e-3*C/cosmo.H0
        freq_obs=F21/(1.+z21)
        for k,kp in enumerate(kperpvals):
            g=lambda x: 1./((1+x)**2.*cosmo.Ez(x))\
            *(arccoeff*params['ARCADE_BIAS']\
            *RADIO_BACKGROUND_MODELS['emiss_arc_const'](x,params)\
            *(freq_obs*(1+x)/1e9)**(-params['ARCADE_POW'])\
            +RADIO_BACKGROUND_MODELS['emiss_agn_fast'](z21,x)\
            *bias(params['MEFF_AGN'],x)\
            +RADIO_BACKGROUND_MODELS['emiss_sfg_fast'](z21,x)\
            *bias(params['MEFF_SFG'],x))**2.\
            *power_lin(kp,x)
            psqvals[k]=integrate.quad(g,z21,MAXZ)[0]*prefactor
        SPLINE_DICT[splkey]\
        =interp.interp1d(np.log(kperpvals),np.log(psqvals))
    zf=1.
    if type(kperp)==np.ndarray:
        kmax=kperp>KPERP_INTERP_MAX
        kmin=kperp<KPERP_INTERP_MAX
        kint=np.logical_and(kperp>=K_PERP_INTERP_MIN,
        kperp<=KPERP_INTERP_MAX)
        output=np.zeros_like(kperp)
        output[kmax]=0.
        output[kmin]=\
        np.exp(SPLINE_DICT[splkey](np.log(KPERP_INTERP_MIN)))
        output[kint]=\
        np.exp(SPLINE_DICT[splkey](np.log(kperp[kint])))
        return output
    else:
        if kperp>K_PERP_INTERP_MIN:
            zf=0.
            kperp=K_PERP_INTERP_MAX
        if kperp<K_PERP_INTERP_MIN:
            kperp=K_PERP_INTERP_MIN
        return zf*np.exp(SPLINE_DICT[splkey](np.log(kperp)))

def corr_radio(rperp,z21,params):
    '''
    3d background quasar correlation function at redshift z.
    which is assumed to be only a function of rperp (no r_parallel)
    Args:
        rperp, float, perpindicular to line-of-sight distance (Mpc/h)
        z21, float,redshift
        params, parameters
    '''

    splkey=('corr_radio',z21)+dict2tuple(params)
    if splkey not in SPLINE_DICT.keys():
        g=lambda x: ps_radio(x,z,params)
        #two factors of x to integrate by dlogx
        #don't multiply by 2*PI since this is divided out by fourier convention
        raxis,cvals=fftlogbessel(g,nrperp,K_PERP_INTERP_MIN-1,K_PERP_INTERP_MAX+1,tdir=-1)
        #print raxis.min()
        #print(minlog)
        #print('raxis.min='+str(raxis.min()))
        SPLINE_DICT[splkey]=interp.interp1d(np.log(raxis),cvals)
    if isinstance(r,np.ndarray):
        rmin=r<10**minlog
        rmax=r>10**maxlog
        rint=np.logical_and(r>=10.**K_PERP_INTERP_MIN,r<=10.**K_PERP_INTERP_MAX)
        output=np.zeros_like(r)
        output[rmin]=SPLINE_DICT[splkey](np.log(10**K_PERP_INTERP_MIN))
        output[rmax]=0.
        output[rint]=SPLINE_DICT[splkey](np.log(r[rint]))
        return output
    else:
        if r<10.**K_PERP_INTERP_MIN:
            zf=1.
            r=10.**K_PERP_INTERP_MIN
        elif r>10.**K_PERP_INTERP_MAX:
            zf=0.
            r=10.**K_PERP_INTERP_MAX
        else:
            zf=1.
        return zf*SPLINE_DICT[splkey](np.log(r))

#************************************************
#21-cm power spectrum
#************************************************

def ps_21_1h(k,z21,params,singlek=False):
    '''
    1-halo 21-cm power-spectrum
    Args:
        k: wavenumber (h/Mpc)
        z21: redshift
        params, parameters,
        singlek, bool. If True, avoid generating interpolated ps
        for entire redshift.
    Returns:
        1-halo term of 21cm emission power spectrum
        K^2 Mpc^3/h^3
    '''
    if not singlek:
        splkey=('ps_21_1h',z21)+dict2tuple(params)
        if not splkey in  SPLINE_DICT:
            kaxis=np.logspace(K_PERP_INTERP_MIN,K_PERP_INTERP_MAX,NPTS)
            psvals=np.zeros_like(kaxis)
            for knum,kval in enumerate(kaxis):
                psvals[knum]=ps_integral(kval,z21,params,2.,0.)
            SPLINE_DICT[splkey]=interp.interp1d(np.log(kaxis),np.log(psvals))
        return np.exp(SPLINE_DICT[splkey](np.log(k)))
    else:
        return ps_integral(k,z21,2.,0.,params)

def ps_21_2h(k,z21,params,singlek=False):
    '''
    2-halo 21-cm power-spectrum
    Args:
        k,wavenumber (h/Mpc)
        z21, redshift
        params, parameters,
        singlek, bool. If True, avoid generating interpolated ps
        for entire redshift.
    Returns:
        2-halo term for 21cm emission power spectrum
        K^2 Mpc^3/h^3
    '''
    if not singlek:
        splkey('ps_21_2h',z21)+dict2tuple(params)
        if not splkey in SPLINE_DICT:
            kaxis=np.logspace(K_PERP_INTERP_MIN,K_PERP_INTERP_MAX,NPTS)
            psvals=np.zeros_like(kaxis)
            for knum,kval in enumerate(kaxis):
                psvals[knum]=np.abs(ps_integral(kval,z21,params,1.,0.))**2.\
                *power_lin(kval,z21)
            SPLINE_DICT[splkey]=interp.interp1d(np.log(kaxis),np.log(psvals))
        return np.exp(SPLINE_DICT[splkey](np.log(k)))
    else:
        return np.abs(ps_integral(k,z21,1.,0.,params))**2.*power_lin(k,z21)
def ps_21_iso(k,z21,params,singlek=False):
    '''
    returns the isotropic 21-cm emission power-spectrum
    Args:
        k, wavenumber (h/Mpc)
        z21, redshift
        params, parameters
        singlek, bool, If True, avoid generating interpolation splines.
    Returns:
        isotropic 21cm power spectrum (Mpc/h)^3 K^2
    '''
    return ps_21_1h(k,z21,params,singlek)+ps_21_2h(k,z21,params,singlek)

def ps_21(kperp,kpara,z21,params,singlek=False):
    '''
    anisotropic 21-cm powers-spectrum
    Args:
        kperp, wavenumber perp to LoS (h/Mpc)
        kpara, wavenumber para to LoS (h/Mpc)
        z21, redshift
        params, parameter dictionary
        singlek, on-the-fly calculation (no precomputing splines)
    Returns:
        21-cm power spectrum (Mpc/h)^3 K^2
        with anisotropic effects
    '''

    k=np.sqrt(kperp**2.+kpara**2.)
    mu=np.arcos(kpara/k)
    return ps_21_iso(k,z21,params,singlek)*rz_distortion(bias_hi(z21,params),
    mu,z21,params)**2.*virial_suppression(k,mu,z21,params)

def ps_tau_1h(k,z21,params,singlek=False):
    '''
    1-halo power-spectrum of tau
    Args:
        k,wavenumber (h/Mpc)
        z21,redshift
        params, parameters
        singlek, bool. If True, avoid generating interpolated ps
        for entire redshift.
    Returns:
        1-halo term fo tau power-spectrum
        Mpc^3/h^3
    '''
    if not singlek:
        splkey=('ps_tau_1h',z21)+dict2tuple(params)
        if not SPLINE_DICT.has_key(splkey):
            kaxis=np.logspace(K_PERP_INTERP_MIN,K_PERP_INTERP_MAX,NPTS)
            psvals=np.zeros_like(kaxis)
            for knum,kval in enumerate(kaxis):
                psvals[knum]=ps_integral(kval,z21,2.,2.,params)
            SPLINE_DICT[splkey]=interp.interp1d(np.log(kaxis),np.log(psvals))
        return np.exp(SPLINE_DICT[splkey](np.log(k)))
    else:
        return ps_integral(k,z21,2.,2.,params)

def ps_tau_2h(k,z21,params,singlek=False):
    '''
    2-halo power-spectrum for tau
    Args:
        k,wavenumber (h/Mpc)
        z21, redshift
        params, parameters
        singlek, bool. If True, avoid generating interpolated ps
        for entire redshift
    Returns:
        2-halo term for tau power-spectrum (Mpc^3/h^3)
    '''
    if not singlek:
        splkey=('ps_tau_2h',z21)+dict2tuple(params)
        if not SPLINE_DICT.has_key(splkey):
            kaxis=np.logspace(K_PERP_INTERP_MIN,K_PERP_INTERP_MAX,NPTS)
            psvals=np.zeros_like(kaxis)
            for knum,kval in enumerate(kaxis):
                psvals[knum]=np.abs(ps_integral(kval,z21,1.,1.,params))**2.\
                *power_lin(kval,z21)
            SPLINE_DICT[splkey]=interp.interp1d(np.log(kaxis),np.log(psvals))
        return np.exp(SPLINE_DICT[splkey](np.log(k)))
    else:
        return p.abs(ps_integral(k,z21,1.,1.,params))**2.\
        *power_lin(k,z21)

def ps_tau_iso(k,z21,params,singlek=False):
    '''
    returns the isotropic optical-depth power-spectrum
    Args:
        k, wavenumber (h/Mpc)
        z21, redshift
        params, parameters
        singlek, bool, If True, avoid generating interpolation splines.
    Returns:
        isotropic optical-depth power spectrum (Mpc/h)^3
    '''
    return ps_tau_1h(k,z21,params,singlek)+ps_tau_2h(k,z21,params,singlek)

def ps_tau(kperp,kpara,z21,params,singlek=False):
    '''
    anisotropic 21-cm powers-spectrum
    Args:
        kperp, wavenumber perp to LoS (h/Mpc)
        kpara, wavenumber para to LoS (h/Mpc)
        z21, redshift
        params, parameter dictionary
        singlek, on-the-fly calculation (no precomputing splines)
    Returns:
        21-cm power spectrum (Mpc/h)^3 K^2
        with anisotropic effects
    '''

    k=np.sqrt(kperp**2.+kpara**2.)
    mu=np.arcos(kpara/k)
    return ps_tau_iso(k,z21,params,singlek)*rz_distortion(bias_hits(z21,params),
    mu,z21,params)**2.*virial_suppression(k,mu,z21,params)

def ps_21tau_1h(k,z21,params,singlek=False):
    '''
    1-halo power-spectrum for 21-tau cross term
    Args:
        k,wavenumber (h/Mpc)
        z21, redshift
        params,parameters
        singlek, bool. If True, avoid generating interpolated ps
        for entire redshift
    Returns:
        1-halo term for 21-tau cross-power spectrum (K Mpc^3/h^3)
    '''
    if not singlek:
        splkey=('ps_21tau_1h',z21)+dict2tuple(params)
        if not SPLINE_DICT.has_key(splkey):
            kaxis=np.logspace(K_PERP_INTERP_MIN,K_PERP_INTERP_MAX,NPTS)
            psvals=np.zeros_like(kaxis)
            for knum,kval in enumerate(kaxis):
                psvals[knum]=np.abs(ps_integral(kval,z21,2.,1.,params))
            SPLINE_DICT[splkey]=interp.interp1d(np.log(kaxis),np.log(psvals))
        return np.exp(SPLINE_DICT[splkey](np.log(k)))
    else:
        return np.abs(ps_integral(k,z21,2.,1.,params))

def ps_21tau_2h(k,z21,params,singlek=False):
    '''
    2-halo power-spectrum for 21-tau cross term
    Args:
        k,wavenumber (h/Mpc)
        z21,redshift
        params, parameters
        singlek, bool, if True, only perform calculation for single
                       k at fixed redshift.
    Returns:
        2-halo term for 21-tau cross-power spectrum (K Mpc^3/h^3)
    '''
    if not singlek:
        splkey=('ps_21tau_2h',z21)+dict2tuple(params)
        if not SPLINE_DICT.has_key(splkey):
            kaxis=np.logspace(K_PERP_INTERP_MIN,K_PERP_INTERP_MAX,NPTS)
            psvals=np.zeros_like(kaxis)
            for knum,kval in enumerate(kaxis):
                psvals[knum]=ps_intergral(kval,z21,1.,1.,params)\
                *ps_integral(kval,z21,1.,0.,params)*power_lin(kval,z21)
            SPLINE_DICT[splkey]=interp.interp1d(np.log(kaxis),np.log(psvals))
        return np.exp(SPLINE_DICT[splkey](np.log(k)))
    else:
        return ps_integral(k,z21,1.,1.,params)*ps_integral(k,z21,1.,0.,params)\
        *power_lin(k,z21)

def ps_21tau_iso(k,z21,params,singlek=False):
    '''
    returns the isotropic 21cm-optical-depth
    cross power-spectrum
    Args:
        k, wavenumber (h/Mpc)
        z21, redshift
        params, parameters
        singlek, bool, If True, avoid generating interpolation splines.
    Returns:
        isotropic 21cm-optical-depth power spectrum K (Mpc/h)^3
    '''
    return ps_21tau_1h(k,z21,params,singlek)+ps_21tau_2h(k,z21,params,singlek)

def ps_21tau(kperp,kpara,z21,params,singlek=False):
    '''
    anisotropic 21-cm tau cross powers-spectrum
    Args:
        kperp, wavenumber perp to LoS (h/Mpc)
        kpara, wavenumber para to LoS (h/Mpc)
        z21, redshift
        params, parameter dictionary
        singlek, on-the-fly calculation (no precomputing splines)
    Returns:
        21-cm-optical-depth cross power spectrum K (Mpc/h)^3 K^2
        with anisotropic effects
    '''

    k=np.sqrt(kperp**2.+kpara**2.)
    mu=np.arcos(kpara/k)
    return rz_distortion(bias_hits(z21,params),mu,z21,params)\
    *rz_distortion(bias_hi(z21,mu,z21,params))\
    *virial_suppression(k,mu,z21,params)*ps_21tau_iso(k,z21,params,singlek)
