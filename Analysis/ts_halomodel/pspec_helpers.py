'''
helper functions for power spectrum
'''

import numpy as np
from settings import COSMO
from _hi_helpers import HI_HELPERS
from _ts_helpers import TS_HELPERS
from settings import M_INTERP_MAX,M_INTERP_MIN,SPLINE_DICT,MAXZ,NINTERP_Z
from utils import massfunc
from utils import bias
import hi_helpers
import scipy.integrate as integrate
import scipy.interpolate as interp
from _ts_hi_models import TS_HI_MODELS
from t21cm import mTb_emiss
from utils import dict2tuple

#************************************
#dark-matter power-spectrum functions
#************************************

def growth_function(z21):
    '''
    growth factor Omega_m(z)^0.55
    Args:
        z21, redshift
    Returns:
        dimensionless growth function
    '''
    hi_helpers.omz(z21)**.55

def virial_suppression(k,mu,z21,params):
    '''
    fingers of god effect causes suppression of
    small scales (see Anderson 2017) due to
    velocity smearing.
    Args:
        k, comoving wave-number (h/Mpc)
        z21, redshift
        params, parameter dictionary
    Returns:
        power-spectrum suppression factor arising from
        virial mostions of galaxies in a spectral
        line survey.
    '''
    mfactor=k*mu*params['SIGMA_V']/COSMO.H0
    return 1./(1.+(mfactor*k)**2.)

def rz_distortion(bvalue,mu,z21,params):
    '''
    redshift space distortions
    (see Andersion 2017)
    Args:
        bias, float, bias depending on structure tracer.
        z21, redshift
        params, parameter dictionary
    Returns:
        dimensionless redshift space distortion
        factor.
    '''
    return (1.+growth_function(z21)*mu**2./bvalue)

def power_lin(k,z):
    '''
    linear matter power spectrum
    Args:
        k, float, wavenumber (comoving h/Mpc)
        z, float, redshift
    Returns:
        linear matter power spectrum (Mpc/h)^3
    '''
    return COSMO.matterPowerSpectrum(k,z)

#*****************************************
#Density function for HI
#*****************************************
def rho_hi(z21,params):
    '''
    comoving density of HI in Msolar*h^3/Mpc^3
    Args:
        z21, float, redshift
        params, dictionary of parameters
    Returns:
        comoving density of HI in Msolar*h^3/Mpc^3
    '''
    splkey=('rho','hi')+dict2tuple(params)
    if not SPLINE_DICT.has_key(splkey):
        zvals=np.linspace(0,MAXZ,NINTERP_Z)
        rhovals=np.zeros_like(zvals)
        for znum, zval in enumerate(zvals):
            g=lambda x: HI_HELPERS[params['MHIFUNC']](10.**x,zval,params)\
            *massfunc(10.**x,zval)
            rhovals[znum]=integrate.quad(g,M_INTERP_MIN,M_INTERP_MAX)[0]
        SPLINE_DICT[splkey]=interp.interp1d(zvals,np.log(rhovals))
    return np.exp(SPLINE_DICT[splkey](z21))


#*****************************************
#bias functions
#*****************************************

def bias_hi(z21,params):
    '''
    bias of HI
    Args:
        z21, float, redshift
        params, dictionary of parameters.
    Returns:
        bias of HI (unitless)
    '''
    splkey=('bias','hi')+dict2tuple(params)
    if not SPLINE_DICT.has_key(splkey):
        zvals=np.linspace(0,MAXZ,NINTERP_Z)
        biasvals=np.zeros_like(zvals)
        for znum,zval in enumerate(zvals):
            g=lambda x:HI_HELPERS[params['MHIFUNC']](10.**x,zval,params)\
            *massfunc(10.**x,zval)*bias(10.**x,zval)
            biasvals[znum]=integrate.quad(g,M_INTERP_MIN,M_INTERP_MAX)[0]\
            /rho_hi(zval,params)
        SPLINE_DICT[splkey]=interp.interp1d(zvals,np.log(biasvals))
    return np.exp(SPLINE_DICT[splkey](z21))

def bias_hits(z21,params):
    '''
    bias of MHI/(Vhalo*Ts)
    Args:
        z21, float, redshift
        params, dictionary of parameters.
    Returns:
        bias of \int_vhalo rho_HI(r)/T_s(r)
    '''
    splkey=('bias_hits')+dict2tuple(params)
    if not SPLINE_DICT.has_key(splkey):
        zvals=np.linspace(0,MAXZ,INTERP_Z)
        biasvals=np.zeros_like(zvals)
        for znum,zval in enumerate(zvals):
            def g(x,bp=1.):
                rv=hi_helpers.rVir(10.**x,zval)
                rs=rv/HI_HELPERS['CSHIFUNC'](10.**x,zval,params)\
                *(1.+zval)/1e3
                rt=params['RT']*(1.+zval)/1e3
                rht=rs*rt/(rs+rt)
                volratio=hi_helpers.expvol(rv,rht)/hi_helpers.expvol(rv,rs)
                return HI_HELPERS['MHIFUNC'](10.**x,zvals,params)\
                /TS_HELPERS['TSFUNC'](10.**x,zvals,params)*volratio\
                *massfunc(10.**x,zval)*bias(10.**x,zval)**bp
            numer=integrate.quad(lambda x: g(x,1.),M_INTERP_MIN,M_INTERP_MAX)[0]
            denom=integrate.quad(lambda x: g(x,0.),M_INTERP_MIN,M_INTERP_MAX)[0]
            biasvals[znum]=numer/denom
        SPLINE_DICT[splkey]=interp.interp1d(zvals,np.log(biasvals))
    return np.exp(SPLINE_DICT[splkey](z21))



def ps_integral(k,z21,power_hi,power_hits,params):
    '''
    Evaluate integral
    \frac{3 h_p^2 c^3 A_{10} (1+z)^2)}{32 \pi \nu_{21}^2 m_p H(z)} \times
    (1+z)^power_hits \times \int_{m_min}^{m_max} dm n(m)\langle
     (\tilde{\rho}_{HI}^power_{HI}(k|m,z)^power_hi) \star
     (\tilde{T}_s(k|m,z))^-power_hits
     where \star denotes convolution.
     this integral appears in most 21-cm halo model terms.
    Args:
        k, wavenumber (h/Mpc)
        z21, redshift
        power_hi, power to raise FT(rho_HI) in convolution
        power_hits, power to raise FT(Ts) in convolution
        params, dictionary of parameters, including name of convolution function
    Returns:
        Integral in K^(2-power_hits) (Mpc/h)^(6-3*power_hi)
    '''
    g=lambda x: massfunc(10.**x,z21)\
    *TS_HI_MODELS[params['RHOTSMODEL']+'_k'](k,10.**x,z21,power_hi,power_hits,params)
    return integrate.quad(g,M_INTERP_MIN,M_INTERP_MAX)[0]\
    *mTb_emiss(z21)**(power_hi+power_hits)
    #*(1.+z21)**power_hits\#!!!! NEED TO MOVE (1+z) FACTOR
    #somewhere else!!!!
