'''
helper functions for HI models
'''
import numpy as np
from settings import *
#****************************************************
#General cosmology functions relevant for HI distributions
##****************************************************

def omz(z):
    '''
    Omemga_m as a function of redshift.
    Args:
        z, redshift, float
    Returns:
        float, fraction of critical density of all ordinary matter.
    '''
    x=COSMO.Om0*(1+z)**3.
    x=x/(COSMO.Om0*(1+z)**3.+COSMO.Ode0+COSMO.Ok(0)*(1+z)**2.)
    return x

def delta(z):
    '''
    critical density from Bryan and Norman 1998
    Args:
         z, float, redshift
    Returns:
         collapsed virial density relative to mean density
    '''
    x=omz(z)-1.
    return 18.*PI*PI+82.*x-39.*x*x

def vVir(m,z):
    '''
    virial velocity of halo in km/sec
    from Maller and Bullock 2004
    Args:
         m,float,mass of halo in solar masses/h
         z,float,redshift
    Returns:
        float, circular virial velocity of dark
        matter halo (km/sec)
    '''
    return 144.*(omz(z)*delta(z)/97.2)**(1./6.)\
    *(m/1e12)**(1./3.)*(1.+z)**.5

def rVir(m,z):
    '''
    virial radius (coMpc/h)
    from Maller and Bullock 2004
    Args:
        m,float,mass from halo in msolar/h
        z,float,redshift
    Returns:
        virial radius of DM halo. (coMpc/h)
    '''
    return 206.*(omz(z)*delta(z)/97.2)**(-1./3.)\
    *(m/1e12)**(1./3.)*(1.+z)**-1.*1e-3*(1.+z)


#****************************************************
#Functions from Padmanabhan Refregier and Amara 2017.
#****************************************************

def mHI_padmanabhan17_b(m,z,params):
    '''
    HI mass - Halo mass relation from
    Padmanabhan Refregier and Amara 2017.
    Args:
        m, float, mass in msolar/h
        z, redshift
        params, dictionary of HI parameters
    Returns:
        hi mass for virial mass m at redshift z
    '''
    return params['ALPHA']*(1.-YP)\
    *COSMO.Ob0/COSMO.Om0*m*(m/1e11)**(params['BETA'])\
    *np.exp(-(10.**params['LOGV0']/vVir(m,z))**3.)

def concentration_padmanabhan17_b(m,z,params):
    '''
    Concentration parameter from
    Padmanabhan Refregier and Amara 2017.
    Args:
        m, float, virial mass in msolar/h
        z, redshift
        params, dictionary of HI model parameters.
    Returns:
        float, concentration parameter.
    '''
    return params['C0']*(m/LITTLEH/1e11)**-.109*4.\
    /(1.+z)**params['GAMMA']


def exp_vol(rv,rs):
    '''
    Volume of exponential density exp(-r/rs) extended to rv
    '''
    return 8.*PI*rs**3.-4.*np.exp(-rv/rs)*(2.*rs*rs+2*rv*rs+rv**2.)

def vHI_padmanabhan17_b(m,z,params):
    '''
    Volume of HI gas in exponential profile in proper kpc^3/h^3
    Args:
    '''
    rv=rVir(m,z)
    rs=rv/concentration_padmanabhan17_b(m,z,params)
    return exp_vol(rv,rs)
