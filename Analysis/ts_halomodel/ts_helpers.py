'''
Helper functions for spin-temperature model
These include prescriptions for mass metallicity relation
'''
from settings import COSMO,KPC,YEAR,LITTLEH
import numpy as np
#********************************************* *
#Functions for stellar-mass halo-mass relation
#**********************************************
def fmstar(x,alpha,gamma,delta):
    '''
    equation 3b from Behroozi et al. 2013 with
    best fit parameters.
    '''
    return -np.log10(10**(alpha*x)+1.)+delta\
    *(np.log10(1.+np.exp(x)))**gamma/(1.+np.exp(10**-x))

def mstellar(m,z):
    '''
    equation 3a from Behroozi et al. 2013 with
    best fit parameters
    Args:
        m, float, stellar mass in msolar/h
        z, float, redshift.
    Returns:
        float, stellar mass associated with halo.
    '''
    mvir=m/LITTLEH
    a=1./(1.+z)
    nu=np.exp(-4.*a*a)
    logeps=-1.777+(-0.006*(a-1.)+(-.000*z))*nu-0.119*(a-1.)#
    logm1=11.514+(-1.793*(a-1)+(-0.251)*z)*nu
    alpha=-1.412+(0.731*(a-1.))*nu
    delta=3.508+(2.608*(a-1)+(-0.043)*z)*nu
    gamma=0.316+(1.319*(a-1)+0.279*z)*nu
    f0=fmstar(0.,alpha,gamma,delta)
    x=np.log10(mvir)-logm1
    fm=fmstar(x,alpha,gamma,delta)
    logmstar=logeps+logm1+fm-f0
    return 10.**logmstar*LITTLEH

#**************************************************
#Functions for metallicity-stellar mass relation
#**************************************************
def zmstar_lilly13(mstar,z,params):
    '''
    equation 40 from Lilly et al. 2013
    Args:
        mstar,float,stellar mass (msolar/h)
        z, float, redshift
        params, dictionary of parameters
    Returns:
        float, metallicity of galactic gas following the MZR from
        Lilly et al. 2013.
    '''
    y=10.**params['LILLY13']['LOGYIELD']
    z0=params['LILLY13']['Z0DY']*y
    lam=params['LILLY13']['LOADING10']\
    *(mstar/LITTLEH/1e10)**params['LILLY13']['LOADINGPOW']
    eps=params['LILLY13']['DEPLETION10']**-1.\
    *(mstar/LITTLEH/1e10)**params['LILLY13']['DEPLETIONPOW']\
    *(1+z)**params['LILLY13']['DEPLETION_EV_POW']
    r=params['LILLY13']['RETURN_FRAC']
    alphaeps=params['LILLY13']['DEPLETION_EV_POW']
    if z<=2:
        ssfr=1./(1.-r)*0.07*(mstar/LITTLEH/10**10.5)**(-0.1)*(1.+z)**3.
        dlogmudz=(3.-alphaeps)/(1+z)
    else:
        ssfr=1./(1.-r)*0.3*(mstar/LITTLEH/10**10.5)**(-0.1)*(1.+z)**(5./3.)
        dlogmudz=(5./3.-alphaeps)/(1+z)
    mu=ssfr/eps
    dzdt=-COSMO.Hz(z)*(1.+z)/KPC*1e9*YEAR#convert H(z) to Gyr^-1
    zeq=np.log10(z0+y/(1.+lam/(1.-r)+mu+dlogmudz*dzdt/eps/(1.-r)))
    return zeq
