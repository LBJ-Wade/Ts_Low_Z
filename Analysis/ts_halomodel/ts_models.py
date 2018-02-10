'''
functions for spin temperature
'''
from settings import *
from ts_helpers import *

def Ts0_MZR_Lilly2013(m,z,params):
    '''
    Predicts spin-temperature at center (star forming region) of
    a halo using the MZR from Lilly 2013.
    '''
    mstar=mstellar(m,z)
    return 10.**(params['ZSLOPE']*zmstar_lilly13(mstar,z,params)\
    +params['ZCONST'])

def Ts_MZR_Lilly2013_r(r,m,z,params):
    '''
    Compute the spin temperature as a function of radius in a dark
    matter halo.
    Args:
        r: radius (co-Mpc/h)
        m: virial mass
        z: redshift
        params, list of parameters
    Returns:
        Spin Temperature (K)

    '''
    return Ts0_MZR_Lilly2013(m,z,params)\
    *np.exp(r/params['RT']*1e3/(1.+z))
