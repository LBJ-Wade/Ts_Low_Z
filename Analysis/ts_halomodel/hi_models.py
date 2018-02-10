'''
file with models describing HI distribution within a dark matter halo.
all of these helper functions will be
'''
from settings import *
from hi_helpers import *

def rho_HI_padmanabhan2017_b_r(r,m,z,params):
    '''
    r, radius in comoving Mpc/h
    radial density function from Padmanabhan, Refregier, Amara
    in msolar/h/(Mpc/h)^3
    m: virial mass in Msol/h
    '''
    rho0=mHI_padmanabhan17_b(m,z,params)/vHI_padmanabhan17_b(m,z,params)
    rs=rVir(m,z)/concentration_padmanabhan17_b(m,z,params)
    return rho0*np.exp(-r/rs)



def rho_HI_padmanabhan2017_b_k(k,m,z,params):
    '''
    radial density function from Padmanabhan, Refregier, Amara
    2017
    m: virial mass in Msol/h
    z: redshift
    '''
    mhi=mHI_padmanabhan17_b(m,z,params)
    rs=rVir(m,z)/concentration_padmanabhan17_b(m,z,params)*(1.+z)/1e3
    return mhi/(1.+(k*rs)*2.)**2.
