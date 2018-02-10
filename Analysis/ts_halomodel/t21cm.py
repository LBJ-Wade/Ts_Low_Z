from settings import COSMO,F21,C,KBOLTZMANN,HPLANCK,M_INTERP_MAX,M_INTERP_MIN
from settings import MAXZ,MSOL,MP,KPC,PI,A10,LITTLEH,NINTERP_Z,TCMB
from _hi_helpers import HI_HELPERS
from _ts_models import TS_MODELS
import hi_helpers
import radio_background as rb
from utils import dict2tuple
from settings import SPLINE_DICT
import scipy.integrate as integrate
import scipy.interpolate as interp
from utils import massfunc
import numpy as np
#*****************************************
#functions for 21cm observables.
#*****************************************
def mTb_emiss(z21,rhohi=1.):
    '''
    mean emission brightness temperature (K)
    Args:
        z21, float,redshift
        rhohi, hi density in Msol Mpc^-3 h^2
    Returns:
        brightness temperature (K)
    '''
    return 3.*HPLANCK*1e-3*C**3.*A10/(32.*PI*KBOLTZMANN*MP*F21**2.)\
    *(1+z21)**2./COSMO.Hz(z21)/(KPC*1e3)**2.*MSOL*rhohi\
    *LITTLEH**2.

def mTb(z21,params):
    '''
    mean brightness temperature from 21cm emission
    Args:
        z21, float, redshift
        params, dictionary of parameters
    Returns:
        mean brightness temperature (K)
    '''
    splkey=('mTb','21cm')+dict2tuple(params)
    if splkey not in SPLINE_DICT.keys():
        zvals=np.linspace(0.,MAXZ,NINTERP_Z)
        mtbvals=np.zeros_like(zvals)
        if params['INCLUDE_ARCADE']:
            acoeff=1.
        else:
            acoeff=0.
        for znum,zval in enumerate(zvals):
            def g(x):
                rv=hi_helpers.rVir(10.**x,zval)
                rs=rv/HI_HELPERS[params['CSHIFUNC']](10.**x,zval,params)\
                *(1.+zval)/1e3
                rt=params['RT']*(1.+zval)/1e3
                rht=rs*rt/(rs+rt)
                volratio=hi_helpers.exp_vol(rv,rht)/hi_helpers.exp_vol(rv,rs)
                return massfunc(10.**x,zval)\
                *(HI_HELPERS[params['MHIFUNC']](10.**x,zval,params)-(1.+zval)\
                *volratio*(TCMB+rb.tb_agn_fast(zval)+rb.tb_sfg_fast(zval)\
                +acoeff*rb.tb_arc_const(F21/(1.+zval),params))\
                /TS_MODELS[params['TS0_FUNC']](10.**x,zval,params))
            mtbvals[znum]=mTb_emiss(zval,1.)*integrate.quad(g,M_INTERP_MIN,
            M_INTERP_MAX)[0]
        SPLINE_DICT[splkey]=interp.interp1d(zvals,np.log(mtbvals))
    return np.exp(SPLINE_DICT[splkey](z21))
