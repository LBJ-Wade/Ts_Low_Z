'''
file with models describing Ts/HI distribution within a dark matter halo.
all of these helper functions will be
'''
from settings import PI
from _hi_helpers import HI_HELPERS
import _ts_models
import numpy as np

def rhoHI_over_Ts_exp_k(k,m,z,power_hi,power_hits,params):
    '''
    returns rhoHI^hi (rhoHI/Ts)^hits for the special case of
    both Ts and rhoHI being given by exponentials.
    Args:
        k, float, comoving wavenumber (h/Mpc)
        m, float, halo-mass (Msolar/h)
        z, float, redshift,
        power_hi, float, power to rais rhoHI to
        power_hits, float, power to raise rhoHI/Ts to
        params, parameters
    '''
    powertot=power_hits+power_hi
    powercorr=power_hi-params['ZSLOPE']*params['CORR_Z_HI']*power_hits
    rv=HI_HELPERS['rVir'](m,z)
    rs=rv/HI_HELPERS[params['CSHIFUNC']](m,z,params)
    rt=params['RT']*1e-3*(1.+z)
    rht=rs*rt/(rs+rt)
    ts=_ts_models.TS_MODELS[params['TS0_FUNC']](m,z,params)
    mhi=HI_HELPERS[params['MHIFUNC']](m,z,params)
    #moment for log normal is \exp(n*mu+n^2*sigma^2)
    hivol=HI_HELPERS['exp_vol'](rv,rs)
    tsvol=HI_HELPERS['exp_vol'](rv,rht)
    output=mhi**powertot*(tsvol/hivol)**power_hits\
    /((1.+(k*rs)**2.)**power_hi*(1.+(k*rht)**2.)**power_hits)**2.\
    *np.exp(.5*((-power_hits+power_hits**2.)*params['SIGMA_LOG_TS0'])**2.)\
    *np.exp(.5*(-powercorr+powercorr**2.)*params['SIGMA_LOG_HI']**2.)
    return output


def rhoHI_over_Ts_exp_r(r,m,z,params,include_temp=True):
    '''
    returns (rhoHI/Ts) for the special case of
    both Ts and rhoHI being given by exponentials.
    Args:
        r, float, comoving radius (h/Mpc)
        m, float, halo-mass (Msolar/h)
        z, float, redshift,
        include_temp, if True, divide by Ts. Otherwise,
        onyl return rhoHI
        params, parameters
    Returns:
        chomoving rhoHI/Ts in Msolar h^2/Mpc^3/K
    '''
    if include_temp: powert=1.
    else: powert=0.
    rv=HI_HELPERS['rVir'](m,z)
    rs=rv/HI_HELPERS[params['CSHIFUNC']](m,z,params)
    rt=params['RT']*1e-3*(1.+z)
    rht=1./(1./rs+powert/rt)
    ts=_ts_models.TS_MODELS[params['TS0_FUNC']](m,z,params)**powert
    mhi=HI_HELPERS[params['MHIFUNC']](m,z,params)
    hivol=HI_HELPERS['exp_vol'](rv,rs)
    output=mhi/ts*np.exp(-r/rht)/hivol
    return output
