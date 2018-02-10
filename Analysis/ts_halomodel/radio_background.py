from settings import *
import scipy.integrate as integrate
import scipy.interpolate as interp
from radio_background_helpers import *

def emiss_sfg(freq,z):
    '''
    the comoving emissivity (W/Hz/Mpc^3) from SFGs at
    observed frequency, freq at redshift z
    Returns: Joules/sec/Herz/Mpc^3
    Args:
        freq, float, frequency in Hz.
        z, float, reshift of emission
    Returns:
        emissivity in W*h^3/Hz/Comoving Mpc^3
    '''
    g=lambda x: lfunc_sfg(np.exp(x),freq,z)*np.exp(x)/np.log(10.)
    return integrate.quad(g,-500.,500.)[0]

def emiss_sfg_fast(zabsorbed,zemitted):
    '''
    the comoving emissivity (W/Hz/Mpc^3) from SFG at
    freq at redshift z
    Returns: Joules/sec/Herz/Mpc^3
    In this fast version, we pre-compute emissivity of 21-cm
    radiation observed at redshift z-obs emitted from redshift z.
    Args:
        zabsorbed, float, redshift of absorption
        zemitted, float, redshift of emission
    Returns:
        comoving emissivity (W/Hz/Mpc^3*h^3)
    '''
    emissvals=np.zeros((NINTERP_Z,NINTERP_Z))
    zaxis=np.linspace(0.,MAXZ,NINTERP_Z)
    if not 'SFG_EMISS' in SPLINE_DICT.keys():
        for znuma,zabs in enumerate(zaxis):
            for znume,zem in enumerate(zaxis):
                freq=F21/(1+zabs)*(1+zem)
                g=lambda x: lfunc_sfg(np.exp(x),freq,zem)*np.exp(x)/np.log(10.)
                emissvals[znuma,znume]=integrate.quad(g,-500.,500.)[0]
        SPLINE_DICT['SFG_EMISS']=\
        interp.RectBivariateSpline(zaxis,zaxis,np.log(emissvals))
    if zabsorbed>zemitted:
        return 0.
    else:
        return np.exp(SPLINE_DICT['SFG_EMISS'].ev(zabsorbed,zemitted))

def tb_sfg_fast(zabs,zmax=MAXZ):
    '''
    compute the brightness temperature of star forming galaxies.
    Args:
        zabs, float, absorption redshift
        zmax, float, Maximum redshift out to which background emission is integrated.
    Returns:
        float, brightness temperature (Kelvin)
    '''
    g=lambda x: emiss_sfg_fast(zabs,x)/(1+x)/COSMO.Ez(x)
    gint=integrate.quad(g,zabs,zmax)[0]
    freq_obs=F21/(1.+zabs)
    return  gint*(C/freq_obs)**2./(2*KBOLTZMANN*4.*PI)\
    *1e-3*C/COSMO.H0/(1e6*KPC**2.)*LITTLEH

def tb_sfg(freq_obs,zmin=0,zmax=MAXZ):
    '''
    compute the brightness temperature of SFGs
    Args:
        freq_obs, float, frequency of observation (Hz)
        zmin, minimum redshift to integrate sources from
        zmax, maximum redshift to integrate sources to
    Returns:
        brightness temperature (Kelvin)
    '''
    g=lambda x: emiss_sfg(freq_obs*(1.+x),x)/(1+x)/COSMO.Ez(x)
    gint=integrate.quad(g,zmin,zmax)[0]
    return  gint*(C/freq_obs)**2./(2*KBOLTZMANN*4.*PI)\
    *1e-3*C/COSMO.H0/(1e6*KPC**2.)*LITTLEH

def emiss_agn(freq,z):
    '''
    the comoving emissivity (W/Hz*h^3/Mpc^3) from AGN at
    freq at redshift z
    Returns: Joules/sec/Herz/Mpc^
    Args:
        freq, float, observed frequency (Hz)
        z, float, redshift
    Returns:
        comoving emissivity (W/Hz*h^3/Mpc^3)
    '''
    g=lambda x: (lfunc_mass_ssagn(np.exp(x),z)*(freq/1.4e9)**-0.8\
                 +lfunc_mass_bllac(np.exp(x),z)*(freq/1.4e9)**-0.1\
                +lfunc_mass_FSRQ(np.exp(x),z)*(freq/1.4e9)**-0.1)*np.exp(x)
    return integrate.quad(g,-500.,500.)[0]*ERG/np.log(10.)

def emiss_agn_fast(zabsorbed,zemitted):
    '''
    the comoving emissivity (W/Hz*h^3/Mpc^3) from AGN at
    freq at redshift z
    Returns: Joules/sec/Herz*h^3/Mpc^3
    In this fast version, we pre-compute emissivity of 21-cm
    radiation observed at redshift z-obs emitted from redshift z.
    Args:
        zabsorbed, redshift at which radiation is observed at redshifted 21cm
        zemitted, float, redshift at which radiation is emitted
    Returns:
        comoving emissivity from AGN (W/Hz*h^3/Mpc^3)
    '''
    emissvals=np.zeros((NINTERP_Z,NINTERP_Z))
    zaxis=np.linspace(0,MAXZ,NINTERP_Z)
    if not 'AGN_EMISS' in SPLINE_DICT.keys():
        for znuma,zabs in enumerate(zaxis):
            for znume,zem in enumerate(zaxis):
                freq=F21/(1+zabs)*(1+zem)
                g=lambda x: (lfunc_mass_ssagn(np.exp(x),zem)*(freq/1.4e9)**-0.8\
                             +lfunc_mass_bllac(np.exp(x),zem)*(freq/1.4e9)**-0.1\
                            +lfunc_mass_FSRQ(np.exp(x),zem)*(freq/1.4e9)**-0.1)*np.exp(x)

                emissvals[znuma,znume]=integrate.quad(g,-500.,500.)[0]*ERG/np.log(10.)
        SPLINE_DICT['AGN_EMISS']=interp.RectBivariateSpline(zaxis,zaxis,np.log(emissvals))
    if zabsorbed>zemitted:
        return 0.
    else:
        return np.exp(SPLINE_DICT['AGN_EMISS'].ev(zabsorbed,zemitted))

def tb_agn_fast(zabs,zmax=MAXZ):
    '''
    compute the brightness temperature of AGN
    Args:
        zabs, redshift of observation
    Returns:
        brightness temperature of AGN (kelvin)
    '''
    g=lambda x: emiss_agn_fast(zabs,x)/(1+x)/COSMO.Ez(x)
    gint=integrate.quad(g,zabs,zmax)[0]
    freq_obs=F21/(1.+zabs)
    return  gint*(C/freq_obs)**2./(2*KBOLTZMANN*4.*PI)\
    *1e-3*C/COSMO.H0/(1e6*KPC**2.)*LITTLEH

def tb_agn(freq_obs,zmin=0,zmax=MAXZ):
    '''
    compute the brightness temperature of AGN
    Args:
        freq_obs, float, frequency of observation (Hz)
        zmin, float, minimum redshift of sources to integrate from
        zmax, float, maximum redshif of sources to integrate too
    Returns:
        brightness temperature of AGN (kelvin)
    '''
    g=lambda x: emiss_agn(freq_obs*(1.+x),x)/(1+x)/COSMO.Ez(x)
    gint=integrate.quad(g,zmin,zmax)[0]
    return  gint*(C/freq_obs)**2./(2*KBOLTZMANN*4.*PI)\
    *1e-3*C/COSMO.H0/(1e6*KPC**2.)*LITTLEH

#******************************************************
#AGN and SFG count functions
#******************************************************
def dn_dlogs_dv(s,freq_obs,zem):
    '''
    comoving number density of AGN with observed fluxes
    s at freq_obs
    Args:
        s, observed flux (Jy)
        freq_obs, observed frequency (Hz)
        zem, agn redshift
    Returns: Comoving number density of AGN per logarithmic flux
    (h/Mpc)^3
    '''
    return lfunc_mass_bllac(l_agn(s,freq_obs,zem,-0.1)/ERG,zem)\
    +lfunc_mass_FSRQ(l_agn(s,freq_obs,zem,-0.8)/ERG,zem)\
    +lfunc_mass_ssagn(l_agn(s,freq_obs,zem,-0.1)/ERG,zem)
def dn_dlogs_domega(s,z,singles=False):
    '''
    Number of AGN per steradian with observed fluxes s
    at freq_obs=F21/(1.+z) with redshift greater than z
    Args:
        s, flux (Jy)
        z, redshift for counting sources above.
    Returns:
        number of sources per steradian with flux between S and S+dS
    '''
    if not singles:
        splkey=('dn_dlogs_domega')
        if not SPLINE_DICT.has_key(splkey):
            svals=np.logspace(S_INTERP_MIN,S_INTERP_MAX,NINTERP_S)
            zaxis=np.linspace(0.,MAXZ,NINTERP_Z)
            ngrid=np.zeros((NINTERP_S,NINTERP_Z))
            for znum,zval in enumerate(zaxis):
                for snum,sval in enumerate(svals):
                    g=lambda x: dn_dlogs_dv(sval,F21/(1.+zval),x)\
                    *dvc_dzdomega(x)
                    ngrid[snum,znum]=integrate.quad(g,zval,MAXZ)[0]
            SPLINE_DICT[splkey]=\
            interp.RectBivariateSpline(np.log(svals),zaxis,ngrid)
        return SPLINE_DICT[splkey].ev(np.log(s),z)
    else:
        g=lambda x: dn_dlogs_dv(s,F21/(1.+z),x)*dvc_dzdomega(x)
        return integrate.quad(g,z,MAXZ)

def var_flux(z):
    '''
    Calculate the variance in flux from all
    radio point sources with redshifts greater than z
    Args:
        z, redshift
    Returns:
        \int dS dN(>z)/dS/dOmega s^2 (Jy^2/Sr)
    '''
    return integrate.quad(lambda x: dn_dlogs_domega(np.exp(x))*np.exp(2.*x),
    S_INTERP_MIN*np.log(10.),S_INTERP_MAX*np.log(10.))[0]


#******************************************************
#Functions describing the brightness temperature of the
#extragalactic ARCADE excess
#******************************************************

def tb_arcade(freq_obs,params):
    '''
    observed arcade spectrum from Fixsen 2011
    Args:
        freq_obs, float, frequency of observation (Hz)
        alpha, float, power law of emission (beyond natural freq^-2
                      temperature evolution)
    Returns:
        brightness temperature of an extragalactic ARCADE excess (Kelvin)
    '''
    return (1.26-.23)*(freq_obs/1e9)**-(2.+params['ARCADE_POW'])
def emiss_arcade_delta(zc,alpha=.6):
    '''
    comoving emissivity of arcade excess at 1GHz
    assuming it is coming from a
    delta function at redshift zc (W/Hz*h^3/Mpc^3)
    Args:
        zc, float, center redshift.
        alpha, power law of emission.
    '''
    y=tb_arcade(1e9)*2.*KBOLTZMANN*4.*PI*COSMO.Ez(zc)*(KPC*1e3)**2.
    y/=(C*C*1e-18*1e-3*C*LITTLEH/COSMO.H0*(1+zc)**-(alpha+1))
    return y
def tb_arc_delta(freq_obs,zc=6.,alpha=.6):
    '''
    brightness temperature of arcade excess (freq_obs)
    assuming it is coming from a
    delta function at redshif zc
    Args:
        freq_obs, float, observed frequency (Hz)
        zc, float, center redshift of emission
        alpha, power law of emission
    Returns:
        brightness temperature (K) at observed frequency freq_observed.


    '''
    print emiss_arcade(zc,alpha)
    return C*C*1e-18*C*1e-3*(1+zc)**-(alpha+1)*\
    emiss_arcade(zc,alpha)*(freq_obs/1e9)**-(2.+alpha)\
    /(2.*(1e3*KPC)**2.*COSMO.Ez(zc)*4.*PI*KBOLTZMANN*COSMO.H0)*LITTLEH

def emiss_arc_const(z,params):
    '''
    comoving emissivity of arcade excess assuming it is uniformly distributed
    between redshifts zmin and zmax
    Args:
        z, redshift of emissivity
        params, dictionary
    Returns:
        comoving emissivity at 1GHz (W/Hz/(Mpc/h)^3)
    '''
    if z>=params['ZMIN_ARCADE'] and z<=params['ZMAX_ARCADE']:
        zint=integrate.quad(lambda x:1./((1+x)**(params['ARCADE_POW']+1.)\
        *COSMO.Ez(x)),params['ZMIN_ARCADE'],params['ZMAX_ARCADE'])[0]
        return (1e3*KPC)**2.*tb_arcade(1e9,params)*2.*KBOLTZMANN*4.*PI/(C*C*1e-18\
        *1e-3*C/COSMO.H0*LITTLEH)/zint
    else:
        return 0.

def tb_arc_const(freq_obs,params):
    '''
    brightness temperature at observed frequeency freq_obs of ARCADE 2 excess
    assuming that it is uniformly distributed between the redshifts of zmin and
    zmax.
    Args:
        freq_obs, float, observed frequeency
        params, dictionary
    Returns:
        brightness temperature at freq_obs
    '''
    zint=integrate.quad(lambda x:1./((1+x)**(params['ARCADE_POW']+1.)*COSMO.Ez(x))
    ,params['ZMIN_ARCADE'],params['ZMAX_ARCADE'])[0]
    z=F21/freq_obs-1.
    return C*C*1e-18*C*1e-3*emiss_arc_const(z,params)\
    *(freq_obs/1e9)**-(2.+params['ARCADE_POW'])\
    /(2.*(1e3*KPC)**2.*4.*PI*KBOLTZMANN*COSMO.H0)*zint*LITTLEH
