'''
Code for calculating rates for absorbing systems.
'''
from utils import massfunc
import scipy.integrate as integrate
from hi_helpers import rVir, vVir
from settings import KPC, SPLINE_DICT, LITTLEH,COSMO,KBOLTZMANN
from settings import MSOL, MP, F21, NHIDT2TAU,PI,C,S_INTERP_MAX
from settings import M_INTERP_MAX,M_INTERP_MIN,TAU_INTERP_MIN,TAU_INTERP_MAX
from settings import SNR_INTERP_MIN,SNR_INTERP_MAX,N_INTERP_SNR,JY
from settings import N_INTERP_RTAU,N_INTERP_M,R_INTERP_MIN, N_INTERP_TAU
from _ts_hi_models import TS_HI_MODELS
import numpy as np
from utils import dict2tuple
from scipy.special import erfinv
import scipy.interpolate as interp
import scipy.optimize as op
import radio_background as rb

def nhidt(b,m,z,params,include_temp=True):
    '''
    Compute NHI/Ts as a function of impact parameter.
    Args:
        b, impact parameter (comoving Mpc/h)
        m, halo virial mass (msolar/h)
        z, redshift
        params, dictionary of parameters.
    Returns:
        proper number density of HI atoms over spin temperature
        in N/m^2/K
    '''
    rv=rVir(m,z)#rvirial in Mpc/h
    g=lambda x: TS_HI_MODELS[params['RHOTSMODEL']+'_r']\
    (np.sqrt(b**2.+x**2.),m,z,params,include_temp)
    if b<rv:
        output=2.*integrate.quad(g,0.,np.sqrt(rv**2.-b**2.))[0]
    else:
        output=0.
    #Now convert to NHI/(Mpc/h)^2, h below converts from Msol/h to Msol
    output=output*MSOL/MP/LITTLEH
    #Now convert to NHI/m^2 (proper)
    output=output/(KPC*1e3/LITTLEH)**2.*(1.+z)**2.
    return output

def sigma_line(m,z,params):
    '''
    Give 21-cm line width for halo mof mass m at redshift z.
    Using Haehneldt 2014 relationship
    Args:
        m, halo mass (Msolar/h)
        z, redshift
        params, dictionary of simulation parameters
    returns:
        standard deviation of line (km/sec)
    '''
    return vVir(m,z)*params['XPEAK']/(2.*np.sqrt(2.)*erfinv(0.9))#km/sec


def tau_gaussian(b,m,z,offset,params,dlogr=0,integrated=False,
                units='Hz',singler=False):
    '''
    optical depth profile for a halo of mass m
    Args:
        b (impact parameter in comoving Mpc/h)
        m, virial mass of halo (msol/h)
        offset, halo rest-frame frequency (Hz)
                or velocity (km/sec)
        params, dictionary of parameters
        singler, evaluate at single r value
        (rather than evaluating interp tables)
    Returns:
        optical depth
    '''
    sigma_v=sigma_line(m,z,params)
    if not singler:
        splkey=('tau_gaussian',z)+dict2tuple(params)
        if splkey not in SPLINE_DICT:
            maxis=np.logspace(M_INTERP_MIN,M_INTERP_MAX,N_INTERP_M)
            raxis=np.logspace(R_INTERP_MIN,0,N_INTERP_RTAU)
            taumaxvals=np.zeros((N_INTERP_M,N_INTERP_RTAU))
            for mnum,mval in enumerate(maxis):
                for rnum,rval in enumerate(raxis):
                    taumaxvals[mnum,rnum]\
                    =NHIDT2TAU\
                    *nhidt(rval*rVir(mval,z),
                    mval,z,params)
            SPLINE_DICT[splkey]\
            =interp.RectBivariateSpline(np.log(maxis),
                                        np.log(raxis),
                                        taumaxvals)
        r=b/rVir(m,z)
        sf=1.
        if isinstance(r,np.ndarray):
            r[r<10.**R_INTERP_MIN]=10.**R_INTERP_MIN
            output=np.zeros_like(r)
            output[r<=1.]=SPLINE_DICT[splkey].ev(np.log(m),np.log(r[r<=1.]),
            dy=dlogr)
            output[r>1.]=0.
        else:
            if r<10.**R_INTERP_MIN:
                r=10.**R_INTERP_MIN
            elif r>1.:
                r=1.
                sf=0.
            output=sf*SPLINE_DICT[splkey].ev(np.log(m),np.log(r),
            dy=dlogr)
    else:
        output=NHIDT2TAU*nhidt(b,m,z,params)
    if integrated:
        return output
    else:
        if units=='Hz':
            freq=offset
            velocity=np.abs(freq/F21-1.)*C*1e-3
        else:
            velocity=offset
        output=output*np.exp(-(velocity/sigma_v)**2./2.)
        output=output*C*1e-3/np.sqrt(2.*PI)/sigma_v/F21
    return output


<<<<<<< HEAD
def r_tau(tau,m,z,params,dtau=0.,recompute=False):
=======
def r_tau(tau,m,z,params,recompute=False,dtau=0.):
>>>>>>> 8b4ee19c000dacdb634c7ffde7a05e3ef4ff1dc1
    '''
    Compute the radius of a halo within which all LoSs subtend a maximum.
    tau'>tau.
    Args:
        tau: minimum tau within radius r.
        m: virial mass of halo.
        z, redshift
        params, dictionary of model parameters
        dlogtau, order of derivative with respect to logtau
    Returns:
        r, radius in comoving Mpc/h within which all otpical depths are
        greater than tau.
    '''
    assert dtau in [0.,1.,0,1]
    splkey=('r_tau',z)+dict2tuple(params)
    if not splkey in SPLINE_DICT or recompute:
        maxis=np.logspace(M_INTERP_MIN,M_INTERP_MAX,N_INTERP_M)
        tauaxis=np.logspace(TAU_INTERP_MIN,TAU_INTERP_MAX,N_INTERP_TAU)[::-1]
        rvals=np.ones((N_INTERP_M,N_INTERP_RTAU))*10.**R_INTERP_MIN
        for mnum,mval in enumerate(maxis):
            r0=R_INTERP_MIN
            rv=rVir(mval,z)
            if tau_gaussian(rv*10.**r0,mval,z,F21,params)>10.**TAU_INTERP_MIN:
                for taunum,tauval in enumerate(tauaxis):
                    g=lambda x:np.abs(np.log(tau_gaussian(rv*10.**x,
                    mval,z,F21,params)/tauval/rv/10.**x))
                    res=op.minimize(g,x0=[r0],bounds=[[r0,0.]],method='SLSQP').x
                    r0=res[0]
                    rvals[mnum,taunum]=10.**r0
        SPLINE_DICT[splkey]=\
        interp.RectBivariateSpline(np.log(maxis),
        np.log(tauaxis[::-1]),np.log(rvals[:,::-1]))
    if dtau==1.:
        #return dr/dtau
        return SPLINE_DICT[splkey].ev(np.log(m),
        np.log(tau),dy=dtau)/tau*r_tau(tau,m,z,params)
    else:
        return np.exp(SPLINE_DICT[splkey].ev(np.log(m),
        np.log(tau)))*rVir(m,z)

def sigma_tau(tau,m,z,params,recompute=False):
    '''
    Area subtended by a halo, mass m, redshift z,
    within which all lines of sight have an optical
    depth to 21cm that is greater than tau.
    Args:
        tau, optical depth
        m, halo mass (msolar/h)
        z, redshift
        params, dictionary of simulation parameters.
    Returns: Subtended area in comoving Mpc/h^2
    '''
    return PI*r_tau(tau,m,z,params,recompute=recompute)**2.

def dsigma_dtau(tau,m,z,params,recompute=False):
    '''
    derivative of area subtended by halo of mass m, redshift z
    with respect to optical depth at radius where optical depth equals tau
    Args:
        tau, float, optical depth at which derivative is evaluated
        m, halo-mass (msolar/h)
        z, redshift
        params, dictionary of parameters
    Returns:
        dsigma/dtau (Mpc/h)^2
    '''
    return 2.*PI*r_tau(tau,m,z,params,recompute=recompute)\
    *r_tau(tau,m,z,params,dtau=1.,recompute=recompute)


def dsigma_dtau(tau,m,z,params):
    '''
    derivative of area subtended by halo of mass m, redshift z
    with respect to optical depth at radius where optical depth equals tau
    Args:
        tau, float, optical depth at which derivative is evaluated
        m, halo-mass (msolar/h)
        z, redshift
        params, dictionary of parameters
    Returns:
        dsigma/dtau (Mpc/h)^2
    '''
    return 2.*PI*r_tau(tau,m,z,params)*r_tau(tau,m,z,params,dtau=1.)


def rline(m,z,channel_width):
    '''
    ratio between effective standard deviation in gaussian line-profile
    for HI absorption convolved with channelization function
    that is a Gaussian with width channel_width, divided by channel width.
    Args:
        m, mass of host halo.
        z, redshift
        channel-width,width of channel (Hz)
    Returns:
        ratio between convolved velocity width and channel width
    '''
    return np.sqrt(1.+(channel_width*C/F21/2e3/sigma_line(m,z,params))**2.)


def d_tau_variance_dm(m,z,params):
    '''
    \int dsigma/dtau (tau|m,z) tau^2 dtau
    Args:
        m, mass (msolar/h)
        z, redshift
        channel_width, channel size in Hz
    Returns:
        second moment of optical depth weighted by area at optical depth
        (Mpc/h)^2
    '''
    splkey=('tau_variance_int',z)+dict2tuple(params)
    if not SPLINE_DICT.has_key(splkey):
        maxis=np.logspace(M_INTERP_MIN,M_INTERP_MAX,N_INTERP_M)
        varvals=np.zeros_like(maxis)
        for mnum,mval in enumerate(maxis):
            g=lambda x: np.abs(dsigma_dtau(10.**x,m,z))*10.**(3.*x)\
            *sigma_line(10.**x,z,params)/channel_width*2.
            varvals[mnum]=integrate.quad(g,TAU_INTERP_MIN,TAU_INTERP_MAX)[0]
        SPLINE_DICT[splkey]=interp.interp1d(np.log(maxis),np.log(varvals))
    return np.exp(SPLINE_DICT[splkey](np.log(m)))


def tau_variance(z,channel_width,params):
    '''
    \int dm dtau dsigma/dtau(tau,m,z) tau^2
    Args:
        z, redshift
        channel_width, width of channels used to measure lines.
    Returns:
        d \sum \tau^2/d r_comoving (h/Mpc)
    '''
    g=lambda x: massfunc(10.**x,z)*d_tau_variance_dm(10.**x,z,params)
    return integrate.quad(g,M_INTERP_MIN,M_INTERP_MAX)

def v_obs_variance(z,dz,survey_area,rms,channel_width,params):
    '''
    Determine variance for
    '''


def tau_limit(snr,m,z,channel_width,params):
    '''
    compute limiting tau in dark matter halo of mass
    m at reshift z detectable against a source with
    observed at significance snr.
    Args:
        snr, snr of observed source
        m, mass of halo (Msol/h).
        channel_width (Hz)
        params, dictionary of params
    Returns:
        limiting detectable optical depth
    '''
    rchan=rline(m,z,channel_width)
    output=-np.log(np.max([1.-rchan/snr,1e-100]))
    return output

def dn_dr(snr,channel_width,z,params,singlesnr=False):
    '''
    number of absorbers per comoving LoS distance interval where fractional flux
    difference is larger than snr times the noise level.
    Args:
        snr, signal to noise ratio
        z, redshift
        params, dictionary of model parameters
    Returns:
        comoving number density of absorbers per line-of-sight interval (h/Mpc)
    '''
    if not singlesnr:
        splkey=('dn_dr',channel_width,z)+dict2tuple(params)
        if not SPLINE_DICT.has_key(splkey):
            snrvals=np.logspace(SNR_INTERP_MIN,SNR_INTERP_MAX,N_INTERP_SNR)
            dndrvals=np.zeros_like(snrvals)
            for snrnum,snrval in enumerate(snrvals):
                g = lambda x: massfunc(10.**x,z)\
                *sigma_tau(tau_limit(snrval,10.**x,z,channel_width,params),
                10.**x,z,params)
                dndrvals[snrnum]=integrate.quad(g,M_INTERP_MIN,M_INTERP_MAX)[0]
            SPLINE_DICT[splkey]=interp.interp1d(np.log(snrvals),dndrvals)
        sf=1.
        if type(snr)==np.ndarray:
            snr[snr<10.**SNR_INTERP_MIN]=SNR_INTERP_MIN
            snr[snr>=10.**SNR_INTERP_MAX]=10.**SNR_INTERP_MAX
            output=SPLINE_DICT[splkey](np.log(snr))
            output[snr<10.**SNR_INTERP_MIN]=0.
        else:
            if snr<10.**SNR_INTERP_MIN:
                snr=10.**SNR_INTERP_MIN
                sf=0.
            elif snr>=10.**SNR_INTERP_MAX:
                sf=1.
                snr=10.**SNR_INTERP_MAX
            output=SPLINE_DICT[splkey](np.log(snr))
        return sf*output
    else:
        g = lambda x: massfunc(10.**x,z)\
        *sigma_tau(tau_limit(snr,10.**x,z,channel_width,params),
        10.**x,z,params)
        return integrate.quad(g,M_INTERP_MIN,M_INTERP_MAX)[0]



def sefd(z,r_ant,trx,eta=0.5):
    '''
    sefd for antenna observing 21cm at redshift z
    Args:
        redshift of 21cm line.
        r_ant, radius of antenna (float)
        trx, receiver temperature
        eta, aperture efficiency
    Returns: SEFD in Jy
    '''
    freq_obs=F21/(1.+z)
    aeff=eta*r_ant**2.*PI
    tsys=60.*(C/freq_obs)**2.55+trx
    return KBOLTZMANN*tsys/2./JY/aeff/eta

def rms_chan(z,nant,r_ant,chan_w,t_int,trx,npol=2,eta=0.5):
    '''
    rms in Jy/beam
    Args:
        redshift of 21cm line.
        nant, number of antennas (float)
        r_ant, radius of antenna (float)
        channel_width, width of integration channel (Hz)
        integration_time, length of integration per field (seconds)
        trx, receiver temperature
        eta, aperture efficiency
    Returns: Amplitude of noise (Jy/beam)
    '''
    return sefd(z,r_ant,trx,eta=eta)/\
    np.sqrt(npol*nant*(nant-1.)*t_int*chan_w)

def instrument_counts(z,survey_area,smin,smax,channel_width,rms,params):
    '''
    Calculate the number of absorbers per comoving redshift interval
    detected at the nsigma level.
    Args:
        z, redshift
        survey_area, area of survey (steradians)
        channel_width, specral resolution (Hz)
        rms, noise level/beam per spectral channel (Jy/bm)
        params, dictionary of model parameters
        nsigma, detection threshold above rms.
    Returns:
        Number of detected absorbers per redshift interval.
    '''
    s_min=rms*10.**SNR_INTERP_MIN
    g = lambda x: rb.dn_dlogs_domega(10.**x,z)\
    *dn_dr(10.**x/rms,channel_width,z,params)
    output=integrate.quad(g,np.max([np.log10(smin),s_min]),np.log10(smax))[0]#dN/dr_com
    output=survey_area*output*C*1e-3/COSMO.Hz(z)*LITTLEH#dN/dz
    return output



def varFlux(s,k=3000.,s0=.88,gamma=-1.75,smin=0.):
    '''
    return variance in flux (Jy^2/sr) from dn/ds=k(s/s0)^gamma power law
    Args:
        k, power-law coefficient (Jy/Sr)
        s0, reference flux (Jy)
        gamma, power law index
        smin, minimum flux to integrate from
    '''
    return (k*s0**(gamma)*s**(-gamma+3))/(3-gamma)-(k*s0**(gamma)*smin**(-gamma+3))/(3-gamma)

def confusionLimit(sigPSF,k,s0,gamma,nsigma=5.):
    sAngle=2.*PI*sigPSF*sigPSF/(gamma-1)
    return (nsigma**(3-gamma)/(3-gamma))**(1/(gamma-1))*(k*s0**gamma*sAngle)**(1/(gamma-1))
def confusionVar(k,s0,gamma,sigPSF,nsigma=5):
    sConf=nsigma*confusionLimit(sigPSF,k,s0,gamma,nsigma)
    return sConf
#    return varFlux(sConf,k,s0,gamma),sConf
