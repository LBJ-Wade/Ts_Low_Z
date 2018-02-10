'''
utils.py utility functions
'''

from settings import COSMO
from colossus.lss import mass_function
from colossus.lss import bias as col_bias
import numpy as np
import pyfftlog

#*********************************************************
#utility functions
#*********************************************************
def dict2tuple(dictionary):
    '''
    Convert a dictionary to a tuple of alternating key values
    and referenced values.
    Example: {'apple':0,'beta':2000} -> ('apple',0,'beta',20000)
    Args:
        dict, dictionary
    Returns:
        tuple of (key0, value0, key1, value1, key2, value2, ...)
    '''
    output=()
    for key in dictionary.keys():
        if isinstance(dictionary[key],dict):
            output=output+dict2tuple(dictionary[key])
        else:
            output=output+(key,dictionary[key])
    return output
#**********************************************************
#cosmology utility functions
#**********************************************************
def massfunc(m,z,model='tinker08',mdef='200m'):
    '''
    dark matter mass function dN/dlog_10M/dV (h^3 Mpc^-3)
    Args:
        mass, float, msolar/h
        z, float, redshift
        model, string, dark matter prescription (see colossus documentation)
        mdef, string, mass definition.
    Returns:
        dark matter mass function dN/dlog_10M/dV (h^3 Mpc^-3)
    '''
    return mass_function.massFunction(m,z,mdef=mdef,model=model,
    q_out='dndlnM')*np.log(10.)


def bias(mvir,z,mode='tinker10',mdef='200m'):
    return col_bias.haloBias(mvir, model = 'tinker10', z = z,mdef=mdef)


def fftlogbessel(func,npts,mini,maxi,tdir=1,dim=2):
    '''
    Bessel function transform.
    Uses pyfftlog
    Args:
        func, function from R1->R1
         npts, int, number of points to evaluate the bessel function transform
         mini,float, minimum log_10 value for input (i)
         maxi,float, maximum log_10 value for output (i)
         tdir,int, direction of the transform (forward or backward)
         dim,int, 2d or 3d transform
    Returns:
        oaxis, npts numpy array,  fourier-space axis.
        ao, npts numpy array, log-fft'd function function values
    '''
    assert dim==2 or dim ==3
    assert tdir==1 or tdir==-1
    if tdir==1:
        onormal=1.
    elif tdir==-1:
        onormal=2.*PI
    if dim==3:
        mu=.5
    else:
        mu=0.
    logic=(mini+maxi)/2.
    nc=(npts+1)/2.
    dlogi=(maxi-mini)*1./npts
    dlni=dlogi*np.log(10.)
    iaxis=10.**(logic+(np.arange(1,npts+1)-nc)*dlogi)
    kr,xsave=pyfftlog.fhti(npts,
                           mu,
                           dlni,
                           0,
                           1.,
                           1)
    logoc=np.log10(kr)-logic
    oaxis=10.**(logoc+(np.arange(1,npts+1)-nc)*dlogi)
    ai=func(iaxis)*iaxis
    if dim==3:
        ai*=2.*iaxis*np.sqrt(PI/2./iaxis)
    #print ai
    ao=2.*PI*pyfftlog.fht(ai.copy(),xsave,tdir=1)/onormal**2.
    if dim==3:
        ao/=oaxis**.5
    ao=ao/oaxis
    #print ao
    return oaxis,ao
