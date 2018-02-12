'''
settings.py
This file is where global variables governing interpolation ranges and physical
constants are defines.
'''
from colossus.cosmology import cosmology
import numpy as np
'''
Physical constants
'''
global PI;PI=np.pi
global MSOL;MSOL=1.99e30#solar mass in kg
global MP;MP=1.67e-27#proton mass in kg
global KPC;KPC=3.068e19#kpc in meters
global LSOL;LSOL=3.828e26#solar luminosity in watts.
global C;C=299792458.#the speed of light in m/s
global G;G=6.67e-11#Gravitational constant. kg^-1 m^3/s^2
global GC;GC=G*MSOL/(1e3*KPC)**3.#G in Mpc^3/solar masses/sec^2
global YEAR;YEAR=365.*24.*3600. #1 year in seconds
global KBOLTZMANN;KBOLTZMANN=1.28e-23#Boltzmann Constant in Joules/Kelvin
global TCMB;TCMB=2.73#CMB temperature in Kelvin
global JY;JY=1e-26#1 Jy in Watts/m^2/Hz
global F21;F21=1420405751.7667 #HI hyperfine line frequency (Hz)
global A10;A10=2.85e-15#Hz (Einstein A coeff for HI hyerfine transition)
global HPLANCK;HPLANCK=6.62607004e-34# Joule Seconds (Planck constant)
global YP;YP=0.24
global NHIDT2TAU
NHIDT2TAU=3.*C*C*A10*HPLANCK/32./PI/KBOLTZMANN/F21 #factor converting
                                                          #from NHI/T*phi(nu)
                                                          #to optical depth
ERG=1e-7#ERG in Joules
'''
Cosmological constants
'''
global COSMOSTR;COSMOSTR='planck15'
global COSMO;COSMO=cosmology.setCosmology(COSMOSTR)
LITTLEH=COSMO.H0/100.#dimensionless hubble constant
'''
interpolation settings
'''
global KMIN;KMIN=-5#minimum k for interpolation
global KMAX;KMAX=5#maximum k for interpolation
global MINLOG;MINLOG=-3
global MAXLOG;MAXLOG=3
global NINTERP_KPARA;NINTERP_KPARA=30
global NINTERP_KPERP;NINTERP_KPERP=30
global KPERP_INTERP_MIN;KPERP_INTERP_MIN=-5
global KPERP_INTERP_MAX;KPERP_INTERP_MAX=3
global KPARA_INTER_MIN;KPARA_INTERP_MIN=-5
global KPARA_INTERP_MAX;KPARA_INTERP_MAX=4
global KPERP_FOREST_MIN;KPERP_FOREST_MIN=-3
global KPERP_FOREST_MAX;KPERP_FOREST_MAX=3
global NINTERP_Z;NINTERP_Z=100
global MAXZ;MAXZ=8.
global M_INTERP_MIN;M_INTERP_MIN=6.
global M_INTERP_MAX;M_INTERP_MAX=15.
global N_INTERP_M;N_INTERP_M=100
global SNR_INTERP_MIN; SNR_INTERP_MIN=np.log10(5.)
global SNR_INTERP_MAX; SNR_INTERP_MAX=12.
global N_INTERP_SNR;N_INTERP_SNR=100
global TAU_INTERP_MIN; TAU_INTERP_MIN=-12.
global TAU_INTERP_MAX; TAU_INTERP_MAX=0.
global N_INTERP_TAU;N_INTERP_TAU=100
global R_INTERP_MIN;R_INTERP_MIN=-9
global R_INTERP_MAX;R_INTERP_MAX=5
global N_INTERP_RTAU; N_INTERP_RTAU=1000
global S_INTERP_MIN; S_INTERP_MIN=-6
global S_INTERP_MAX; S_INTERP_MAX=6.
global NINTERP_S; NINTERP_S=100
'''
integration settings.
'''
global EPSABS;EPSABS=1e-8
global EPSREL;EPSREL=1e-8
global NPTS;NPTS=256
'''
global SPLINE dictionary
'''
global SPLINE_DICT;SPLINE_DICT={}
