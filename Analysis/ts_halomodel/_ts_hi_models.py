'''
models describing rho_HI/Ts distribution in dark-matter halo
'''
import inspect
import ts_hi_models
'''
TS_MODELS is a dictionary with string keys referring to functions
describing the distribution of (rho_HI/spin-temperature)^\alpha
 in dark-matter halos in both real and Fourier space.
'''

global TS_HI_MODELS;TS_HI_MODELS={}
'''
populate TS_HI_MODELS
'''
for model in inspect.getmembers(ts_hi_models,inspect.isfunction):
    TS_HI_MODELS[model[0]]=model[1]
