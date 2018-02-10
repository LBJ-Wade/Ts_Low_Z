'''
file with models describing HI distribution within a dark matter halo.
'''
import inspect
import hi_models

'''
HI_MODELS is a dictionary with strings for keys
that contains HI density functions for dark matter halos
in both real and Fourier space.
The format of a density function should be
rho(r (Mpc/h), virial mass (msolar/h), redshift,model_parameters)
OR
rho(k (h/Mpc), virial mass (msolar/h), redshift, model_parameters)

where model_parameters is a dictionary of parameter names
and values.
'''
global HI_MODELS;HI_MODELS={}
'''
populate HI_MODELS
'''
for model in inspect.getmembers(hi_models,inspect.isfunction):
    HI_MODELS[model[0]]=model[1]
