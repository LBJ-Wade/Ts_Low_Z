'''
This file indexes functions for computing spin temperature
'''
import ts_models
import inspect

global TS_MODELS;TS_MODELS={}
for model in inspect.getmembers(ts_models,inspect.isfunction):
    TS_MODELS[model[0]]=model[1]
