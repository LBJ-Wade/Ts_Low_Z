'''
Indexed functions for computing the radio background
'''

import inspect
import radio_background

global RADIO_BACKGROUND_MODELS;RADIO_BACKGROUND_MODELS={}

for model in inspect.getmembers(radio_background,inspect.isfunction):
    RADIO_BACKGROUND_MODELS[model[0]]=model[1]
