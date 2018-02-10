'''
functions that help with HI distributions.
'''
import inspect
import hi_helpers

global HI_HELPERS;HI_HELPERS={}
'''
populate HI_HELPERS
'''
for helper in inspect.getmembers(hi_helpers,inspect.isfunction):
    HI_HELPERS[helper[0]]=helper[1]
