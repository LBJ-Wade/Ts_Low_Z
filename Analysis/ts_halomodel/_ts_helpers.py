'''
functions that help with Tspin calculations.
'''
import inspect
import ts_helpers

global TS_HELPERS;TS_HELPERS={}
'''
populate TS_HELPERS
'''
for helper in inspect.getmembers(ts_helpers,inspect.isfunction):
    TS_HELPERS[helper[0]]=helper[1]
