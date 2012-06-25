# -*- coding: utf-8 -*-
# Created on Fri Jun 15 12:25:37 2012
# @author: Merlijn van Deen <deen@physics.leidenuniv.nl>

"""
Functions to load packing files.

Supported formats:
    - Jo's format, single packing
    - Jo's format, multiple packings (e.g. shear simulations)
    

Jo's format
-----------
Example:
    
```
N = 16 ,L = 9.5144469664731326 ,L1= { 9.8178283514271385 , 0.0000000000000000 }  ,L2= { -0.5721514491604729 , 9.2204403904322638 }  ,P = 0.0000001000000000 ,P0= 0.0000001000000000 ,
{
9.2855569466200685 ,	4.6248765843676147 ,	1.0000000000000000 ,	
( ... more (x, y, r) lines ... )
1.7389654101234185 ,	3.8493233607950302 ,	1.3999999999999999 ,	
}
```
"""

import os
import numpy as np
from numpy import float64

from load_log import getUniqueParsedLogLines
from fs_tools import getPrefix

def loadPackings(filename):
    try:
        f = open(filename)
    except Exception:
        import time
        time.sleep(30)
        f = open(filename)
        
    raw = f.read().strip()
    f.close()
    
    packings = raw.split("\n\n")
    return [loadPackingData(p) for p in packings]
    
def loadPacking(filename):
    data = loadPackings(filename)
    assert(len(data) == 1)
    return data[0]
    
def loadPackingData(raw):
    retval = {}
    
    raw_as_list = "[" + raw.replace("=", ",").replace("{", "(").replace("}", ")").replace("\r", "") + "]"
    variables = ['N', 'L', 'L1', 'L2', 'P', 'P0']
    f = float64
    dtypes    = [int, f  , f   , f   , f  , f   ]
    
    data = eval(raw_as_list, dict((a,a) for a in variables))
    # data = ["N", 64, "L1", [...], ..., [particle positions]]        
    
    # first process the variables
    vardata = data[:-1] # ["N", 64, ..., "P0", 1e-7]
    keys = vardata[0::2] # ["N", ..., "P0"]
    values = vardata[1::2] # [64, ..., 1e-7]
    
    for key, dtype, value in zip(keys, dtypes, values):
        retval[key] = dtype(value)
                
    # now process the particles
    particles = np.array(data[-1]) # [x1, y1, r1, x2, ...]
    particles = particles.reshape((particles.size/3, 3)) # [[x1, y1, r1], [x2, ...]] 
    particles = particles.transpose()
    
    particles.dtype = [('x', float64), ('y', float64), ('r', float64)]
    
    particles = particles[0] # hack around structured array creation; for some
                             # reason there is a wrapping [ ]; this removes that.
    
    #particles.dtype.names = ['x', 'y', 'r']
    retval['particles'] = particles
    
    return retval
   

class PackingWarning(UserWarning):
    pass

def getPackings(folder):
    prefix = getPrefix(folder)
    for logline in getUniqueParsedLogLines(folder):
        try:
            packingfn = "%s~%04i.txt" % (prefix, logline['PackingNumber'])
            packing = loadPacking(os.path.join(folder, packingfn))
            packing.update(logline) # copy information from log line
            yield packing
        except Exception, e:
            import warnings
            warnings.simplefilter("always", PackingWarning)
            warnings.warn(str(e), PackingWarning)

def getShear(folder):
    prefix = getPrefix(folder)