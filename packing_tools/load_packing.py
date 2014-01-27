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
import re
import numpy as np
import itertools
#import jamBashbulk

from numpy import float64, array, float64 as float128
from numpy import loadtxt
from load_log import getUniqueParsedLogLines
from fs_tools import getPrefix

if os.name == 'posix':
    from parser import parser
    def loadPackings(filename):
        return [p for p in parser.read_packings(filename)]
else:
    def loadPackings(filename):
        try:
            f = open(filename)
        except Exception:
            import time
            time.sleep(30)
            f = open(filename)
            
        raw = f.read().strip().lstrip('\x00')
        f.close()
        
        packings = raw.split("\n\n")
        return [loadPackingData(p) for p in packings]
        
def loadPacking(filename):
    if os.path.split(filename)[1].lower().startswith('xconf'):
        return read_xconf(filename)
    data = loadPackings(filename)
    assert(len(data) == 1)
    return data[0]

def savePacking(packing, filename, mode='w'):
    f = open(filename, mode)
    f.write(strPacking(packing))
    f.close()

def strPacking(packing):
    ret = []
    ret.append("N = {N} ,L = {L:.16f} ,L1= {{ {L1[0]:.16f} , {L1[1]:.16f} }}  ,L2= {{ {L2[0]:.16f} , {L1[0]:.16f} }}  ,P = {P:.16f} ,P0= {P0:.16f} ,".format(**packing))
    ret.append("{")
    for x, dx, y, dy, r in packing["particles"]:
        ret.append("{x:.16f} ,\t{y:.16f} ,\t{r:.16f} ,\t".format(x=x+dx, y=y+dy, r=r))
    ret.append("}")
    ret.append("")
    ret.append("")
    return "\n".join(ret)

def mirror(packing, x=True, y=False):
    if y:
        raise Exception('not implemented')

    import copy
    packing = copy.deepcopy(packing)
    
    if x:
        # flip alpha
        packing['L2'][0] = -packing['L2'][0]   
        
        # flip X values
        packing['particles']['x'] = -packing['particles']['x']
        packing['particles']['x_err'] = -packing['particles']['x_err']

        # shift particles back into unit cell
        packing['particles']['x'] += packing['L1'][0]

    return packing
def read_xconf(filename):
    """Load Packing from the ancient XConf format"""
    f = open(filename)
    precomments = list(itertools.takewhile(lambda x: x.startswith("#"), f))
    lx, ly = array([[float(fl) for fl in x.split("=")[1].split()] for x in precomments[-2:]])
    data = loadtxt(itertools.takewhile(lambda x: not x.startswith("#"), f))
    postcomments = list(itertools.takewhile(lambda x: x.startswith("#"), f))
    P0 = float(postcomments[0].split("=")[1].strip().split()[0])
    
    particles = data.view([('r', '<f8'), ('x', '<f8'), ('y', '<f8')])[:,0]
    particles['x'], particles['y'] = (particles['x'] * lx[:,np.newaxis] + particles['y'] * ly[:,np.newaxis])
    
    
    f.close()    
    packing = {'L1': lx, 'L2': ly, 'P0': P0,
               'particles': particles}

    return packing

def parse_longdouble(value):
    import ctypes
    from ctypes import c_longdouble, pointer

    val = c_longdouble()
    val_ptr = pointer(val)
    try:
        ctypes.cdll.LoadLibrary("libc.so.6").sscanf(value, "%Le", val_ptr)
    except WindowsError:
        return float(value)

    return np.ctypeslib.as_array(val_ptr, (1,))[0]
    
def loadPackingData(raw):
    retval = {}
    raw_as_list = "[" + raw.replace("=", ",").replace("{", "(").replace("}", ")").replace("\r", "") + "]"
    raw_as_list = re.sub("([0-9]+\.[0-9]+)", r"ld('\1')", raw_as_list)

    variables = ['N', 'L', 'L1', 'L2', 'P', 'P0']
    f = np.longdouble
    dtypes    = [np.int32, f  , f   , f   , f  , f   ]
    
    eval_environ = dict((a,a) for a in variables)
    eval_environ['ld'] = parse_longdouble
    data = eval(raw_as_list, eval_environ)
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
    x = particles[0]
    y = particles[1]
    r = particles[2]

    #convertedVectors = jamBashbulk.convertLvectors(retval["L1"], retval["L2"])
    #retval['alpha'] = convertedVectors['alpha']
    #retval['delta'] = convertedVectors['delta']
    #for key, value in jamBashbulk.get_packing_data(retval["N"], retval["P0"], x, y, r, **convertedVectors).iteritems():
    #    retval[key+"_calc"] = value


    x_major = float64(x); x_minor = float64(x-float128(x_major))
    y_major = float64(y); y_minor = float64(y-float128(y_major))
    r = float64(r)

    particles = array(zip(x_major, x_minor, y_major, y_minor, r), dtype=[('x', float64), ('x_err', float64),
                                                                         ('y', float64), ('y_err', float64),
                                                                         ('r', float64)])

#    particles = particles[0] # hack around structured array creation; for some
                             # reason there is a wrapping [ ]; this removes that.
   
    # split particles in two float64 parts
    #particles.dtype.names = ['x', 'y', 'r']
    retval['particles'] = particles
    
    return retval
   

class PackingWarning(UserWarning):
    pass

import tempfile
import shutil
import subprocess
def mkPackingLog(packingFile):
    path, fn = os.path.split(packingFile)
    N,P,num = fn.split(".")[0].split("~")

    tmpdir = tempfile.mkdtemp(prefix='mkPacking-Merlijn')
    try:
        os.makedirs(tmpdir + "/Packings/%s~%s" % (N,P))
        shutil.copy(packingFile, tmpdir + "/Packings/%s~%s/" % (N,P) + fn)
        j2d = subprocess.Popen(["/home/valhallasw/src/jambashbulk/bin/jam2D"], stdin=subprocess.PIPE, cwd=tmpdir)
        j2d.communicate("\n".join(["Packings/%s~%s" % (N,P), "1", "1", num, "3"]))
        shutil.copy(tmpdir + "/Packings/%s~%s/" % (N,P) + "log"+fn, path + "/log" + fn)
    finally:
        shutil.rmtree(tmpdir)
        pass

def getPackings(folder):
    import warnings
    warnings.simplefilter("always", PackingWarning)
    prefix = getPrefix(folder)
    import glob
    packingfns = set(glob.glob(os.path.join(folder, "%s~????.txt") % prefix))
    loglines = list(getUniqueParsedLogLines(folder))
    packingswithlogs = set(os.path.join(folder, packingfn) for packingfn in ("%s~%04i.txt" % (prefix, logline['PackingNumber']) for logline in loglines))
    missing = packingfns - packingswithlogs

    for logline in getUniqueParsedLogLines(folder):
        try:
            packingfn = "%s~%04i.txt" % (prefix, logline['PackingNumber'])
            fullfn = os.path.join(folder, packingfn)
            try:
                packingfns.remove(fullfn)
            except KeyError, e:
                warnings.warn(repr(e), PackingWarning)
            packing = loadPacking(fullfn)
            packing.update(logline) # copy information from log line
            yield packing
        except Exception, e:
            warnings.warn(repr(e), PackingWarning)

    if len(missing) != 0:
         for fn in missing:
             packing = loadPacking(fn)
             packing["PackingNumber"] = int(os.path.split(fn)[1].split(".")[0].split("~")[-1])
             import jamBashbulk
             packing = jamBashbulk.updatePackingParams(packing)
             yield packing
             #print "Generating logfile for %s..." % fn
             #mkPackingLog(fn)
 

def getShear(folder):
    prefix = getPrefix(folder)
