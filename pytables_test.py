# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 13:33:29 2012

@author: deen
"""

import os
import time
import glob

import h5py
import numpy as np
from numpy import float64, array

from load_packing import getPackings 
from fs_tools import getPrefix
 
import tables
def tables_require_group(root, name, *args, **kwargs):
    h5f = root._v_file
    try:
        return root._f_getChild(name)
    except tables.exceptions.NoSuchNodeError:
        return h5f.createGroup(root, name, *args, **kwargs)

def createpyTablesNodeForPacking(root, packing):
    Ngroup = tables_require_group(root, "N%i" % packing['N'])
    Pgroup = tables_require_group(Ngroup, "p%.4e" % packing['P0'])
    
    pkgroup_name = "%04i" % packing['PackingNumber']
    if pkgroup_name in Pgroup:
        raise Exception("Packing already in hdf5 storage")
    pkgroup = tables_require_group(Pgroup, pkgroup_name)
    
    packing = packing.copy()  
    particles = packing.pop('particles')
    
    for key, value in packing.iteritems():
        pkgroup._v_attrs[key] = value
    
    root._v_file.createTable(pkgroup, 'particles', particles, expectedrows=max(particles.shape), chunkshape=particles.shape)
    
    return pkgroup._v_pathname


packing_attr_cache_dtype = np.dtype([
('N - Ncorrected', '<i4'),
('PackingNumber', '<i4'),
('P0', '<f8'),
('current time', '|S19'),
('H', '<f8'),
('L1', '<f8', (2,)),
('L', '<f8'),
('N', '<i4'),
('P', '<f8'),
('L2', '<f8', (2,)),
('sxy', '<f8'),
('delta', '<f8'),
('sxx', '<f8'),
('alpha', '<f8'),
('syy', '<f8'),
('runtime (s)', '<f8'),
('FIRE iteration count', '<i4'),
('Z', '<f8'),
('dH', '<f8'),
('phi', '<f8'),
('CG iteration count', '<i4'),
('dU', '<f8'),
('maxGrad', '<f8'),
('Uhelper', '<f8'),
('path', '|S128')])

def addTo_packing_attr_cache(cache, packing, path):
    packing = packing.copy()
    packing['path'] = path
    
    row = cache.row
    for key in cache.colnames:
        if key in packing:
            row[key] = packing[key]
    row.append()

def recursor(parent):
    yield parent
    try:
        for children in parent.itervalues():
            for item in recursor(children):
                yield item
    except AttributeError:
        pass

 
bbase = r"U:\novamaris\simu\Packings\N256"

print os.getcwd()

f = tables.openFile(getPrefix(bbase) + "_tables.h5", mode = "a")
root = f.root
try:
    attrcache = root.packing_attr_cache
except tables.exceptions.NoSuchNodeError:
    attrcache = f.createTable(root, 'packing_attr_cache', packing_attr_cache_dtype)
try:
    for base in glob.glob(bbase + "*"):
        for packing in getPackings(base):
            path = createpyTablesNodeForPacking(root, packing)
            addTo_packing_attr_cache(attrcache, packing, path)
finally:
    f.flush()
    f.close()    
    
    