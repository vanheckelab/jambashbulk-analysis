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

from load_packing import loadPacking
from load_log import getUniqueParsedLogLines
from fs_tools import getPrefix

        
def getPackings(folder):
    prefix = getPrefix(folder)
    for logline in getUniqueParsedLogLines(folder):
        packingfn = "%s~%04i.txt" % (prefix, logline['PackingNumber'])
        packing = loadPacking(os.path.join(folder, packingfn))
        packing.update(logline) # copy information from log line
        yield packing
 
def createh5pyNodeForPacking(root, packing):
    Ngroup = root.require_group("N%i" % packing['N'])
    Pgroup = Ngroup.require_group("p%.4e" % packing['P0'])
    
    pkgroup_name = "%04i" % packing['PackingNumber']
    if pkgroup_name in Pgroup:
        raise Exception("Packing already in hdf5 storage")
    pkgroup = Pgroup.require_group(pkgroup_name)
    
    packing = packing.copy()  
    particles = packing.pop('particles')
    
    for key, value in packing.iteritems():
        pkgroup.attrs[key] = value

    pkgroup.create_dataset('particles', data=particles)
    #for dim, values in particles.iteritems():
    #    pkgroup.create_dataset(dim, data=values)

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

run_h5py = False
if run_h5py:
    start = time.time()
    f = h5py.File(getPrefix(bbase) + ".h5",'a')
    try:
        for base in glob.glob(bbase + "*"):
            for packing in getPackings(base):
                createh5pyNodeForPacking(f, packing)
    finally:
        f.flush()
        f.close()
    print "insertion: %i s" % (time.time() - start)
    
    start = time.time()
    f = h5py.File(getPrefix(bbase) + ".h5",'a')
    Z = filter(None, [item.attrs.get('Z') for item in recursor(f)])
    print "Getting all Z: %i s" % (time.time() - start)
    print np.median(Z)
    
else:
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
    
    