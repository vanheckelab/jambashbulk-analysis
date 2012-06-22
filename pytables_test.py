# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 13:33:29 2012

@author: deen
"""

import warnings
from tables.path import NaturalNameWarning
warnings.simplefilter('ignore', NaturalNameWarning)

import os
import time
import glob

import numpy as np
from numpy import float64, array

from load_packing import getPackings 
from fs_tools import getPrefix
 
import tables
from pytables_tools import require_group, require_table, add_to_table, store_table

def create_group_for_packing(root, packing):
    Ngroup = require_group(root, "N%i" % packing['N'])
    Pgroup = require_group(Ngroup, "p%.4e" % packing['P0'])
    
    pkgroup_name = "%04i" % packing['PackingNumber']
    if pkgroup_name in Pgroup:
        raise Exception("Packing already in hdf5 storage")
    pkgroup = require_group(Pgroup, pkgroup_name)  
    
    return pkgroup

def insert_packing(node, packing):
    packing = packing.copy()  
    particles = packing.pop('particles')
    
    for key, value in packing.iteritems():
        node._v_attrs[key] = value
    
    store_table(node, 'particles', particles)    

def createpyTablesNodeForPacking(root, packing):
    pkgroup = create_group_for_packing(root, packing)
    insert_packing(pkgroup, packing)
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

import sys
bbase = sys.argv[1] #r"/home/merlijn/simu/Packings/N1024"

print os.getcwd()

ff = os.path.join('/data/misc/granulargroup/h5/', getPrefix(bbase) + "_tables.h5")

f = tables.openFile(ff, mode = "a")
root = f.root

attrcache = require_table(root, 'packing_attr_cache', packing_attr_cache_dtype)
try:
    for base in glob.glob(bbase + "*"):
        for packing in getPackings(base):
            try:
                path = createpyTablesNodeForPacking(root, packing)
		print path
                add_to_table(attrcache, packing, path=path)
            except Exception, e:
                warnings.warn(str(e), stacklevel=5)
finally:
    f.flush()
    f.close()    
    
    
