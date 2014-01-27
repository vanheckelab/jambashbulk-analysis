# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 13:33:29 2012

@author: deen
"""

import warnings
from tables.path import NaturalNameWarning
warnings.simplefilter('ignore', NaturalNameWarning)

import sys
import os
sys.path.append(os.path.split(os.path.split(__file__)[0])[0])
import os
import time
import glob

import numpy as np
from numpy import float64, array
try:
    from numpy import float128 as bigfloat
except Exception:
    from numpy import float64 as bigfloat
    
from packing_tools.load_packing import getPackings 
 
import tables
from pytables_tools import require_group, require_table, add_to_table, store_table

def create_group_for_packing(root, packing):
    Ngroup = require_group(root, "N%i" % packing['N'])
    Pgroup = require_group(Ngroup, "p%.4e" % packing['P0'])
    
    pkgroup_name = "%04i" % packing['PackingNumber']
    if pkgroup_name in Pgroup:
        #raise Exception("Packing already in hdf5 storage")
        return None
    pkgroup = require_group(Pgroup, pkgroup_name)  
    
    return pkgroup

def insert_attribute(node, key, value):
    if getattr(value, 'dtype', type(value)) == bigfloat:
        #long double. split in two entries
        major = np.float64(value)
        minor = np.float64(value - bigfloat(major))
        insert_attribute(node, key, major)
        insert_attribute(node, key+"_err", minor)
        return
    node._v_attrs[key] = value

def insert_packing_parameters(node, packing):
    for key, value in packing.iteritems():
        if key == 'particles':
            continue
        insert_attribute(node, key, value)

def insert_packing_particles(node, packing):
    store_table(node, 'particles', packing['particles'])   
    
def insert_packing(node, packing):
    insert_packing_parameters(node, packing)
    insert_packing_particles(node, packing)

def createpyTablesNodeForPacking(root, packing):
    pkgroup = create_group_for_packing(root, packing)
    if pkgroup:
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
('Z_calc', '<f8'),
('path', '|S128')])


def main(*argv):
    if len(argv) != 2:
        print argv
        raise Exception("Usage: %s <glob base> <output hdf file>")

    bbase = argv[0]

    f = tables.openFile(argv[1], mode = "a")
    root = f.root

    attrcache = require_table(root, 'packing_attr_cache', packing_attr_cache_dtype)
    try:
        for base in glob.glob(bbase + "*"):
            print base
            for packing in getPackings(base):
                print packing['N'], packing['P'], packing['PackingNumber']
                path = createpyTablesNodeForPacking(root, packing)
                if path:
                    add_to_table(attrcache, packing, path=path)
    finally:
        f.flush()
        f.close()

if __name__=="__main__":
    main(*sys.argv[1:])
