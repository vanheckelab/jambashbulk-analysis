# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 13:33:29 2012

@author: deen
"""

import os
import time
import glob
import itertools
import pandas
import h5py
import numpy as np
from numpy import float64, array

from load_packing import loadPackings 
from fs_tools import getPrefix
 
import tables
from pytables_tools import require_group, require_table, add_to_table, store_table
from pytables_test import create_group_for_packing, insert_packing

packing_attr_cache_dtype = np.dtype([
('L', dtype('float64')),
('L1', dtype('float64'), (2,)),
('L2', dtype('float64'), (2,)),
('N', dtype('int64')),
('P', dtype('float64')),
('P0', dtype('float64')),
('eta', dtype('float64')),
('s_xy', dtype('float64')),
('Ncontacts', dtype('int64')),
('Nchanges', dtype('int64')),
('N+', dtype('int64')),
('N-', dtype('int64')),
('Z', dtype('float64')),
('step#', dtype('int64')),
('alpha', dtype('float64')),
('delta', dtype('float64')),
('phi', dtype('float64')),
('#rattler', dtype('int64')),
('s_xx', dtype('float64')),
('s_yy', dtype('float64')),
('U', dtype('float64')),
('dU', dtype('float64')),
('H', dtype('float64')),
('dH', dtype('float64')),
('t_run', dtype('int64')),
('#FIRE', dtype('int64')),
('#CG', dtype('int64')),
('gg', dtype('float64')),
('creation-date', dtype('|S19')),
('path', '|S128')])

if len(sys.argv) != 3:
    print "Usage: %s <base from which to add shear data from> <h5 file to store shear data in>"
    exit(1)

bbase = sys.argv[1] #r"U:\ilorentz\simulations\Packings\N256~P1e-3"
outfile = sys.argv[2] #r"U:\ilorentz\simulations\Packings\N256~P1e-3\data.h5"

f = tables.openFile(outfile, mode = "a")
root = f.root

attrcache = require_table(root, 'shear_attr_cache', packing_attr_cache_dtype)

def shear_type(st):
    if st.startswith("SR"):
        return "SR"
    else:
        return st

def groupkey(args):
    N,P,st,step,number = args
    return N,P,number,shear_type(st)
    
def sortkey(args):
    return list(groupkey(args)) + args

def read_csv(fn):
    f = open(fn)
    comments = ""
    
    while True:
        oldpos = f.tell()
        line = f.readline()
        
        if line.startswith('#'):
            comments += line
        else:
            f.seek(oldpos)
            break
    return pandas.read_csv(f, sep="[ \t]"), comments

def process_measurement(f, attrcache, key, base, m):
    
    group = require_group(f,'/'.join(key))
    
    spec = "~".join(m) + ".txt"
    print key, ":", spec
    
    log, comments = read_csv(os.path.join(base, "log" + spec))
    data, comments = read_csv(os.path.join(base, "data" + spec))
    group._v_attrs['comments'] = comments
    particles = pandas.DataFrame(loadPackings(os.path.join(base, "particles" + spec)))

    packings = data.join(particles, rsuffix="_").join(log, rsuffix="__")

    for rowno, row in packings.iterrows():
        subgroup = require_group(group, "%04i" % rowno)
        path = insert_packing(subgroup, row)
        add_to_table(attrcache, row, path=subgroup._v_pathname)

try:
    shear_measurements = [os.path.split(fn)[1][4:-4].split("~")
                            for fn
                            in glob.glob(bbase + r"\data*.txt")]
    
    for key, measurements in itertools.groupby(sorted(shear_measurements, key=sortkey), groupkey):
        for m in measurements:
            pass
        
        process_measurement(root, attrcache, key, bbase, m)        
finally:
    f.flush()
    f.close()
#try:
#    for base in glob.glob(bbase + "*"):
#        for packing in getPackings(base):
#            path = createpyTablesNodeForPacking(root, packing)
#            add_to_table(attrcache, packing, path=path)
#finally:
#    f.flush()
#    f.close()    
#    
#    