# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 13:33:29 2012

@author: deen
"""

import os, sys
import time
import glob
import itertools
import pandas
import numpy as np
from numpy import float64, array, dtype

from load_packing import loadPackings 
from fs_tools import getPrefix
 
import tables
from pytables_tools import require_group, require_table, add_to_table, store_table
from pytables_test import create_group_for_packing, insert_packing_parameters, insert_packing_particles

packing_attr_cache_dtype = np.dtype([
('L', dtype('float64')),
('L1', dtype('float64'), (2,)),
('L2', dtype('float64'), (2,)),
('N', dtype('int64')),
('P', dtype('float64')),
('P0', dtype('float64')),
('gamma', dtype('float64')),
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

def shear_type(st):
    if st.startswith("SR"):
        return "SR"
    else:
        return st

def groupkey(args):
    N,P,st,step,number = args[0]
    return N,P,number,shear_type(st)
    
def sortkey(args):
    return list(groupkey(args)) + args[0]

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

def process_measurement(f, key, base, m, spec):
    group = require_group(f,'/'.join(key))
    
#    spec = "~".join(m) + ".txt"
    print key, ":", spec,
    
    log, comments = read_csv(os.path.join(base, "log" + spec))
    print "log",
    sys.stdout.flush()
    
    data, comments = read_csv(os.path.join(base, "data" + spec))
    print "data",
    sys.stdout.flush()


    # 'eta' (applied strain) was renamed to gamma
    if 'eta' in data and 'gamma' not in data:
        data['gamma'] = data['eta']
        data.pop('eta')
        
    data["gamma"][log["step#"] == 0] = 0 # work around gamma=10^-9 bug
    
    group._v_attrs['comments'] = comments

    print "particles",
    sys.stdout.flush()
    
    try:
        open(os.path.join(base, "particles" + spec)).close()
	particles = loadPackings(os.path.join(base, "particles" + spec))
	particles = [particles[0]] + particles[1::2]
        particles = pandas.DataFrame(particles)
	print len(particles)
	print len(data)
        packings = data.join(particles, rsuffix="_").join(log, rsuffix="__")
	print len(packings)
    except IOError, e:
        print "(failed: %r)" % e
        sys.stdout.flush()
        packings = data.join(log, rsuffix="__")

    attrcache = require_table(group, 'data', packing_attr_cache_dtype,
                              expectedrows=len(log), chunkshape=(len(log),))
    for rowno, row in packings.iterrows():
        subgroup = require_group(group, "%04i" % row["step#"])
        
        insert_packing_parameters(subgroup, row)
        add_to_table(attrcache, row, path=subgroup._v_pathname)
        
        if 'particles' in packings:
            insert_packing_particles(subgroup, row)
    
    attrcache.flush()
    f.flush()
        

def key_for_fn(fn):
  part1 = os.path.split(os.path.split(fn)[0])[1].split("~") 
  part2 = os.path.split(fn)[1][4:-4].split("~")[2:]
  return list(part1) + list(part2)


try:
    shear_measurements = [(key_for_fn(fn), fn)
                            for fn
                            in glob.glob(bbase + "*/data*.txt")]
    for key, measurements in itertools.groupby(sorted(shear_measurements, key=sortkey), groupkey):
        for m in measurements:
            pass
        print key
	try:
	  print root, key, os.path.split(m[1])[0], m[0]
          process_measurement(root, key, os.path.split(m[1])[0], m[0], os.path.split(m[1])[1][4:])        
	except Exception, e:
	  raise
	  print key, e
finally:
    f.flush()
    f.close()
print "rrr"
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
