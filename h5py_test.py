# -*- coding: utf-8 -*-
# Created on Fri Jun 15 13:23:01 2012
# @author: Merlijn van Deen <deen@physics.leidenuniv.nl>

"""
Basic example of how to use h5py to save packing files.
"""
import glob
import time
import numpy as np

import h5py

from fs_tools import getPrefix
from load_packing import getPackings

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


bbase = r"U:\novamaris\simu\Packings\N256"

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