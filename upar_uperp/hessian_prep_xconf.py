# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 12:07:12 2013

@author: Merlijn van Deen
"""
print "Importing packages: this takes *ages* on the cluster!"
import sys
sys.path.append('/home/merlijn/analysis/src/phdlib')
sys.path.append('C:/Users/Merlijn van Deen/Documents/GitHub/phd-library')

from packing_tools import make_shear_graphs
from packing_tools import load_packing
from packing_tools import V_harm
import os, glob
import time
import numpy as np


from hpc import HessianPackingCalculator

loadfolder = "/clusterdata/merlijn/fixedbox/"


import socket
if socket.gethostname() == "Isilwen":
    loadfolder = "d:/h5/xconf/"

storefolder = loadfolder + "uperp/"

def runner(fn):
    start = time.time()
    pack = load_packing.loadPacking(fn)
    N = pack['particles'].shape[0]
    
    if N != 1024:
        return
    else:
        print fn,
    
    Pstr = ("%.3e" % (pack["P0"])).replace(".", "").replace("000", "")
    npfn = ("N%i~P%s~" % (N,Pstr))+ os.path.split(fn)[1].replace(".txt", "") + "~FIXEDBOX~part_u_s.npz"

    try:
        globb = storefolder + "bad/*~%s~*" % os.path.split(fn)[1].split(".")[0]
        print globb,
        lfn = glob.glob(globb)[0]
        data = np.load(lfn)
        K0_inv = data["K0_inv"]
        print " cache loaded"
    except IndexError, e:
        print " cache not loaded (%r)" % e
        K0_inv = None
#    print "Plotting..."
#    make_shear_graphs.doPackingStuff(pack)
#    break
    #
    #print "Calculating G et al..."
    #conts = V_harm.get_contacts(pack)
    #K = V_harm.Hess(conts, pack)
    #elc = V_harm.El_Con(K, conts['rattlers'], pack)
    #V_harm.Pmod(elc)

    HPC = HessianPackingCalculator(None, pack, K0_inv=K0_inv)
    u_parr, u_parr_scaled, u_perp, u_perp_scaled = HPC.get_scaled_uparperps()
    G = HPC.deform_packing()["energy_cost"]
    print "Energy cost in shear: ", G
    
    print "saving %s to %s... runtime HPC: %f s" % (fn, npfn, time.time() - start)
    np.savez(storefolder + npfn,
             u_parr = u_parr, u_parr_scaled = u_parr_scaled,
             u_perp = u_perp, u_perp_scaled = u_perp_scaled,
             K0_inv = HPC.K0_inv,
             rij = HPC.contacts['rij'], dij = HPC.contacts['dij'], G=G)
             
print "Starting run / starting glob"
files = glob.glob(loadfolder + "Xconf*.txt")
p = None
if len(sys.argv) > 1 and sys.argv[1] == "-parr":
    import multiprocessing
    nprocs = multiprocessing.cpu_count() - 1
    print "Starting multiprocess pool on %i CPUs" % (nprocs)
    p = multiprocessing.Pool(processes=nprocs)
    map = p.map
map(runner, files)

if p:
    p.close()