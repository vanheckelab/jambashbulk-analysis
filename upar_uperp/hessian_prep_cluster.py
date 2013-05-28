# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 13:51:38 2013

@author: Merlijn van Deen


TO DO checks:
    
    1) is simple shear hetzelfde als pure shear (misschien met een factor 2 ertussen ofzo)
    1b) maakt de magnitude echt niet uit? (het is linalg, dus je mag hopen van niet...)
    2) zijn de krachten op de deeltjes in de eindsituatie echt nul?
    3) wat is de bewegingsvector dan? plot dit...
"""
from __future__ import division
from pylab import shuffle
import sys, glob, os
sys.path.append('/home/merlijn/analysis/src/phdlib')
import tables
import itertools
import numpy as np


path = "/clusterdata/merlijn/h5/"

from hpc import HessianPackingCalculator
import time

def worker((fn, num)):
    try:
        f = tables.File(fn)
        groups = iter(iter(iter(f).next()).next()).next()
        itr = iter(groups)
        for i in range(num+1):
            group = itr.next()
    
        start = time.time()
        npfn = path + os.path.split(fn)[1].split(".")[0] + "~" + group._v_name + "~part_u_s"
        print group._v_pathname, "; ", npfn, start
        sys.stdout.flush()
        
        
        if group._v_pathname in ["/N512/p1.0000e-06/9050"]:
            print "<SKIP DUE TO ERROR>"
            return
    
        HPC = HessianPackingCalculator(group)
        
        for deformation in [2,3]:
            u_parr, u_parr_scaled, u_perp, u_perp_scaled = HPC.get_scaled_uparperps(deformation=2)
        
            np.savez(npfn + str(deformation) + ".npz",
                     u_parr = u_parr, u_parr_scaled = u_parr_scaled,
                     u_perp = u_perp, u_perp_scaled = u_perp_scaled,
                     rij = HPC.contacts['rij'], dij = HPC.contacts['dij'])
                    
        print "done, runtime: ", time.time()-start
    except Exception, e:
        print "FFFFUUUUU", e

if __name__ == "__main__":
    import sys
    
    paths = []
    print "BASES: ", sys.argv[1:]
    for base in sys.argv[1:]:
        print "PROCESSING BASE " + path + base + "_tables/*.h5"
        paths += glob.glob(path + base + "_tables/*.h5")
    print paths
    #if len(sys.argv) > 2:
    #    filt = sys.argv[2]
    #    paths = [p for p in paths if filt in p]
    tasks = []
    for fn in paths:
        print "Processing ", fn
        try:
            f = tables.File(fn)
            groups = iter(iter(iter(f).next()).next()).next()
            groupnums = [i for i,g in enumerate(groups)]
            
            print "Processing ", fn, "; groups ", min(groupnums), "...", max(groupnums)
            tasks.extend([(fn, num) for num in groupnums])
            f.close()
        except Exception, e:
            print fn, e

    import multiprocessing   
    print "Starting %i processes to process %i tasks" % (multiprocessing.cpu_count()-1, len(tasks))
    p = multiprocessing.Pool(processes=multiprocessing.cpu_count()-1)
    p.map(worker, tasks)
