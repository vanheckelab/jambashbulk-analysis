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
from pylab import *
from packing_tools import V_harm
from packing_tools.make_shear_graphs import plotparticles, plotunitcell
reload(V_harm)
import tables
import itertools

def subset(dictionary, keys):
    return {k: v for (k,v) in dictionary.iteritems() if k in keys}

from hpc import HessianPackingCalculator

if __name__=="__main__":
    import glob
    from PyQt4.QtGui import QApplication
    
    close('all')
    import time
    t0 = time.time()
    counter = 0
    total = 6000
    
    for fn in glob.glob("d:/h5/N1024_tables/*~P1e-1*.h5"):
        f = tables.File(fn)
        groups = iter(iter(iter(f).next()).next()).next()
        
        u_parrs = []
        u_perps = []
        
        npgfn = "d:/h5/" + os.path.split(fn)[1].split(".")[0] + "~u_parr_perps.npz"
        print "------> ", npgfn
        try:
            data = np.load(npgfn)
            locals().update(data)
            data.close()
            save = False
        except IOError:
            for i,group in enumerate(groups):
                QApplication.processEvents()
                
                npfn = "d:/h5/" + os.path.split(fn)[1].split(".")[0] + "~" + group._v_name + "~part_u_s.npz"
                print group._v_pathname, "; ",npfn, "; ",
                if group._v_pathname in ["/N512/p1.0000e-06/9050"]:
                    print "<SKIP DUE TO ERROR>"
                    continue
                try:
                    data = np.load(npfn)
                    locals().update(data)
                    data.close()
                    save = False
                    raise IOError
                except IOError:
        
                    HPC = HessianPackingCalculator(group)
                    QApplication.processEvents()
                    u_parr, u_parr_scaled, u_perp, u_perp_scaled = HPC.get_scaled_uparperps()
                    save = True
                
                u_parrs.append(u_parr_scaled)
                u_perps.append(u_perp_scaled)
                
                t1 = time.time()
                counter += 1        
                tpp = (t1-t0) / counter
                timeleft = datetime.timedelta(seconds = (total-counter) *  tpp)
                eta = (datetime.datetime.now() + timeleft).strftime("%d %b %Y %H:%M:%S")
    
                if save:
                    np.savez(npfn,
                             u_parr = u_parr, u_parr_scaled = u_parr_scaled,
                             u_perp = u_perp, u_perp_scaled = u_perp_scaled,
                             K0_inv = HPC.K0_inv,
                             rij = HPC.contacts['rij'], dij = HPC.contacts['dij'])
                print "%i done; %f s per packing; %s left (eta: %s)" % (counter, tpp, timeleft, eta) 
                #np.save(npfn + ".K0.np", HPC.K0)
                raise Exception
                
            u_parrs = reduce(np.append, u_parrs)
            u_perps = reduce(np.append, u_perps)
            np.savez(npgfn, u_parrs=u_parrs, u_perps=u_perps)
        
        
        #from cdf import densplot
        figure("u_parrs"); hist(u_parrs, bins=linspace(-10,10,500), normed=True, histtype="step", label=fn)
        #figure(fn + " u_parrs"); hist(u_parrs, bins=linspace(-10,10,500), normed=True, histtype="step")
        figure("u_perps"); hist(u_perps, bins=linspace(0,6,500), normed=True, histtype="step", label=fn)
        #figure(fn + " u_perps"); hist(u_perps, bins=linspace(0,6,500), normed=True, histtype="step")