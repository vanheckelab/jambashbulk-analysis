# -*- coding: utf-8 -*-
"""
Superseded by dodo_HPC.py.

Created on Thu May 16 16:32:05 2013

@author: Merlijn van Deen
"""
from pylab import *

import sys
sys.path.append('/home/merlijn/analysis/src/phdlib/')

import random
import warnings

from hdf_tools.getcc import get_first_cc, get_first_ccs
from packing_tools import V_harm

import tables
import glob,os



basepath = "/clusterdata/merlijn/h5/"

groups = [os.path.split(x)[1].split(".")[0].split("~") for x in glob.glob(basepath + "*_tables/*.h5")] #[("N16", "P1e-5")] #[("N22", "P1e-4"), ("N22", "P1e-6"), ("N64", "P1e-2"), ("N512", "P1e-2")]
groups = [g for g in groups if len(g) == 2]
tasks = []

import datetime
outfileprefix = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
os.makedirs(outfileprefix)

def gettasks():
    for (Ns, Ps) in groups:
        print "Parsing", (Ns, Ps),
        P = float(Ps[1] + "." + Ps[2:])
        if P > 2e-1 or P < 8e-7:
            print "skipping due to pressure"
            continue
        if int(Ns[1:]) > 1024:
            print "skipping large packings"
            continue
        sys.stdout.flush()
        postfix = Ns+Ps
        ftn = "%(Ns)s_tables/%(Ns)s~%(Ps)s.h5" % locals()
        fsn = "%(Ns)s_shear/%(Ns)s~%(Ps)s.h5" % locals()
        try:
            fs = tables.File(basepath + fsn)
        except IOError, e:
            print e
            continue
    
        ltasks = []
        for Nsub in fs.root:
            if Nsub._v_name[0] != "N":
                continue
            for Psub in Nsub:
                if Psub._v_name[0] != "P":
                    continue
                for grp in Psub:
                    ltasks.append((ftn, fsn, int(grp._v_name)))
        tasks.extend(ltasks)
        print len(ltasks), "added"
                    
        fs.close()
    return tasks
 
def runner((ftn, fsn, num)): 
    oldstderr = sys.stderr
    try:
        ft = tables.File(os.path.join(basepath, ftn))
        fs = tables.File(os.path.join(basepath, fsn))
        nums = "%04i" % num
        
        outfile = open(os.path.join(outfileprefix, str(os.getpid()) + ".out"), "a")
        errfile = open(os.path.join(outfileprefix, str(os.getpid()) + ".err"), "a")  

        sys.stderr = errfile

        def customwarn(message, category, filename, lineno, file=None, line=None):
            errfile.write(warnings.formatwarning(message, category, filename, lineno))
        
        warnings.showwarning = customwarn
        
        basepack = iter(iter(ft.root).next()).next().__getattr__(nums)

        Ns, Ps = os.path.split(ftn)[1].split(".")[0].split("~")
                
        # bepalen G van base packings
        packing = dict((x, basepack._v_attrs[x]) for x in basepack._v_attrs._v_attrnames)
        packing['particles'] = basepack.particles.read()
        
        contacts = V_harm.get_contacts(packing)
        K = V_harm.Hess(contacts, packing)
        el_mod = V_harm.El_Con(K, contacts['rattlers'], packing)
        
        G_base = el_mod["c"][2]
        
        # nu van packing vlak vóór en vlak ná contact change
        # voor moet gelijk zijn aan G_base, na zou anders zijn
        
        
        sgroup = iter(iter(fs.root).next()).next().__getattr__(nums)
        
        #before, after = get_first_cc(sgroup)
        bset, aset = get_first_ccs(sgroup)
        
        # First, determine the closest BEFORE
        for step in reversed(bset):
        #    print "before - trying %i" % step["step#"]
            
            node = fs.getNode(step["path"])
            if "particles" in node._v_children:
                break
    
        packing = dict((x, node._v_attrs[x]) for x in node._v_attrs._v_attrnames)
        packing['particles'] = node.particles.read()
        
        contacts = V_harm.get_contacts(packing)
        K = V_harm.Hess(contacts, packing)
        el_mod = V_harm.El_Con(K, contacts['rattlers'], packing)
        
        G_bef = el_mod["c"][2]
            
        # Next, determine the closest AFTER
        for step in aset:
        #    print "after - trying %i" % step["step#"]
            
            node = fs.getNode(step["path"])
            if "particles" in node._v_children:
                break
    
        packing = dict((x, node._v_attrs[x]) for x in node._v_attrs._v_attrnames)
        packing['particles'] = node.particles.read()
        
        contacts = V_harm.get_contacts(packing)
        K = V_harm.Hess(contacts, packing)
        el_mod = V_harm.El_Con(K, contacts['rattlers'], packing)
        
        G_aft = el_mod["c"][2]
        
        
        outfile.write("%s;%.2f;%i;%f;%f;%f\n" % (Ns[1:], around(log10(float(Ps[1:])),2), num, G_base, G_bef, G_aft))
    except Exception, e:
        import traceback
        traceback.print_exc(None, errfile)
    finally:
        sys.stderr = oldstderr
        errfile.flush()
        errfile.close()
        outfile.flush()
        outfile.close()
        
def runtasks(tasks):
    import progressbar
    import multiprocessing
    import time

    widgets = [progressbar.SimpleProgress(), ' (',  progressbar.Percentage(), ')  ', progressbar.Bar(), ' ', progressbar.ETA()]
    pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(tasks)).start()

    p = multiprocessing.Pool(processes=multiprocessing.cpu_count()-1)
    
    for i, _ in enumerate(p.imap(runner, tasks)):
        pbar.update(i)

    pbar.finish()

if __name__=="__main__" and "--run" in sys.argv:
    tasks = gettasks()
    random.seed(0)
    random.shuffle(tasks)
    print len(tasks), "tasks to process"
    if raw_input("Run all? ").lower() == "y":
        runtasks(tasks)
    runnum = int(raw_input("How many? "))
    if runnum:
        runtasks(tasks[:runnum])
