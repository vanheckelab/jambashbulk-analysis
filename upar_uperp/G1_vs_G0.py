# -*- coding: utf-8 -*-
"""
Created on Thu May 16 16:32:05 2013

@author: Merlijn van Deen
"""

from hdf_tools.getcc import get_first_cc, get_first_ccs
from packing_tools import V_harm

import tables
close('all')
import glob,os

basepath = "d:/h5/"

groups = [os.path.split(x)[1].split(".")[0].split("~") for x in glob.glob(basepath + "*_tables/*.h5")] #[("N16", "P1e-5")] #[("N22", "P1e-4"), ("N22", "P1e-6"), ("N64", "P1e-2"), ("N512", "P1e-2")]

groups = groups[:2]
tasks = []

def gettasks():
    for (Ns, Ps) in groups:
        print "Parsing", (Ns, Ps)
        postfix = Ns+Ps
        ftn = "%(Ns)s_tables\%(Ns)s~%(Ps)s.h5" % locals()
        fsn = "%(Ns)s_shear\%(Ns)s~%(Ps)s.h5" % locals()
        fs = tables.File(basepath + fsn)
    
        for Nsub in fs.root:
            if Nsub._v_name[0] != "N":
                continue
            for Psub in Nsub:
                if Psub._v_name[0] != "P":
                    continue
                for grp in Psub:
                    tasks.append((ftn, fsn, int(grp._v_name)))
                    
        fs.close()
    return tasks
 
def runner((ftn, fsn, num)): 
    print os.getpid()
    try:
        ft = tables.File(os.path.join(basepath, ftn))
        fs = tables.File(os.path.join(basepath, fsn))
        nums = "%04i" % num
        
        outfile = open(str(os.getpid()) + ".out", "a")
        errfile = open(str(os.getpid()) + ".err", "a")  
        
        basepack = iter(iter(ft.root).next()).next().__getattr__(nums)
                
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
            print "before - trying %i" % step["step#"]
            
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
            print "after - trying %i" % step["step#"]
            
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
        raise
    finally:
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

          
#            node = fs.getNode(after["path"])
#            packing = dict((x, node._v_attrs[x]) for x in node._v_attrs._v_attrnames)
#            packing['particles'] = node.particles.read()
#            
#            contacts = V_harm.get_contacts(packing)
#            K = V_harm.Hess(contacts, packing)
#            el_mod = V_harm.El_Con(K, contacts['rattlers'], packing)
#            
#            G_aft = el_mod["c"][2]
#            
#            figure("G_base_G_bef"+postfix)
#            plot(G_base, G_bef, ".")
#            
#            figure("G_base_G_aft"+postfix)
#            plot(G_base, G_aft, ".")
#            print nums, G_base, G_bef, G_aft
#        except tables.exceptions.NoSuchNodeError, e:
#            print postfix, n, e
#        except V_harm.PackingException, e:
#            print postfix, n, e
            
            
#    figure("G_base_G_bef"+postfix)
#    title("N=%(Ns)s, P=%(Ps)s" % locals())
#    xlabel("$G_0$")
#    ylabel("$G_-$")
#    
#    figure("G_base_G_aft"+postfix)
#    title("N=%(Ns)s, P=%(Ps)s" % locals())
#    xlabel("$G_0$")
#    ylabel("$G_1$")
#    
#    xmin, xmax, ymin, ymax = axis()
#    axis([0, xmax, ymin, ymax])
#    
#    axhline(0, color="k")
#    
#    mx = min(xmax, ymax)
#    plot([0,mx],[0,mx], "k-")
#    
#    l2 = array([0,0])
#    trans_angle = gca().transData.transform_angles(array((45,)),
#                                                   l2.reshape((1,2)))[0]
#    text(mx/2, mx/2, "shear thickening\n\nshear thinning", rotation=trans_angle,horizontalalignment='center',
#         verticalalignment='center' )
#         
#    text(mx/2, 0, "\nunstable", horizontalalignment="center", verticalalignment="top")