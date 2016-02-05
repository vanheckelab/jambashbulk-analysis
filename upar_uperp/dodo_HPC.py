# coding: utf-8
import os
from hpc import HessianPackingCalculator
import tables
from hdf_tools import getcc
import numpy as np
from numpy import mean
import traceback


def RunOnH5File(source, target, targetdir):
    """ For all packings in an H5 shear file, calculate:
    
        On the BASE packing:
           * γmk, γbk in linear response
           * <δ>.
           * c1...c6
           * Galpha (simple shear direction)

        On the AFTER CC packing:
           * c1...c6
           * Galpha (simple shear direction)

        NB. There should also be a _tables version in the same directory...
        
        source = "/mnt/user/valhallasw/auto/h5/N32~P2154e-6_shear.h5"
        target = "/mnt/user/valhallasw/auto/linres/N32~P2154e-6_linres.npy"
    """
    try:
        os.makedirs(targetdir)
    except OSError:
        pass
        
    try:
        f = tables.File(source).root.__iter__().next().__iter__().next()
        tabf = tables.File(source.replace("_shear", "_tables")).root.__iter__().next().__iter__().next()
    except StopIteration:
        return
    elements = []
    for i,item in enumerate(f):
        if item._v_name in ["0000", "0001"]:
            continue

        try:
            d = item.SR.data.read()
            staticpack = getattr(tabf, item._v_name)
           
            # get static packing stuff
            zHPC = HessianPackingCalculator(staticpack)
            ccGammas = zHPC.find_first_ccs()
            delta = mean(zHPC.contacts["dij"][zHPC.contacts["dij"] > 0])
            zElCon = zHPC.get_el_con()
            zElCon.update(ccGammas)
            
            try:
                # now we also want to store the linear response particle data!
                u_par, u_perp = zHPC.get_uparrs(2, True)
                deltas = zHPC.contacts["dijfull"]
                Rijs = zHPC.contacts["rij"]
                
                np.save(targetdir + "/%s.npy" % item._v_name, {'u_par': u_par, 'u_perp': u_perp, 'deltas': deltas, 'Rijs': Rijs})
            except Exception as e:
                errf = open('errors', 'a')
                errf.write("\n\n{source} #{num}\n===============================\n".format(source=source, num=item._v_name))
                traceback.print_exc(file=errf)
                errf.write("\n")
                errf.write(repr(e))
                errf.write("\n")
                errf.close()
            del zHPC
            
            # before shear packing attrs
            bef, aft = getcc.get_first_cc(item, d)

            ag = bef["gamma"]

            sortg = np.argsort(d["gamma"])

            for gi in sortg[d[sortg]["gamma"] <= ag]:
                path = d[gi]["path"]
                if d[gi]["gamma"] < ag/2:
                    continue
                SR = item._v_file.getNode(path)
                if hasattr(SR, 'particles'):
                    aftpath = path
                    break
            else:
                aftpath = None
                print item._v_name, "before not found"

            if aftpath:
                HPC = HessianPackingCalculator(f._v_file.getNode(aftpath))
                bElCon = HPC.get_el_con()
                del HPC
            else:
                bElCon = {'c%i'%i: np.nan for i in [1,2,3,4,5,6]}

            # after shear packing attrs
            bef, aft = getcc.get_first_cc(item, d)

            ag = aft["gamma"]

            sortg = np.argsort(d["gamma"])

            for gi in sortg[d[sortg]["gamma"] >= ag]:
                path = d[gi]["path"]
                if d[gi]["gamma"] > 2*ag:
                    continue
                SR = item._v_file.getNode(path)
                if hasattr(SR, 'particles'):
                    aftpath = path
                    break
            else:
                aftpath = None
                print item._v_name, "before not found"

            if aftpath:
                HPC = HessianPackingCalculator(f._v_file.getNode(aftpath))
                aElCon = HPC.get_el_con()
                del HPC
            else:
                aElCon = {'c%i'%i: np.nan for i in [1,2,3,4,5,6]}
          
            # build return dict
            values = [delta] + zElCon.values() + bElCon.values() + aElCon.values()
            keys = ['mean_delta_base'] + \
                   [k + "_base" for k in zElCon.keys()] + \
                   [k + "_min" for k in bElCon.keys()] + \
                  [k + "_plus" for k in aElCon.keys()]
            
            dtypes = [(name, np.float64) for name in keys]
            
            values = [d["N"][0], d["P"][0], int(item._v_name)] + values
            dtypes = [("N", np.int32), ("P", np.float64), ("num", np.int32)] + dtypes
            rd = np.array([tuple(values)], dtype=dtypes)
                         
            elements.append(rd)
            


        except Exception as e:
            errf = open('errors', 'a')
            errf.write("\n\n{source} #{num}\n===============================\n".format(source=source, num=item._v_name))
            traceback.print_exc(file=errf)
            errf.write("\n")
            errf.write(repr(e))
            errf.write("\n")
            errf.close()
            continue
 
    if elements:
        np.save(target, np.concatenate(elements))
        
