# coding: utf-8
from hpc import HessianPackingCalculator
import tables
from hdf_tools import getcc
import numpy as np
from numpy import mean


def RunOnH5File(source, target):
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
        f = tables.File(source).root.__iter__().next().__iter__().next()
        tabf = tables.File(source.replace("_shear", "_tables")).root.__iter__().next().__iter__().next()
    except StopIteration:
        return
    elements = []
    for i,item in enumerate(f):
        if item._v_name == "0000":
            continue

        try:
            d = item.SR.data.read()
        except tables.NoSuchNodeError:
            raise
            continue
            
        try:
            staticpack = getattr(tabf, item._v_name)
        except tables.NoSuchNodeError:
            raise
            continue
            
        try:
            # get static packing stuff
            zHPC = HessianPackingCalculator(staticpack)
            gmk, gbk = zHPC.find_first_ccs()
            delta = mean(zHPC.contacts["dij"][zHPC.contacts["dij"] > 0])
            zElCon = zHPC.get_el_con()
            del zHPC

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
                raise Exception("No after state found")

            aHPC = HessianPackingCalculator(f._v_file.getNode(aftpath))
            aElCon = aHPC.get_el_con()
            del aHPC
            
            # build return dict
            values = [delta, gmk, gbk] + zElCon.values() + aElCon.values()
            keys = ['mean_delta_base', 'gamma_lr_mk_base', 'gamma_lr_bk_base'] + \
                   [k + "_base" for k in zElCon.keys()] + \
                   [k + "_plus" for k in aElCon.keys()]
            
            dtypes = [(name, np.float64) for name in keys]
            
            values = [d["N"][0], d["P"][0], int(item._v_name)] + values
            dtypes = [("N", np.int32), ("P", np.float64), ("num", np.int32)] + dtypes
            rd = np.array([tuple(values)], dtype=dtypes)
                         
            elements.append(rd)

        except Exception as e:
            raise
            errf = open('errors', 'a')
            errf.write("%s raised %r" % (item._v_pathname, e))
            errf.close()
            continue
 
    if elements:
        np.save(target, np.concatenate(elements))
        
