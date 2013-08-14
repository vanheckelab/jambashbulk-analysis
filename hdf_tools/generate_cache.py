# -*- coding: utf-8 -*-
# Created on Wed Jul 04 11:47:14 2012
# @author: Merlijn van Deen <deen@physics.leidenuniv.nl>

"""

"""

import socket
if socket.gethostname() == "Isilwen":
    basepath = "d:/h5/"
else:
    basepath = r"/data/misc/granulargroup/h5/"

from getcc import get_first_cc
from pylab import *
import numpy as np
import pandas
import glob
import tables

import itertools
import pytables_tools

def get_data_row(packing, packing_base=None):
    null,N,P,num = packing._v_pathname.split('/')
    num = int(num)
    
    try:
        data = packing.SR.data.read()
    except tables.exceptions.NoSuchNodeError:
        return None

    try:
        before, after = get_first_cc(packing, data)   
    except Exception, e:
        if e.message in ["no contact change at all!", "step ordering", "no zero-change state", ]:
            print packing._v_pathname, e.message
            return None
        raise

    i_min = before['step#']
    i_plus = after['step#']
    
    N = data["N"][0]
    P = data["P0"][0]
    
    retval = dict(zip(['Nchanges_min', 'alpha_min', 'dH_base', 'dU_base', 'alpha_plus', 'delta_base', 'phi_min', 'sxy_calc_base', 'Uhelper_base', 'syy_calc_base', 'Nchanges_plus', 'L_base', 'num', 'i_min', 'P0_base', 's_xy_plus', 'L_min', 's_yy_plus', 'gamma_plus', 'N-_min', 'i_plus', 'phi_plus', 'maxGrad_base', 'P_calc_base', 'gg_min', 'syy_base', 'Z_base', 'N+_min', 'Neff_plus', 'P_min', 'Ncontacts_plus', 's_yy_min', 's_xy_min', 'dH_min', 'PackingNumber_base', 'Z_plus', 'gg_plus', 'phi_base', 'dU_min', 'Neff_min', 'P_base', 'Z_calc_base', 'alpha_base', 'N', 'P', 'U_calc_base', 'dU_plus', 'U_plus', 'P_plus', 'delta_plus', 'runtime (s)_base', 'L_plus', 'H_base', 'H_calc_base', 'N-_plus', 'gg_calc_base', 'H_plus', 'N - Ncorrected_base', 'delta_min', 's_xx_plus', 'N+_plus', 'sxx_calc_base', 'sxx_base', 'U_min', 'H_min', 'Z_min', 'gamma_min', 'Ncontacts_min', 'phi_calc_base', 'sxy_base', 'dH_plus', 's_xx_min'], itertools.repeat(nan)))    
    
    retval = {'N': N, 'P': P, 'num': num,
              'i_min': i_min, 'i_plus': i_plus}

    retval.update(Neff_min= N-before["#rattler"],
                  Neff_plus=N-after["#rattler"])

    if packing_base is not None:        
        for key in packing_base._v_attrs._v_attrnames:
            if isinstance(packing_base._v_attrs[key], np.float64):
                retval[key+"_base"] = packing_base._v_attrs[key]

    for v in ['L', 'P', 
              'gamma', 's_xy', 'Ncontacts',
              'Nchanges', 'N+', 'N-',
              'Z', 'alpha', 'delta',
              'phi', 's_xx', 's_yy',
              'U', 'dU',
              'H', 'dH', 
              'gg']:
                  # NOT L1 and L2 as those are vectors and the saving algorithm
                  # doesnt work for that...
        retval.update({
            v + "_min" : before[v],
            v + "_plus" : after[v],
        })
        
    return retval


def get_data_rows(path = basepath + r'/N*_shear/N*~P*.h5'):
    paths = glob.glob(path)
    print paths
    for fn in paths:
          print fn, 
          sff = tables.File(fn)
          try:
              pff = tables.File(fn.replace("_shear", "_tables"))
          except IOError, e:
              print e
              continue
          
          # ASSUME: only one N, P...
          try:
              sf = sff.root.__iter__().next().__iter__().next()
          except StopIteration:
              print "could not find shear data"
              continue
          
          try:
              pf = pff.root.__iter__().next().__iter__().next()
          except StopIteration:
              print "could not find packing data"
              continue
          
          print sf._v_pathname, pf._v_pathname
          
          for shearpack in sf:
              try:
                  basepack = pf.__getattr__(shearpack._v_name)
              except tables.exceptions.NoSuchNodeError:
                  basepack = None
              yield get_data_row(shearpack, basepack) 

          sff.close()
          pff.close()
          
def store_data_rows():
    it = get_data_rows()
    first = it.next()

    f = tables.File(basepath + r"/shear_summary_noparticles_20130730.h5", "w")

    cols = first.keys()
    df = pandas.DataFrame([first.values()], columns=cols)
    df = df.to_records()
    df = df[cols]
    
    table = f.createTable(f.root, "data", df, 
                          expectedrows = 50000)

    try:
        for row in it:
            if row is not None:
                pytables_tools.add_to_table(table, row)
        print "Done!"
    finally:
        print "Flushing!"
        f.flush(); f.close()
    
if __name__=="__main__":
    store_data_rows()
