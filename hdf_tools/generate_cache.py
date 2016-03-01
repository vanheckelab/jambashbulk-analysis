# -*- coding: utf-8 -*-
# Created on Wed Jul 04 11:47:14 2012
# @author: Merlijn van Deen <deen@physics.leidenuniv.nl>

"""

"""

import socket
import os
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
from scipy.optimize import curve_fit

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

    data_reduced = data[:i_min+1]
    data_reduced = data_reduced[data_reduced['gamma'] <= before['gamma']]

    if len(data_reduced) < 5:
        Glin = Glinerr = Gquad = Gquaderr = l = lerr = np.nan
    else:
        # linear fit through (0,0)
        func = lambda gamma, G: gamma*G
        (Glin,), pcov = curve_fit(func, data_reduced['gamma'], data_reduced['s_xy'])
        try:
            Glinerr = sqrt(pcov[0][0])
        except TypeError:
            Glinerr = np.nan

        # quadratic fit through (0,0)
        func = lambda gamma, l, G: gamma*G + gamma*gamma*l
        (l, Gquad), pcov = curve_fit(func, data_reduced['gamma'], data_reduced['s_xy'])
        try:
            lerr = sqrt(pcov[0][0])
            Gquaderr = sqrt(pcov[1][1])
        except TypeError:
            lerr = Gquaderr = np.nan

    retval = dict(zip(['Nchanges_min', 'alpha_min', 'dH_base', 'dU_base', 'alpha_plus', 'delta_base', 'phi_min', 'sxy_calc_base', 'Uhelper_base', 'syy_calc_base', 'Nchanges_plus', 'L_base', 'num', 'i_min', 'P0_base', 's_xy_plus', 'L_min', 's_yy_plus', 'gamma_plus', 'N-_min', 'i_plus', 'phi_plus', 'maxGrad_base', 'P_calc_base', 'gg_min', 'syy_base', 'Z_base', 'N+_min', 'Neff_plus', 'P_min', 'Ncontacts_plus', 's_yy_min', 's_xy_min', 'dH_min', 'PackingNumber_base', 'Z_plus', 'gg_plus', 'phi_base', 'dU_min', 'Neff_min', 'P_base', 'Z_calc_base', 'alpha_base', 'N', 'P', 'U_calc_base', 'dU_plus', 'U_plus', 'P_plus', 'delta_plus', 'runtime (s)_base', 'L_plus', 'H_base', 'H_calc_base', 'N-_plus', 'gg_calc_base', 'H_plus', 'N - Ncorrected_base', 'delta_min', 's_xx_plus', 'N+_plus', 'sxx_calc_base', 'sxx_base', 'U_min', 'H_min', 'Z_min', 'gamma_min', 'Ncontacts_min', 'phi_calc_base', 'sxy_base', 'dH_plus', 's_xx_min', 'Glin', 'Glin_err', 'Gquad', 'Gquad_err', 'lambdaquad', 'lambdaquad_err'], itertools.repeat(nan)))    
    
    retval = {'N': N, 'P': P, 'num': num,
              'i_min': i_min, 'i_plus': i_plus,
              'Glin': Glin, 'Glinerr': Glinerr,
              'Gquad': Gquad, 'Gquaderr': Gquaderr,
              'lambdaquad': l, 'lambdaquaderr': lerr}

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

STEPORDERS = []
ignore_errnos = [2,22] #Invalid Argument due to empty file
def get_data_rows(path, linrespath):
    paths = sorted(glob.glob(path))
    print paths
    for fn in paths:
        print fn, 
        sff = tables.File(fn)
        try:
            pff = tables.File(fn.replace("_shear", "_tables"))
        except IOError, e:
            print e
            continue
        try:
            linresnpy = np.load(os.path.join(linrespath,
                                             os.path.split(fn)[-1].replace("_shear.h5", "_linres.npy"))
                               )
        except IOError, e:
            if e.errno in ignore_errnos:
                continue
            raise
        
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
        
        print sf._v_pathname, pf._v_pathname,
        
        i = 0
        for shearpack in sf:
            try:
                if shearpack._v_name in ['0000', '0001']:
                    continue # SKIP packing 1 - because of borked simulations continueing with #1
                basepack = pf.__getattr__(shearpack._v_name)
            except tables.exceptions.NoSuchNodeError:
                basepack = None
            try:
                row = get_data_row(shearpack, basepack) 
            except Exception, e:
                if e.args[0][0] == 'step ordering':
                    print e
                    STEPORDERS.append(shearpack._v_pathname)
                    row = None
                else:
                    raise
                    
            if row:
                #update row with linres data
                lrd = linresnpy[linresnpy["num"] == int(shearpack._v_name)]
                for name in lrd.dtype.names:
                    if name in ["N", "P", "num"]:
                        continue
                    row[name] = lrd[name][0] if lrd else nan
                    
            yield row
            i += 1
        print i
        sff.close()
        pff.close()

def store_data_rows(frompath="/mnt/user/valhallasw/auto/h5/*_shear.h5", linrespath="/mnt/user/valhallasw/auto/linres/", topath="/mnt/user/valhallasw/auto/shear_summary_noparticles_20150421.h5"):
    it = get_data_rows(frompath, linrespath)
    first = it.next()

    f = tables.File(topath, "w")

    def sortkey(colname):
        prefix = 100
        if "_base" in colname:
            prefix = 0
        elif "_min" in colname:
            prefix = 1
        elif "_plus" in colname:
            prefix = 2
            
        return (prefix, colname.lower())
    
    cols = sorted(first.keys(), key=sortkey)
    
    df = pandas.DataFrame([[first[x] for x in cols]], columns=cols)
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
    store_data_rows(*sys.argv[1:])
