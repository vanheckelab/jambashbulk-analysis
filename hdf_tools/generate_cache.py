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

from pylab import *
import numpy as np
import pandas
import glob
import tables

import itertools
import pytables_tools

def loadnpy(fn):
    try:
        return np.load(fn)
    except ValueError:
        return np.load(fn + ".reduced.npy")

try:
    G_hess_data = loadnpy(basepath + "Cs.npy")
except Exception:
    print "Not including G_hess data in generated files!"
    G_hess_data = None

def dictlistgroupby(d, keys):
    key = lambda i: [i[key] for key in keys]
    
    d = sorted(d, key=key)
    for group, values in itertools.groupby(d, key):
        yield dict(zip(keys, group)), [v for v in values]


def determine_cs_indices(data, num=0):   
    basechange = array(unique(data["Nchanges"]))[num]
    
    gamma_plus = amin(data['gamma'][data['Nchanges'] > basechange])
    gamma_min = amax(data['gamma'][data['gamma'] < gamma_plus])
    
    return where(data['gamma'] == gamma_min)[0][0], where(data['gamma'] == gamma_plus)[0][0]

def determine_G(data):
    gamma = amin(data["gamma"][data["gamma"] > 0])
    sigma = data["s_xy"][data["gamma"] == gamma][0]
    
    G = sigma/gamma
    return G
    
    
def determine_Glambda2(data, i_min):
    d = data[data["gamma"] <= data["gamma"][i_min]]
    
    lambda2, G2, zero2 = polyfit(d["gamma"], d["s_xy"], 2)
    return {'lambda2': lambda2, 'G2': G2, 'zero2': zero2}
        
def determine_particle_movement(packing, i0, i1):
    particles0 = packing.SR.__getattr__("%04i" % i0).particles.read()
    particles1 = packing.SR.__getattr__("%04i" % i1).particles.read()
    
    x0, y0 = particles0["x"], particles0["y"]
    x1, y1 = particles1["x"], particles1["y"]
    
    return sqrt(sum((x0-x1)**2 + (y0-y1)**2))
    

def get_data_row(packing, packing_base):
    null,N,P,num = packing._v_pathname.split('/')
    num = int(num)
    
    data = packing.SR.data.read()

    G = determine_G(data)
    
    i_min, i_plus = determine_cs_indices(data)
    
    quadrG = determine_Glambda2(data, i_min)
    
    N = data["N"][0]
    P = data["P0"][0]    
    
    try:
        Cs = G_hess_data[(G_hess_data["N"] == N) * (abs(G_hess_data["P"] - P) < abs(P*1e-5)) * (G_hess_data[ "PackingNumber"] == num)][0]
    except Exception, e:
        print e
        Cs = None
    
    retval = {'N': N, 'P': P, 'num': num, 'G_first': G,
              'i_min': i_min, 'i_plus': i_plus}
    retval.update(quadrG)
    for i in range(1,7):
        key = "c%i" % i
        if Cs is None:
            retval[key] = np.nan
        else:
            retval[key] = Cs[key]

    retval.update(Neff_min= N-data["#rattler"][i_min],
                  Neff_plus= N-data["#rattler"][i_plus])
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
            v + "_min" : data[v][i_min],
            v + "_plus" : data[v][i_plus]
        })
        
    return retval


def get_data_rows():
    for fn in glob.glob(basepath + r'/N*_shear/*.h5'):
          sf = tables.File(fn)
          pf = tables.File(fn.replace("_shear", "_tables"))
          for N in sf.root:
              try:
                  Npf = pf.getNode(N._v_pathname)
                  for pressure in N:
                      try:
                          press = float(pressure._v_name[1] + "." + pressure._v_name[2:])
                          Ppf = Npf._v_children[Npf._v_children.keys()[[float(k[1:]) for k in Npf._v_children].index(press)]] 
                          for packing in pressure:
                              try:
                                  packpf = Ppf._v_children[packing._v_name]
                                  yield get_data_row(packing, packpf)
                              except Exception, e:
                                  print packing, e
                                  pass
                      except Exception, e:
                          print "-!!->", pressure, e
              except Exception, e:
                  print "-!!!!!!!!->", N, e
          sf.close()
          pf.close()
          
def store_data_rows():
    it = get_data_rows()
    first = it.next()
    cols = first.keys()
    
    f = tables.File(basepath + r"/shear_summary_noparticles.h5", "w")

    df = pandas.DataFrame([first.values()], columns=cols)
    df = df.to_records()
    df = df[cols]
    
    table = f.createTable(f.root, "data", df, 
                          expectedrows = 50000)

    for row in it:
        pytables_tools.add_to_table(table, row)
    
    f.flush(); f.close()
    
if __name__=="__main__":
    store_data_rows()
    pass
