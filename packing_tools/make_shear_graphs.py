# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 11:42:13 2012

@author: Merlijn van Deen
"""
import sys
sys.path.append("../hdf_tools")

import pylab
import pandas
import load_log, load_packing
from load_packing import loadPackingData
from itertools import izip
from pytables_import_shear import read_csv
from generate_cache import determine_cs_indices
from V_harm import get_contacts

def read_log(name):
    data1, comment = read_csv(base_path + "data" + name + ".txt")
    data2, comment = read_csv(base_path + "log" + name + ".txt")
    
    data = data1.join(data2, rsuffix="_")
    data["gamma"] = data["gamma_alpha"]
    del data["gamma_alpha"]
    return data, comment

def loadPackings(filename):
    try:
        f = open(filename)
    except Exception:
        import time
        time.sleep(30)
        f = open(filename)

    data = ""
    for line in f:
        data += line
        if line.strip() == "" and lastline.strip() == "":
            yield data.strip()
            data = ""
        lastline = line

def determine_all_cs_indices(data):
    num = 0
    try:
        while True:
            yield determine_cs_indices(data, num)
            num += 1
    except ValueError:
        raise StopIteration

base_path = r"D:\h5\N2048-large-shear" + "\\"
#name = "N2048~P1e-6~SR050~step101~9004"
name = "N2048~P1e-2~SS0e-2~step2000~9004"
data, comment = read_log(name)

#indices = [0] + [x[1] for x in determine_all_cs_indices(data)]

def scalefig(pack, extent=0.1):
        xmin = (                min(0, pack['L2'][0]))
        xmax = (pack['L1'][0] + max(0, pack['L2'][0]))
        dx = abs(xmax-xmin)
        xmin = xmin-extent*dx        
        xmax = xmax+extent*dx 
        
        ymin = 0
        ymax = pack["L2"][1]
        dy = abs(ymax-ymin)
        ymin = ymin-extent*dy
        ymax = ymax+extent*dy
        
        axis((xmin, xmax, ymin, ymax))

def plotunitcell(pack, **kwargs):
    if 'color' not in kwargs:
        kwargs['color'] = 'gray'

    steps = [array([0,0])]
    for d in [pack["L1"], pack["L2"], -pack["L1"], -pack["L2"]]:
        steps.append(steps[-1] + d)
    
    plot(*zip(*steps), **kwargs)

def plotparticle(pack, i, **kwargs):
    if 'fc' not in kwargs:
        kwargs['fc'] = "none"

    part = pack['particles'][i]
    
    circ=pylab.Circle((part['x'], part['y']), radius=part['r'], **kwargs)
    ax=gca()
    ax.add_patch(circ)

startnumconts = None
def plotparticles(pack, **kwargs):
    global startnumconts

    conts = get_contacts(pack)
    numconts = sum(conts['connmatrix'], axis=0)    
    
    if startnumconts is None:
        startnumconts = numconts
        
    for i in range(len(pack['particles'])):
        if True:
            if numconts[i] < startnumconts[i]:
                kwargs['fc'] = 'red'
            elif numconts[i] > startnumconts[i]:
                kwargs['fc'] = 'blue'
            else:
                kwargs['fc'] = 'none'
        else:
            rg = 1e-3
            dx = (pack['particles'][i]['x'] - startpack['particles'][i]['x']) - \
                 (((pack["L2"][0] - startpack["L2"][0])/startpack["L2"][1]) * startpack['particles'][i]["y"])
            
            if dx <= -rg:
                dx = -rg + 1e-30
            elif dx >= rg:
                dx = rg - 1e-30
            
            kwargs['fc'] = cm.spectral((dx+rg)/(2*rg))
        plotparticle(pack, i, **kwargs)

def doPackingStuff(pack):
    plotunitcell(pack)
    plotparticles(pack)
#
#ioff()
#startpack = None
#startnumconts = None
#for (pack, (i, subdata)) in izip(loadPackings(base_path + "particles" + name + ".txt"), data.iterrows()):
#    if True: #i in indices:
##        print i, indices.index(i), subdata["gamma"]
#        print i, subdata["gamma"]
#
#        pack = loadPackingData(pack)
#
#        if startpack is None:
#            startpack = pack
#
#        figure()
#        subplot(111, aspect='equal')
#        scalefig(startpack)
#        doPackingStuff(pack)
#        #savefig("dx/" + name + ('~CC%03i.png' % indices.index(i)))
#        savefig(r"D:\h5\N2048-large-shear\dx/" + name + ('~%04i.png' % i))
#        close()