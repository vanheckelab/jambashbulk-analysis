# -*- coding: utf-8 -*-
"""
Functions to plot packings (i.e. plot circles on a matplotlib graph).

Created on Wed Dec 12 11:42:13 2012

@author: Merlijn van Deen
"""
import sys, os
sys.path.append(os.path.split(__file__)[0])
sys.path.append(os.path.join(os.path.split(__file__)[0], "..", "hdf_tools"))

import pylab
from pylab import *
import pandas
import load_log, load_packing
from load_packing import loadPackingData
from itertools import izip
from pytables_import_shear import read_csv
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

def scalefig(pack, extent=0.1):
    """ Scales axis to roughly correspond to the plotted packing.
        extent=0.1 specifies extra padding (0.1 = 10%)
    """
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

def plotunitcell(pack, offset=(0,0), **kwargs):
    """ Draw the unit cell parallelogram. 
        offset specifies the (x,y) position of the lower left corner.
        kwargs are passed to plot().
    """
    if 'color' not in kwargs:
        kwargs['color'] = 'gray'

    steps = [array([0,0])]
    for d in [pack["L1"], pack["L2"], -pack["L1"], -pack["L2"]]:
        steps.append(steps[-1] + d)
    steps = array(steps)
    steps += offset
    
    plot(*zip(*steps), **kwargs)

def plotparticle(pack, i, offset=(0,0), **kwargs):
    """ Plot a single particle pack[i].
        offset (x,y) (data coordinates) is applied to the particle position
        kwargs are passed to pylab.Circle
    """
    if 'fc' not in kwargs:
        kwargs['fc'] = "none"

    part = pack['particles'][i]
    
    circ=pylab.Circle((part['x'] + offset[0], part['y'] + offset[1]), radius=part['r'], **kwargs)
    ax=gca()
    ax.add_patch(circ)

startnumconts = None
def plotparticles(pack, offset=(0,0), coloring=None, **kwargs):
    """ Plot all particles in a packing.
        offset (x,y) (data coordinates) is applied to all particles.
        coloring can be None (default), 'contactchange' or 'movement'.
        the latter two compare the current packing to the initial packing
        and is useful to plot a set of figures at different strains, comparing
        each to the initial state. Not used much, so likely to be buggy/broken.

    """
    global startnumconts

    conts = get_contacts(pack)
    numconts = sum(conts['connmatrix'], axis=0)    
    
    if startnumconts is None:
        startnumconts = numconts
        
    for i in range(len(pack['particles'])):
        if coloring == "contactchange":
            if numconts[i] < startnumconts[i]:
                kwargs['fc'] = 'red'
            elif numconts[i] > startnumconts[i]:
                kwargs['fc'] = 'blue'
            else:
                kwargs['fc'] = 'none'
        elif coloring == "movement":
            rg = 1e-3
            dx = (pack['particles'][i]['x'] - startpack['particles'][i]['x']) - \
                 (((pack["L2"][0] - startpack["L2"][0])/startpack["L2"][1]) * startpack['particles'][i]["y"])
            
            if dx <= -rg:
                dx = -rg + 1e-30
            elif dx >= rg:
                dx = rg - 1e-30
            
            kwargs['fc'] = cm.spectral((dx+rg)/(2*rg))
        plotparticle(pack, i, offset=offset, **kwargs)

def plotcontacts(pack, offset=(0,0), lwfactor=200, **kwargs):
    """ Plot all contacts in a packing, with a width corresponding to the overlap/force.
        offset (x,y; data coordinates) is applied to all contacts,
        lwfactor is used to determine line thickness: linewidth (pt) = overlap delta * lwfactor
        kwargs are passed to plot()
    """
    if 'color' not in kwargs:
        kwargs['color'] = 'black'
    conts = get_contacts(pack)
    Nparts = conts['dij'].shape[0]
    xij, yij = conts['xij'], conts['yij']
    
    for row in range(Nparts):
        for col in range(row):
            dij = conts['dij'][row, col] 
            if dij > 0:
                part = pack['particles'][col]
                x,y = part['x'], part['y']
                dx, dy = xij[row, col], yij[row, col]
                kwargs.setdefault('lw', dij*lwfactor)
                plot([x+offset[0], x+dx+offset[0]], [y+offset[1], y+dy+offset[1]], **kwargs) # *6000 works for N256 P1e-4 
                
                

def doPackingStuff(pack, newfig=True, plot_contacts=True, plot_unitcell=True):
    """ Plot a packing with periodic copies\
           - plot contacts if plot_contacts=True,
           - plot in a new figure and rescale if newfig=True,
           - plot unit cell if plot_unitcell=True
    """
    if newfig:
        figure()
        subplot(111, aspect= 'equal')
        
    if plot_unitcell:
        plotunitcell(pack, lw=4, color="blue")
        
    plotparticles(pack, lw=3, coloring="contactchange")
    
    if plot_contacts:
        plotcontacts(pack)
    
    for offset_num in [(-1,0), (0,-1), (-1,1), (1,-1), (-1,-1), (1,0), (0,1), (1,1)]:
        realoffset = (pack["L1"] * offset_num[0] + pack["L2"] * offset_num[1])
        plotparticles(pack, realoffset, ls="dashed", lw=2, coloring="contactchange") # color="red",
        if plot_contacts:
            plotcontacts(pack, realoffset, color="gray")
    
    if newfig:
        scalefig(pack, extent=0.15)