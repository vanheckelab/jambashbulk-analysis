# -*- coding: utf-8 -*-
"""
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
    if 'color' not in kwargs:
        kwargs['color'] = 'gray'

    steps = [array([0,0])]
    for d in [pack["L1"], pack["L2"], -pack["L1"], -pack["L2"]]:
        steps.append(steps[-1] + d)
    steps = array(steps)
    steps += offset
    
    plot(*zip(*steps), **kwargs)

def plotparticle(pack, i, offset=(0,0), **kwargs):
    if 'fc' not in kwargs:
        kwargs['fc'] = "none"

    part = pack['particles'][i]
    
    circ=pylab.Circle((part['x'] + offset[0], part['y'] + offset[1]), radius=part['r'], **kwargs)
    ax=gca()
    ax.add_patch(circ)

startnumconts = None
def plotparticles(pack, offset=(0,0), coloring=None, **kwargs):
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

def plotcontacts(pack, offset=(0,0), **kwargs):
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
                plot([x+offset[0], x+dx+offset[0]], [y+offset[1], y+dy+offset[1]], lw = dij*200, **kwargs) # *6000 works for N256 P1e-4 
                
                

def doPackingStuff(pack, newfig=True, plot_contacts=True, plot_unitcell=True):
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

if __name__=="__main__":
    CClast = CC = 0
    
    if "p" not in locals():
        #p = load_packing.loadPackings(r'C:\Users\Merlijn van Deen\Documents\Virtualbox Shared Folder\Packings\Packings\N32~P1e-2\particlesN32~P1e-2~SS2e-2~step0100~0017.txt')
        #p = load_packing.loadPackings(r'C:\Users\Merlijn van Deen\Documents\Virtualbox Shared Folder\Packings\Packings\N256~P1e-4\particlesN256~P1e-4~SS1e-3~step0250~9481.txt'); outputdir = "N256~P1e-4\\maxi"
        #p = load_packing.loadPackings(r'C:\Users\Merlijn van Deen\Documents\Virtualbox Shared Folder\Packings\Packings\N256~P1e-4\particlesN256~P1e-4~SS4e-4~step0100~9481.txt')
        p = load_packing.loadPackings(r'C:\Users\Merlijn van Deen\Documents\Virtualbox Shared Folder\Packings\Packings\N256~P1e-4\particlesN256~P1e-4~SS1e-4~step0150~9481.txt'); outputdir = "N256~P1e-4\\miniN256"
        
        # alleen voor shear position files
        p = [p[0]] + p[1::2]

        #p = load_packing.loadPackings(r"C:\Users\Merlijn van Deen\Documents\Virtualbox Shared Folder\Packings\N256~P1e-4~9481~creation.positions")[::2]        
        #p = load_packing.loadPackings(r"C:\Users\Merlijn van Deen\Documents\Virtualbox Shared Folder\Packings\Packings\N16~P1e-3/N16~P1e-3~9516.intermediate.positions")

 
    ioff()   
    def get_alpha(pack):
        return pack["L2"][0] / pack["L2"][1]
    
    def get_volume(pack):
        return pack["L1"][0] * pack["L2"][1]
        
    def get_stress(pack):
        conts = get_contacts(pack)
        volume = get_volume(pack)
        stresssum = sum(nan_to_num((conts['dij'] / conts['rij']) * conts['xij'] * conts['yij']))
        return - stresssum / volume / 2 # double counted bonds
        
    alpha_zero = get_alpha(p[0])
    gammas = []
    sigmas = []
    rearrplots = []
    
#    linslope = 6.7185519883010259e-05/0.0015666020948904083 # MANUAL
    linslope = 3.298000702835313e-07 / 3.4795112498202385e-05
    maxstrain = 4.0000000000000000027e-04 * 0.2
    maxstress = 1e-6
    minstress = 0
    
    for i,pack in enumerate(p):
        animate()
        alpha = get_alpha(pack)
        conts = get_contacts(pack)
        print i, 
        
        gamma = alpha - alpha_zero
        sigma = get_stress(pack)
        print "gamma: ", gamma,
        print "sigma: ", sigma,
        if startnumconts is not None:
            CC = sum(abs(startnumconts - sum(conts['connmatrix'], axis=0)))
            print "CC: ", CC
        else:
            CC = 0
        
        gammas.append(gamma)
        sigmas.append(sigma)

        if True:
            figure(figsize=(12,6))
            subplot(122)
            
            #plot(linspace(0,maxstrain), linspace(0,maxstrain)*linslope, "r--", lw=3)
            plot(gammas, sigmas, color="black", lw=3)
            plot(gamma, sigma, "o", color="red", markersize=15)
            
            if CC != CClast:
                rearrplots.append([gamma, sigma])
            
            import itertools;
            colors = itertools.cycle(["b", "r"]) # manual hack
            for g,s in rearrplots:
                plot([g, g, 0], [minstress, s, s], color=colors.next(), lw=1)
                
            CClast = CC
            axhline(0)
            axis(xmin=0, xmax=maxstrain, ymin=minstress, ymax=maxstress)
            
            xticks([])
            yticks([])
           
            subplot(121, aspect= 'equal')
            
        if False:
            figure()
            subplot(111, aspect='equal')
        
        if True:
            baseoffset = (0.5 * pack["L1"] + 0.5 * pack["L2"])
            plotunitcell(pack, lw=4, color="blue")
                        
            for offset_num in [(0,0), (-1,0), (0,-1),(-1,-1)]:
                realoffset = (pack["L1"] * offset_num[0] + pack["L2"] * offset_num[1]) + baseoffset
                plotparticles(pack, realoffset, lw=3, coloring="contactchange")
                plotcontacts(pack, realoffset)
                            
            #xlabel("step: % 3i" % i, fontname="Consolas", fontsize=30) 
            xticks([])
            yticks([])
            scalefig(p[0], extent=0.15)

            savefig(r"C:\Users\Merlijn van Deen\Documents\Virtualbox Shared Folder\Packings\Packings" + "\\" + outputdir + "\\%05i.png" % i, bbox_inches='tight', pad_inches=0.1)
         
        close()
      
#    for i,pack in enumerate(load_packing.loadPackings(r'C:\Users\Merlijn van Deen\Documents\Virtualbox Shared Folder\Packings\Packings\N256~P1e-4')):
#        print i
#        doPackingStuff(pack)
#        xlabel("step: % 3i" % i, fontname="Consolas", fontsize=30) 
#        xticks([])
#        yticks([])
#        scalefig(p[0], extent=0.15)
#        
#        savefig("d:/h5/rearstep_sim/%05i.png" % i, bbox_inches='tight', pad_inches=0.1)
#        close()

    # do the first one manually. Again.
#    pack = p[0]
#    figure()
#    subplot(111, aspect= 'equal')
#    plotunitcell(pack, lw=4, color="blue")
#    plotparticles(pack, lw=3)
#    
#    for offset_num in [(-1,0), (0,-1), (-1,1), (1,-1), (-1,-1), (1,0), (0,1), (1,1)]:
#        realoffset = (pack["L1"] * offset_num[0] + pack["L2"] * offset_num[1])
#        plotparticles(pack, realoffset, color="red", ls="dashed", lw=3)
#    
#    scalefig(pack, extent=0.15)
#    xlabel("step: % 3i" % 0, fontname="Consolas", fontsize=30) 
#    xticks([])
#    yticks([])
#    scalefig(p[0], extent=0.15)
#    
#    savefig("d:/h5/initial_sim/%05i_manual.png" % 0, bbox_inches='tight', pad_inches=0.1)
#    close()
#
#    for i,pack in enumerate(p):
#        print i
#        doPackingStuff(pack)
#        xlabel("step: % 3i" % i, fontname="Consolas", fontsize=30) 
#        xticks([])
#        yticks([])
#        scalefig(p[0], extent=0.15)
#        
#        savefig("d:/h5/initial_sim/%05i.png" % i, bbox_inches='tight', pad_inches=0.1)
#        close()
#        
#    for i,pack in enumerate(load_packing.loadPackings('d:/h5/N32~P1e-2~shear~FS.positions')):
#        print i
#        doPackingStuff(pack)
#        xlabel("step: % 3i" % i, fontname="Consolas", fontsize=30) 
#        xticks([])
#        yticks([])
#        scalefig(p[0], extent=0.15)
#        
#        savefig("d:/h5/fixedstep_sim/%05i.png" % i, bbox_inches='tight', pad_inches=0.1)
#        close()
#        
#    for i,pack in enumerate(load_packing.loadPackings('d:/h5/N32~P1e-2~shear~REAR.positions')):
#        print i
#        doPackingStuff(pack)
#        xlabel("step: % 3i" % i, fontname="Consolas", fontsize=30) 
#        xticks([])
#        yticks([])
#        scalefig(p[0], extent=0.15)
#        
#        savefig("d:/h5/rearstep_sim/%05i.png" % i, bbox_inches='tight', pad_inches=0.1)
#        close()
#    ion()
        
#    for i,pack in enumerate(load_packing.loadPackings('d:/h5/N32~P1e-2~shear~SmallStep.convergedpositions')):
#        print i
#        doPackingStuff(pack)
#        xlabel("step: % 3i" % i, fontname="Consolas", fontsize=30) 
#        xticks([])
#        yticks([])
#        scalefig(p[0], extent=0.15)
#        
#        savefig("d:/h5/rearstep_sim/%05i.png" % i, bbox_inches='tight', pad_inches=0.1)
#        close()
    
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
