# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 11:35:00 2013

@author: Merlijn van Deen
"""
import os
from pylab import *

def NPmapper(N,P, colormap='Accent'):
    # for N, use symbols:
    marker = {16: "o", 22: "v", 32: "D", 64: "^", 128: "p", 256: "*", 512: "s", 1024: "<", 2048: "d", 4096: ">"}.get(N, ".")
    
    # for P, use the hsv color map
    cmap = get_cmap('Accent')
    
    # make sure we are not outside the bound of the color map
    Pmin, Pmax = 1e-7, 1e-1
    P = amin(amax(P, Pmin), Pmax)
    
    # now normalize P logarithmically to 0..1
    normedP = (log10(P) - log10(Pmin)) / ((log10(Pmax) - log10(Pmin)) * 1+1e-9)
    
    color = cmap(normedP)
    
    return {'marker': marker, 'color': color}
    
def getmappeddrive(substr):
    drives = {}
    for lineset in os.popen("mountvol").read().split("\n\n"):
        if lineset.strip().startswith('\\\\?\\Volume'):
            vol, path = lineset.split("\n")
            vol = vol.strip()
            path = path.strip()
            if path.startswith("***"):
                continue
            if substr in vol:
                return path
    raise Exception('drive not found')
    
ESPCIDrive = lambda: getmappeddrive("{07c40d6f-c442-11e1-b08d-9cb70dd0f78e}")
ToshibaUSB3 = lambda: getmappeddrive("{a2449b94-7082-11e2-a88c-028037ec0200}")