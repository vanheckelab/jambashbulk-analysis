# -*- coding: utf-8 -*-
"""
Smart annotation in MPL plots

Created on Fri Sep 27 14:52:59 2013

@author: Merlijn van Deen
"""
import pylab
from pylab import argmin, array, sum

def pick_handler(event):
    artist = event.artist
    mouseevent = event.mouseevent
    
    i = argmin(sum((array(artist.get_data()).T - array([mouseevent.xdata, mouseevent.ydata]))**2, axis=1))
    
    x,y = array(artist.get_data()).T[i]
    
    Nelements = len(artist.get_data()[0])
    d = str(artist.get_url()[i]) if Nelements > 1 else str(artist.get_url())
    
    ax = artist.axes
    
    if not hasattr(ax, "annotation"):
        ax.annotation = ax.annotate(d,
            xy=(x,y), xycoords='data',
            xytext=(0,30), textcoords="offset points",
            size="larger", va="bottom", ha="center",
            bbox=dict(boxstyle="round", fc="w", alpha=0.5),
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3,rad=-0.2"), 
        )
        
    ax.annotation.xy = (x,y)
    ax.annotation.set_text(d)
    artist.figure.canvas.draw()

def plot(*args, **kwargs):
    ax = kwargs.pop("ax", pylab.gca())
    
    if "ident" in kwargs:
        kwargs["url"] = kwargs.pop("ident")
        kwargs["picker"] = kwargs.get("picker", 10)
        canvas = ax.get_figure().canvas
        if not hasattr(canvas, "annotation_handler_connected"):
            ax.get_figure().canvas.mpl_connect('pick_event', pick_handler)
            canvas.annotation_handler_connected = True
        
    return ax.plot(*args, **kwargs)
    
plot.__doc__ = pylab.plot.__doc__ + """
ident = label (or list of labels, for multiple points) to identify points using
        the mouse annotation picker
"""