# -*- coding: utf-8 -*-
# Created on Thu Jul 05 11:59:42 2012
# @author: Merlijn van Deen <deen@physics.leidenuniv.nl>

"""

"""

from scipy.stats import gaussian_kde
from numpy import std, float16
from pylab import plot, linspace
import warnings

def densplot(data, range_, *args, **kwargs):
    """ data = array of elements to plot the PDF of
        range_ = [xmin, xmax]
        other arguments are passed to the plot function
    """
    if std(data) == 0:
        warnings.warn("std(data) is 0, cannot plot PDF")
        return
    density = gaussian_kde(data)
    xs = linspace(*range_)
    plot(xs, density(xs), *args, **kwargs)
    
def plot_cdf(data, *args, **kwargs):
    if data.size == 0:
        return
    length = len(data)
    maxdatastep = max(data)-min(data) #max(diff(data))
    data = list(data) + list(data)
    x = sorted(data)
    x = [data[0]-maxdatastep] + x + [data[-1]+maxdatastep]
    
    y = [0] + sorted(range(length) + range(length))[1:] + [length, length]
    y = float16(y)/length
    
    plot(x,y, *args, **kwargs)
    
def plot_icdf(data, *args, **kwargs):
    if data.size == 0:
        return
    length = len(data)
    maxdatastep = max(data)-min(data) #max(diff(data))
    data = list(data) + list(data)
    x = sorted(data)
    x = [data[0]-maxdatastep] + x + [data[-1]+maxdatastep]
    
    y = [0] + sorted(range(length) + range(length))[1:] + [length, length]
    y = float16(y)/length
    
    plot(x,1-y, *args, **kwargs)