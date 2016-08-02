# -*- coding: utf-8 -*-
# Created on Thu Jul 05 11:59:42 2012
# @author: Merlijn van Deen <deen@physics.leidenuniv.nl>

"""
Tools to plot CDFs and PDFs
"""

from scipy.stats import gaussian_kde
from numpy import std, array, float_, isfinite
from pylab import plot, linspace
import warnings

def densplot(data, range_, *args, **kwargs):
    """ Plot a density plot (PDF).

        data = array of elements to plot the PDF of
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
    """ Plot a cumulative distribution function (CDF).
    
    data = array of elements to calculate the CDF from.
    
    All other arguments are passed to plot().
    """
    x, y = get_cdf_data(data)
    return plot(x,y, *args, **kwargs)

def get_cdf_data(data, lower=None, upper=None):
    """ Build CDF from data.
    
    data = array of elements to calculate the CDF from,
    upper = maximum x-value to include (will be 1),
    lower = minimum x-value to include (will be 0).

    Returns: (x, y)
    x: numpy array, x values of the CDF
    y: numpy array, y values of the CDF
    """
    data = array(data)
    data = data[isfinite(data)]
    if data.size == 0:
        return
   
    length = len(data)
    maxdatastep = max(data)-min(data) #max(diff(data))
    if lower is None:
        lower = data[0]-maxdatastep
    if upper is None:
        upper = data[-1]+maxdatastep
    data = list(data) + list(data)
    x = sorted(data)
    x = [lower] + x + [upper]
    
    y = [0] + sorted(range(length) + range(length))[1:] + [length, length]
    y = float_(y)/length

    return x, y    
   
def plot_icdf(data, *args, **kwargs):
    """ Plot the inverse CDF (1-CDF) of `data`.
        data = array of elements to calculate the CDF from.
    
        All other arguments are passed to plot().
    """
    if data.size == 0:
        return
    length = len(data)
    maxdatastep = max(data)-min(data) #max(diff(data))
    data = list(data) + list(data)
    x = sorted(data)
    x = [data[0]-maxdatastep] + x + [data[-1]+maxdatastep]
    
    y = [0] + sorted(range(length) + range(length))[1:] + [length, length]
    y = float_(y)/length
    
    return plot(x,1-y, *args, **kwargs)