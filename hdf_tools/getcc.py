# -*- coding: utf-8 -*-
"""
Created on Mon May 06 12:06:04 2013

Tool to get contact changes from an hdf group
"""
import time
from numpy import min, max, sum, diff, log10, abs, amin, amax, where

def get_first_ccs_base(group, data=None):
    if data is None:
        try:
            data = group.SR.data.read()
        except Exception, e:
            print e, " - sleeping to recover?"
            time.sleep(10)
            print "retrying...."
            data = group.SR.data.read()
        
    if sum(diff(data["step#"]) != 1):
        print group._v_pathname, diff(data["step#"])
        raise Exception("step ordering")        
    #data.sort(order="step#") # because this sometimes is not the case O_o.
    
    if min(data['Nchanges']) != 0:
        raise Exception('no zero-change state')
    
    if max(data["Nchanges"]) == 0:
        raise Exception('no contact change at all!')
        
    if round(log10(data[0]["gamma"]),2) == -16:
        # remove incorrect first entry
        data = data[1:]
        
    # find gamma0, sigma0
    try:
        lastconvergenceid = where(diff(log10(abs(diff(log10(data["gamma"]))))) > .3)[0][0]+1
        subdata = data[:lastconvergenceid+1]
    except Exception, e:
        if abs(diff(log10(abs(diff(log10(data["gamma"])))))[-1] + log10(2)) > 1e-6:
            print group._v_pathname, e, diff(log10(abs(diff(log10(data["gamma"]))))), "(38)"
        subdata = data[:]  
        
    return locals()
    
def get_first_cc(group, data=None):
    """Given a HDF group, returns the row (from the data array) just before and just after the CC"""
    exec ""
    locals().update(get_first_ccs_base(group, data))
    
    before = subdata[subdata["gamma"] == amax(subdata[subdata["Nchanges"] == 0]["gamma"])]
    before=before[0];
    after  = subdata[subdata["gamma"] == amin(subdata[subdata["Nchanges"] > 0]["gamma"])]
    after=after[0];

    return before, after  
    
def get_first_ccs(group, data=None):
    exec ""
    locals().update(get_first_ccs_base(group, data))
    subdata.sort(order=["gamma"])
    before = subdata[subdata["gamma"] <= amax(subdata[subdata["Nchanges"] == 0]["gamma"])]
    after  = subdata[subdata["gamma"] >= amin(subdata[subdata["Nchanges"] > 0]["gamma"])]

    return before, after   
