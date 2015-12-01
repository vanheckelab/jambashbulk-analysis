# -*- coding: utf-8 -*-
"""
Created on Mon May 06 12:06:04 2013

Tool to get contact changes from an hdf group
"""
import time
import warnings
from numpy import min, max, sum, diff, log10, abs, amin, amax, where
import numpy as np
import sys

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
        raise Exception(("step ordering", group._v_pathname, diff(data["step#"])))
    #data.sort(order="step#") # because this sometimes is not the case O_o.
    
    if min(data['Nchanges']) != 0 and not missingzeroOK:
        raise Exception('no zero-change state')
    
    if max(data["Nchanges"]) == 0:
        raise Exception('no contact change at all!')
        
    if round(log10(data[0]["gamma"]),2) == -16:
        # remove incorrect first entry
        data = data[1:]
    
    # We have two methods to find the end of this CC convergence. Mid-2014, an extra parameter
    # 'ccnum' was added, which allows us to directly select the data relating to this CC.
    # Before, we had to derive this from the convergence method. Throw a warning if we still have 
    # to do this.
    
    try:
        thiscc = data[0]["ccnum"]
        lastconvergenceid = where(data["ccnum"] == thiscc)[0][-1]
        subdata = data[:lastconvergenceid+1]
        
        return locals()
    except (KeyError, ValueError, IndexError) as e:
        # 'ccnum' not found
        pass
        
    
    warnings.warn("Using convergence derivation instead of ccnum")
    
    try:
        lastconvergenceid = where(diff(log10(abs(diff(log10(data["gamma"]))))) > .3)[0][0]+1
        subdata = data[:lastconvergenceid+1]
    except Exception, e:
        if abs(diff(log10(abs(diff(log10(data["gamma"])))))[-1] + log10(2)) > 1e-6:
            #print >>sys.stdout, group._v_pathname, e, diff(log10(abs(diff(log10(data["gamma"]))))), "(38)"
            pass
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

def get_multi_ccs(group, data=None, error_cutoff=1000):
    if data is None:
        data = get_first_ccs_base(group, data)["data"]
    
    ordata = np.recarray.copy(data)
    data = np.recarray.copy(data)
    
    
    id_ctr = 0
    stop = False
    
    after = None
    for i in range(error_cutoff):
        try:
            result = get_first_ccs_base(group, data)
            #print result
            lcid = id_ctr + result['lastconvergenceid'] + 1
        except KeyError:
            # no convergence found, assume = last id of array
            lcid = id_ctr + len(data)
            stop = True

        subdata = result["subdata"]
        print data["Nchanges"], subdata["Nchanges"]
        
        
        if amin(subdata["N+"] + subdata["N-"]) > 0:
            # The base state is not in our data set! Add the previous after state and try from there.
            if after:
                subdata = np.hstack([after, subdata])
                # and now re-set Nchanges based on our newly-found base state...
                subdata["Nchanges"] -= subdata[0]["Nchanges"]
            else:
                raise Exception("Help! We don't know the base state, and there is no previous after state...")
        
        before = subdata[subdata["gamma"] == amax(subdata[subdata["Nchanges"] == 0]["gamma"])]
        before=before[0];
        before = ordata[ordata["index"] == before["index"]][0]
        
        after  = subdata[subdata["gamma"] == amin(subdata[subdata["Nchanges"] > 0]["gamma"])]
        after=after[0];
        after = ordata[ordata["index"] == after["index"]][0]
        
        yield (before, after)
        
        if not stop:

            data = ordata[lcid:]
            if len(data) == 0:
                # no data to view, so we should also break...
                break
            id_ctr = lcid
            
            # reset Nchanges to new base state
            sub = amin(data["Nchanges"])
            
            data["Nchanges"] = data["Nchanges"] - sub
            # also subtract from after, as we'll use that in the next loop
            after["Nchanges"] = after["Nchanges"] - sub
        else:
            break

    else:
        raise Exception("Error cutoff triggered after %i iterations" % i)