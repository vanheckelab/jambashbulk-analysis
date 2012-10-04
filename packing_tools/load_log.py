# -*- coding: utf-8 -*-
# Created on Fri Jun 15 12:31:25 2012
# @author: Merlijn van Deen <deen@physics.leidenuniv.nl>

"""
Functions to load packing log files.

Two formats are supported, both from Jo's code:
    - Old format (single runtime (s))
    - New format (double runtime (s))
"""

import os
import glob
import numpy as np
from fs_tools import getPrefix

def parseTime(formatted):
    from datetime import datetime
    return datetime.strptime(formatted, "%Y-%m-%d_%H-%M-%S").isoformat()

def parseOldLogLine(line):
    keys = ["PackingNumber", "N", "P0", "P", "alpha", "delta", "L",
              "phi", "Z", "N - Ncorrected", "sxx", "syy", "sxy", "Uhelper",
              "dU", "H", "dH", "runtime (s)", "FIRE iteration count",
              "CG iteration count", "maxGrad", "current time"]
    types = dict([(key, np.float64) for key in keys]) # default format for parameters
    types['current time'] = parseTime
    types['FIRE iteration count'] = np.int
    types['currentPackingNumber'] = np.int
    types['N'] = np.int
    types['CG iteration count'] = np.int
    
    return _parseLogLine(line, keys, types)
    
def parseNewLogLine(line):
    keys = ["PackingNumber", "N", "P0", "P", "alpha", "delta", "L",
              "phi", "Z", "N - Ncorrected", "sxx", "syy", "sxy", "Uhelper",
              "dU", "H", "dH", "runtime (s)", "FIRE iteration count",
              "CG iteration count", "maxGrad", "runtime (s)", "current time"]
    types = dict([(key, np.float64) for key in keys]) # default format for parameters
    types['current time'] = parseTime
    types['FIRE iteration count'] = np.int
    types['currentPackingNumber'] = np.int
    types['N'] = np.int
    types['CG iteration count'] = np.int
    
    return _parseLogLine(line, keys, types)

def _parseLogLine(line, keys, types):
    values = line.strip().split()
    mapped = dict(zip(keys, values))
    
    for key in mapped:
        mapped[key] = types[key](mapped[key]) # convert values
    
    return mapped

def parseLogLine(line):
    try:
        return parseNewLogLine(line)
    except ValueError:
        return parseOldLogLine(line)


def getLogLines(folder):
    prefix = getPrefix(folder)
    
    logfiles = glob.glob(os.path.join(folder, "log%s.txt" % prefix)) + \
               glob.glob(os.path.join(folder, "log%s~????.txt" % prefix))
    
    for logfile in logfiles:
        try:
            f = open(logfile)
        except Exception:
            import time
            time.sleep(30)
            f = open(logfile)
            
        for line in f.xreadlines():
	    if line.startswith('seed'):
	        continue
            yield line
        f.close()

def getParsedLogLines(folder):
    for logline in getLogLines(folder):
        yield parseLogLine(logline)
        
def getUniqueParsedLogLines(folder):
    seen = []
    key = lambda line: line['PackingNumber']
    for logline in getParsedLogLines(folder):
        if key(logline) in seen:
            continue
        seen.append(key(logline))
        yield logline
