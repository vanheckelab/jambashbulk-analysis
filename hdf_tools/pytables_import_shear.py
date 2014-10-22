# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 13:33:29 2012

@author: deen
"""

import os, sys
sys.path.append(os.path.split(os.path.split(__file__)[0])[0])
import time
import glob
import itertools
import pandas
import numpy as np
from numpy import float64, array, dtype, argmin, isnan
import re
from packing_tools.load_packing import loadPackings 
 
import tables
from pytables_tools import require_group, require_table, add_to_table, store_table
from pytables_test import create_group_for_packing, insert_packing_parameters, insert_packing_particles

packing_attr_cache_dtype = np.dtype([
('L', dtype('float64')),
('L1', dtype('float64'), (2,)),
('L2', dtype('float64'), (2,)),
('N', dtype('int64')),
('P', dtype('float64')),
('P0', dtype('float64')),
('gamma', dtype('float64')),
('s_xy', dtype('float64')),
('Ncontacts', dtype('int64')),
('Nchanges', dtype('int64')),
('N+', dtype('int64')),
('N-', dtype('int64')),
('Z', dtype('float64')),
('step#', dtype('int64')),
('alpha', dtype('float64')),
('delta', dtype('float64')),
('phi', dtype('float64')),
('#rattler', dtype('int64')),
('s_xx', dtype('float64')),
('s_yy', dtype('float64')),
('U', dtype('float64')),
('dU', dtype('float64')),
('H', dtype('float64')),
('dH', dtype('float64')),
('t_run', dtype('int64')),
('#FIRE', dtype('int64')),
('#CG', dtype('int64')),
('gg', dtype('float64')),
('creation-date', dtype('|S19')),
('path', '|S128')])

insert_particles = True

def shear_type(st):
    if st.startswith("SR"):
        return "SR"
    else:
        return st

def groupkey(args):
    N,P,st,step,number = args[0]
    return N,P,number,shear_type(st)
    
def sortkey(args):
    return list(groupkey(args)) + args[0]

class CSVReadException(Exception):
    pass

def read_csv(fn, converters={}):
    try:
        f = open(fn)
        comments = ""
        
        while True:
            oldpos = f.tell()
            line = f.readline()
            
            if line.startswith('#'):
                comments += line
            else:
                f.seek(oldpos)
                break
        return pandas.read_csv(f, sep="[ \t]", converters=converters), comments
    except Exception, e:
        raise CSVReadException(e)

def load_measurement(base, spec):
    log, comments = read_csv(os.path.join(base, "log" + spec), converters={'P0': np.float64})
    print "log",
    sys.stdout.flush()
    
    data, comments = read_csv(os.path.join(base, "data" + spec))
    print "data",
    sys.stdout.flush()
    assert (len(log) == len(data)), "log and data have different lengths"

    # 'eta' (applied strain) was renamed to gamma
    if 'eta' in data and 'gamma' not in data:
        data['gamma'] = data['eta']
        data.pop('eta')
    
    # and then to gamma_alpha, but we will keep gamma as name here....
    if 'gamma_alpha' in data and 'gamma' not in data:
        data['gamma'] = data['gamma_alpha']
        data.pop('gamma_alpha')
    
    # gamma=10^-9 is actually gamma=0; fix this
    if any(log["step#"] == 0):
        data["gamma"][log["step#"] == 0] = 0
    
    # and I think there is also some bug with very low gamma (10^-16-ish) but I
    # have no clue what the exact problem was (or the fix)
    # (probably the first strain is not 0 there. Or something.)
    # :-( - Merlijn 4/2/13
    # 10^-16 was actually 10^-8, I think. Incorrect 'small gamma' correction where the data for 10^-8 was written with 10^-16 as listed strain

    if insert_particles:
        try:
            if os.path.exists(os.path.join(base, "particles" + spec)):
                print "particles",
                sys.stdout.flush()
                # new format, where we know the particle positions for every step
                particles = loadPackings(os.path.join(base, "particles" + spec))
    
                if len(particles) > len(data):
                    particles = [particles[0]] + particles[1::2]
                particles = pandas.DataFrame(particles)
                assert (len(data) == len(particles)), "data and particles have different lengths"
    
                #import pdb; pdb.set_trace()
                packings = data.join(particles, rsuffix="_").join(log, rsuffix="__")
                print len(packings)
            elif os.path.exists(os.path.join(base, spec)):
                print "oldparticles",
                sys.stdout.flush()

                # old format where we only know some particle positions
                packings = data.join(log, rsuffix="__")
                
                particles = loadPackings(os.path.join(base, spec))
                particles = pandas.DataFrame(particles)
                
                # create comparison columns: change in L2[0]
                particles["dL20"] = (array([x[0] for x  in particles["L2"]]) - particles["L2"][0][0])
                packings["dL20"] = packings["gamma"] * packings["L"]
                
                # from this, we can determine the correct step#s
                particles["step#"] = array(packings["step#"][argmin(abs(packings["dL20"][:, np.newaxis] - array(particles["dL20"])), axis=0)])
                packings = packings.merge(particles, how="left", left_on="step#", right_on="step#", suffixes=("", "_"))
            else:
                print "noparticles"
                packings = data.join(log, rsuffix="__")
        except IOError, e:
            raise
            print "(failed: %r)" % e
            sys.stdout.flush()
            packings = data.join(log, rsuffix="__")
            
    else:
        packings = data.join(log, rsuffix="__")
        
    return comments, packings, log

def process_measurement(f, key, base, m, spec):
    group = require_group(f,'/'.join(key))
    
#    spec = "~".join(m) + ".txt"
    print key, ":", spec,
      
    comments, packings, log = load_measurement(base, spec)

    group._v_attrs['comments'] = comments
    attrcache = require_table(group, 'data', packing_attr_cache_dtype,
                              expectedrows=len(log), chunkshape=(len(log),))
    for rowno, row in packings.iterrows():
        name = "%04i" % row["step#"]
        if name not in group._v_children.keys():
            subgroup = require_group(group, name)
            
            insert_packing_parameters(subgroup, row)
            add_to_table(attrcache, row, path=subgroup._v_pathname)
        else:
            subgroup = getattr(group, name)
            
        if insert_particles and \
           'particles' in packings and \
           'particles' not in subgroup and \
            not (isinstance(row['particles'], float) and isnan(row['particles'])):
            insert_packing_particles(subgroup, row)
    
    attrcache.flush()

def key_for_fn(fn):
  part1 = os.path.split(os.path.split(fn)[0])[1].split("~") 
  part2 = os.path.split(fn)[1][4:-4].split("~")[2:]
  return list(part1) + list(part2)

def main(*argv):
    argv = list(argv)
    global insert_particles
    if '-noparticles' in argv:
        argv.remove('-noparticles')
        insert_particles = False

    if len(argv) not in (2,3):
        print "Usage: %s <base from which to add shear data from> <h5 file to store shear data in> <start id (counter in file)> [-noparticles]"
        raise Exception()
        
    bbase = argv[0] #r"U:\ilorentz\simulations\Packings\N256~P1e-3"
    outfile = argv[1] #r"U:\ilorentz\simulations\Packings\N256~P1e-3\data.h5"
    if len(argv) == 3:
        start = int(argv[2])
    else:
        start = 0
    
    f = tables.openFile(outfile, mode = "a")
    root = f.root
    
    try:
        shear_measurements = [(key_for_fn(fn), fn)
                                for fn
                                in sorted(glob.glob(bbase + "*/data*.txt"))]
        for i,(key, measurements) in enumerate(itertools.groupby(sorted(shear_measurements, key=sortkey), groupkey)):
            print i, 
            if i < start:
                continue
            for m in measurements:
                pass
            print key

            ignores = {"N512": {"P3162e-1": [3, 18, 29, 53, 64, 82, 94, 95, 110, 122, 131,
                                             146, 176, 194, 196],
                               },
                       "N32":  {"P1e+0": [9605],},
                       "N64":  {"P3162e-1": [9105, 9231, 9521, 9827]}}

            if int(key[2]) in ignores.get(key[0], {}).get(key[1], []):
                print "Skipping (FORCED IGNORE)"
                continue
            try:
                print root, key, os.path.split(m[1])[0], m[0]
                process_measurement(root, key, os.path.split(m[1])[0], m[0], os.path.split(m[1])[1][4:])        
            except (ValueError, ), e:
                if str(e) == "cannot index with vector containing NA / NaN values"  or \
                   re.match(r"Expecting \d+ columns, got \d+ in row \d+", str(e)):
                    print key, repr(e)
                else:
                    raise
            except (AssertionError, ), e:
                if str(e) in ["log and data have different lengths", "data and particles have different lengths"]:
                    print key, repr(e)
                else:
                    raise
            except (IOError, EOFError, StopIteration, CSVReadException), e:
                print key, repr(e)
            except (Exception, ), e:
                if str(e) in ["Reindexing only valid with uniquely valued Index objects"]:
                    print key, repr(e)
                else:
                    raise
    finally:
        f.flush()
        f.close()
        print "rrr"
        
if __name__ == "__main__":
    try:
        sys.argv = argv
    except NameError:
        pass
    main(*sys.argv[1:])


#try:
#    for base in glob.glob(bbase + "*"):
#        for packing in getPackings(base):
#            path = createpyTablesNodeForPacking(root, packing)
#            add_to_table(attrcache, packing, path=path)
#finally:
#    f.flush()
#    f.close()    
#    
#    
