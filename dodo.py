import glob
import os
import sys

import tables
import numpy as np

bp = '/home/valhallasw/src/phd-library/'
sys.path.append(bp)
sys.path.append(bp.replace('valhallasw', 'merlijn')

from hdf_tools.pytables_test import main as import_tables
from hdf_tools.pytables_import_shear import main as import_shear
import hdf_tools.generate_cache
from hdf_tools.generate_cache import store_data_rows as generate_cache
from hdf_tools import getcc
from upar_uperp.hpc import HessianPackingCalculator
from upar_uperp.dodo_HPC import RunOnH5File as calc_linear_response

import itertools


def getNPdirpairs():
    dirs = glob.glob('SOURCES/0??? */N*~P*')
#    dirs = [d for d in dirs if len(d.split("/")[-1].split("~")[0]) == 3]
    k = lambda x: os.path.split(x)[-1]

    for NP, pack_dirs in itertools.groupby(sorted(dirs, key=k), key=k):
        yield NP, sorted(pack_dirs, reverse=True)

def task_read_base_packings():
    globaldeps = [os.path.join(bp, 'hdf_tools', 'pytables_test.py'),
                  os.path.join(bp, 'packing_tools', 'load_packing.py')]

    for NP, pack_dirs in getNPdirpairs():
        dependencies = globaldeps[:]
        for pack_dir in pack_dirs:
            dependencies.extend(glob.glob(os.path.join(pack_dir, NP + "~????.txt")))

        target = os.path.join('auto/h5', NP + "_tables.h5")

        # it would be better to split this in multiple temporary files, then merging them.
        yield {'basename': target,
               'targets': [target],
               'file_dep': dependencies,
               'actions': ["rm -f "+target] + [(import_tables, (pd, target)) for pd in pack_dirs]}

def task_read_shear_packings():
    globaldeps = [os.path.join(bp, 'hdf_tools', 'pytables_import_shear.py'),
                  os.path.join(bp, 'hdf_tools', 'pytables_tools.py'),
                  os.path.join(bp, 'hdf_tools', 'pytables_test.py'),
                  os.path.join(bp, 'packing_tools', 'load_packing.py')]

    for NP, pack_dirs in getNPdirpairs():
        dependencies = globaldeps[:]
        for pack_dir in pack_dirs:
            dependencies.extend(glob.glob(os.path.join(pack_dir, "log" + NP + "~SR*.txt")))
            dependencies.extend(glob.glob(os.path.join(pack_dir, "data" + NP + "~SR*.txt")))
            dependencies.extend(glob.glob(os.path.join(pack_dir, "particles" + NP + "~SR*.txt")))

        dependencies.append(os.path.join('auto/h5', NP + "_tables.h5"))
        target = os.path.join('auto/h5', NP + "_shear.h5")

        # it would be better to split this in multiple temporary files, then merging them.
        yield {'basename': target,
               'targets': [target],
               'file_dep': dependencies,
               'actions': ["rm -f " + target] + [(import_shear, (pd, target)) for pd in pack_dirs]}

def task_summarize_shear():
    dependencies = [os.path.join(bp, 'hdf_tools', 'generate_cache.py'),
                    os.path.join(bp, 'hdf_tools', 'getcc.py'),
                    os.path.join(bp, 'hdf_tools', 'pytables_tools.py')]

    target = 'auto/h5/shear_summary_cache.h5'
    mask = 'auto/h5/*_shear.h5'

    for NP, pack_dirs in getNPdirpairs():
        dependencies.append(os.path.join('auto/h5', NP + "_shear.h5"))
    for NP, pack_dirs in getNPdirpairs():
        dependencies.append(os.path.join('auto/h5', NP + "_tables.h5"))
    for NP, pack_dirs in getNPdirpairs():
        dependencies.append(os.path.join('auto/linres', NP + "_linres.npy"))

    hdf_tools.generate_cache.ignore_errnos.append(2)
    return {'basename': target,
            'targets': [target],
            'file_dep': dependencies,
            'actions': ["rm -f " + target, (generate_cache, (mask, "auto/linres/", target))]}

def task_linear_response():
    globaldeps = [os.path.join(bp, 'packing_tools', 'V_harm.py'),
                  os.path.join(bp, 'upar_uperp', 'hpc.py'),
                  os.path.join(bp, 'upar_uperp', 'dodo_HPC.py')]
    for NP, pack_dirs in getNPdirpairs():
        dependencies = globaldeps[:]
        src = os.path.join('auto/h5', NP + "_shear.h5")
        dependencies.append(src)
        target = os.path.join('auto/linres', NP + '_linres.npy')
        lrdatadir = os.path.join('auto/linres_upps', NP)

        yield {'basename': target,
               'targets': [target],
               'file_dep': dependencies,
               'actions': [(calc_linear_response, (src, target, lrdatadir))]}

#DOIT_CONFIG = {'default_tasks': ['auto/h5/shear_summary_cache.h5', task_summarize_shear]}

x = task_summarize_shear()

print "hdf_tools.generate_cache.ignore_errnos.append(2)"
print "generate_cache", x['actions'][1][1]
