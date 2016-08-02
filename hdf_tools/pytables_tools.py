# -*- coding: utf-8 -*-
# Created on Fri Jun 15 13:59:04 2012
# @author: Merlijn van Deen <deen@physics.leidenuniv.nl>

"""
pyTables helper functions

Most of these are helper functions for `pytables_test.py` and
`pytables_import_shear.py`. The exception to this is `read_packing`,
which loads a packing from the hdf5 file to a dictionary, as used
by most other code.
"""

import tables
import numpy

def require_groups(root, name, *args, **kwargs):
    names = name.split("/")
    return reduce(lambda x,y: require_group(x,y,*args,**kwargs), names, root)

def require_group(root, name, *args, **kwargs):
    if '/' in name:
        return require_groups(root, name, *args, **kwargs)
    h5f = root._v_file
    try:
        return root._f_getChild(name)
    except tables.exceptions.NoSuchNodeError:
        return h5f.createGroup(root, name, *args, **kwargs)

def require_table(root, name, dtype, *args, **kwargs):
    try:
        return root._f_getChild(name)
    except tables.exceptions.NoSuchNodeError:
        return root._v_file.createTable(root, name, dtype, *args, **kwargs)

def store_table(root, name, data, *args, **kwargs):
    data = numpy.array(data)
    return root._v_file.createTable(root, name, data, expectedrows=data.shape[0], chunkshape=data.shape)

def add_to_table(table, data={}, **kwargs):
    row = table.row
    for key in table.colnames:
        if key in data:
            row[key] = data[key]
        elif key in kwargs:
            row[key] = kwargs[key]
        else:
            if table.coldtypes[key] == numpy.string_:
                row[key] = "(unknown)"
            else:
                row[key] = numpy.array(0, dtype=table.coldtypes[key]) * numpy.nan
            
    row.append()

def read_packing(pack):
    """
    in: pack -- pytables node containing a packing
       (e.g. `/N1024/P3.1620e-03/0090` in N1024~P3162e-3_tables.h5,
       or `/N1024/P3.1620e-03/0090/SR/0000` in N1024~P3.162e-3_shear.h5)

    >>> import tables
    >>> f = tables.File(r"N256~P1e-3_tables.h5")
    >>> node = f.root.__getattr__('N256').__getattr__('p1.0000e-03').__getattr__('9000')
    >>> read_packing(node)
    { (...),
        'L': 37.448300000000003,
        'L1': array([ 37.28610578,   0.        ]),
        'L2': array([ -0.19961572,  37.61113632]),
        'N': 256,
        ...
        'particles': array([ (33.51530087911755, 2.706168622523819e-16, 25.13120589237754, -4.2327252813834093e-16, 1.0),
                             ...
                          ]),
       (...)
     }
    """
    packing = dict((x, pack._v_attrs[x]) for x in pack._v_attrs._v_attrnames)
    packing['particles'] = pack.particles.read()
    return packing
