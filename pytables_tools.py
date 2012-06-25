# -*- coding: utf-8 -*-
# Created on Fri Jun 15 13:59:04 2012
# @author: Merlijn van Deen <deen@physics.leidenuniv.nl>

"""
pyTables helper functions
"""

import tables

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
    return root._v_file.createTable(root, name, data, expectedrows=data.shape[0], chunkshape=data.shape)

def add_to_table(table, data={}, **kwargs):
    data = data.copy()
        
    row = table.row
    for key in table.colnames:
        if key in data:
            row[key] = data[key]
        if key in kwargs:
            row[key] = kwargs[key]
            
    row.append()