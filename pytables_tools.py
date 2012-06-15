# -*- coding: utf-8 -*-
# Created on Fri Jun 15 13:59:04 2012
# @author: Merlijn van Deen <deen@physics.leidenuniv.nl>

"""
pyTables helper functions
"""

import tables

def require_group(root, name, *args, **kwargs):
    h5f = root._v_file
    try:
        return root._f_getChild(name)
    except tables.exceptions.NoSuchNodeError:
        return h5f.createGroup(root, name, *args, **kwargs)

def require_table(root, name, dtype):
    try:
        return root._f_getChild(name)
    except tables.exceptions.NoSuchNodeError:
        return root._v_file.createTable(root, name, dtype)

def store_table(root, name, data, *args, **kwargs):
    return root._v_file.createTable(root, name, data, expectedrows=data.shape[0], chunkshape=particles.shape)

def add_to_table(table, data={}, **kwargs):
    data = data.copy()
    data.update(kwargs)
        
    row = table.row
    for key in table.colnames:
        if key in data:
            row[key] = data[key]
    row.append()