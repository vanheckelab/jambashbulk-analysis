# -*- coding: utf-8 -*-
# Created on Fri Jun 15 13:10:54 2012
# @author: Merlijn van Deen <deen@physics.leidenuniv.nl>

"""
Filesystem tools
"""

import os
import numpy as np

def getPrefix(folder):
    folder = os.path.realpath(folder)
    prefix = os.path.split(folder)[1]
    return prefix

