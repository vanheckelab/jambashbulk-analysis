import numpy as np
from numpy import floor, sqrt, sum, tril, where
from copy import copy

def get_contacts(packing):
    """ In: packing={'particles': {'x': [...], 'y': [...], 'r': [...]},
                    {'alpha': ..., 'delta': ..., 'L': ...}}

        Out: {xij, yij, connmatrix (= Cij)}
    """
    particles = packing['particles']
    alpha, delta = packing['alpha'], packing['delta']
    L = packing['L']
       
    xi, xj = np.meshgrid(particles['x'], particles['x'])
    yi, yj = np.meshgrid(particles['y'], particles['y'])

    lyy = L*(1+delta)
    lxx = L/(1+delta)
    lyx = L*alpha
    lxy = 0
    
    xij = xj-xi
    yij = yj-yi
    
    ny = -floor((yij + 0.5 * lyy)/lyy)
    
    xij += ny * lyx
    yij += ny * lyy
    
    nx = -floor((xij + 0.5 * lxx)/lxx)
    
    xij += nx * lxx
    yij += nx * lxy
    
    rij = sqrt(xij**2 + yij**2)

    dij = (particles['r'] + particles['r'][:,np.newaxis]) - rij
        
    connmatrix = dij > 0
    np.fill_diagonal(connmatrix, False)

    while True:
        Ncontacts = sum(connmatrix, axis=0)
        rattlers = where((Ncontacts > 0) * (Ncontacts <= 2))[0]

        if len(rattlers) == 0:
            break

        connmatrix[rattlers,:] = False
        connmatrix[:, rattlers] = False

    dij[~connmatrix] = 0
    return {'xij': xij, 'yij': yij, 'dij': dij, 'connmatrix': connmatrix} 
