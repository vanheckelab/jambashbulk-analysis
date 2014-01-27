# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 11:32:09 2013

@author: Merlijn van Deen
"""

import time
from pylab import *
import numpy as np
from packing_tools import V_harm
from packing_tools.make_shear_graphs import plotparticles, plotunitcell

def subset(dictionary, keys):
    return dict((k, v) for (k,v) in dictionary.iteritems() if k in keys)

class HessianPackingCalculator(object):
    loaded = False
    def __init__(self, group):
        self.loaded = True
        try:
            self.packing = dict((x, group._v_attrs[x]) for x in group._v_attrs._v_attrnames)
            self.packing['particles'] = group.particles.read()
        except AttributeError:
            self.packing = group
        
        self.contacts = V_harm.get_contacts(self.packing)
        self.Hess = np.float64(V_harm.Hess(self.contacts, self.packing))
        self.K_ext = self.Hess[:,:]
        self.rattlers = self.contacts['rattlers']
        
        #V_harm.Pmod(V_harm.El_Con(K_ext,rattlers,packing))
        
        xDOF = lambda num: num
        yDOF = lambda num: num + self.packing['N']
        
        # Remove rattler DOF from K_extended
        self.K_ext = np.delete(self.K_ext, np.append(xDOF(self.rattlers), yDOF(self.rattlers)), 0)
        self.K_ext = np.delete(self.K_ext, np.append(xDOF(self.rattlers), yDOF(self.rattlers)), 1)
        
        # Determine K0 from K_extended: remove the last four DOF = boundary DOF
        self.K0 = self.K_ext[:-4,:-4]
        if self.K0.size == 0:
            raise Exception("0x0 matrix after removing rattlers")
        
        # Diagonalize K: K = Q Λ Q⁻¹
        # Λ = diag(eigenvalues)
        
        [self.eigenvalues, self.Q] = eigh(self.K0)
        
        # K⁻¹ = Q Λ⁻¹ Q⁻¹, but we have λ=0 eigenvalues (global translations, zero energy
        # cost. We want these DOF to stay constant under stress, so we set λ⁻¹ = 0 for
        # them.
        self.eigenvalues_inv = 1/self.eigenvalues
        self.eigenvalues_inv[:2] = 0
        
		# this is *slow* on opterons. For some reason. Or maybe it's the memory. Whatever.
        #Lambda_inv = diag(eigenvalues_inv)
        #self.K0_inv = dot(self.Q, dot(Lambda_inv, self.Q.T))
        
        # Where Q = [v1, v2, v3, v4, ...] (stacked column vectors)
        # and Q⁻¹ = Q.T
        
        self.L1 = self.packing["L1"]
        self.L2 = self.packing["L2"]
        self.L = sqrt(self.L1[0]*self.L2[1])

    def deform_packing(self, deformation=1, gamma=1):
        assert(self.loaded)
        deformation_imposed_ext = zeros(self.K_ext.shape[0])
        
        if deformation == 1:
            # Deformation 1: pure shear 45°    
            L1def = array([self.L1[1], self.L1[0]]) * gamma  # delta L1
            L2def = array([self.L2[1], self.L2[0]]) * gamma  # delta L2
        elif deformation == 2:
            # Deformation 1: simple shear, note the factor 2...
            L1def = array([0,0]) * gamma * 2       # delta L1
            L2def = array([self.L2[1], 0]) * gamma * 2  # delta L2
            
        deformation_imposed_ext[-4:-2] = L1def
        deformation_imposed_ext[-2:]   = L2def
        
        # first, we calculate the forces on the particles using the extended hessian
        forces_ext = dot(self.K_ext, deformation_imposed_ext)
        
        # now, we want to calculate the movement of the particles that will result in
        # zero net force on the particles. For that, we want to know the movement
        # of the particles that will result in -forces_ext.
        forces_particles = forces_ext[:-4]
        
        #movement_particles = dot(self.K0_inv, -forces_particles) # note the -
        movement_particles = self.Q.dot(self.Q.T.dot(-forces_particles) * self.eigenvalues_inv)
        
        delta_x = movement_particles[:len(movement_particles)/2]
        delta_y = movement_particles[len(movement_particles)/2:]
        
        # re-insert rattlers
        rattler_insert_locations = self.rattlers - arange(len(self.rattlers))
        delta_x = insert(delta_x, rattler_insert_locations, zeros_like(rattler_insert_locations))
        delta_y = insert(delta_y, rattler_insert_locations, zeros_like(rattler_insert_locations))    
    
        # We now merge this into the imposed strain
        deformation_imposed_ext[:-4] += movement_particles[:]
        
        # and we can calculate the resulting forces
        forces_ext = dot(self.K_ext, deformation_imposed_ext)
        
        # and the energy cost
        energy_cost = dot(deformation_imposed_ext.T, forces_ext)
        
        return locals()

    def plot_deformation(self, delta_x, delta_y, L1def, L2def):
        assert(self.loaded)
        subplot(111, aspect="equal")
        axis((-self.L*0.3, self.L*1.3, -self.L*0.3, self.L*1.3))
        
        plotparticles(self.packing, fc="#dddddd")
        for offset_mult in unique(itertools.permutations([-1, -1, 0, 1, 1], 2)):
            offset = self.L1*offset_mult[0] + self.L2*offset_mult[1]
            plotparticles(self.packing, offset=offset, color="gray", fc="#eeeeee")
    
        factor = 2/self.L2[1] # scale by the height of the box = movement expected 
                           # due to fixed amount of shear
        
        plotunitcell({'L1': self.L1, 'L2': self.L2})
        plotunitcell({'L1': self.L1 + L1def * factor, 'L2': self.L2 + L2def * factor}, color="black")
        
        for (x,xerr,y,yerr,r), dx, dy in zip(self.packing['particles'], delta_x, delta_y):
            arrow(x,y,dx*factor,dy*factor, head_width=0.3)

    plot_deformation_keys = ["delta_x", "delta_y", "L1def", "L2def"]
    
   
    def plot_pure_shear(self):
        base = self.deform_packing(1)
        self.plot_deformation(**subset(base, self.plot_deformation_keys))
        
    def plot_simple_shear(self):
        base = self.deform_packing(2)
        self.plot_deformation(**subset(base, self.plot_deformation_keys))   
    
    def plot_shear_difference(self):
        base = deform_packing(1)
        other = deform_packing(2)
        
        params = {}
        for key in self.plot_deformation_keys:
            params[key] = base[key] - other[key]
        self.plot_deformation(**params)    

    def get_uparrs(self, deformation=1, including_non_bonds=False):
        assert(self.loaded)
        base = self.deform_packing(deformation)
        
        xij_hat = self.contacts['xij'] / self.contacts['rij']
        yij_hat = self.contacts['yij'] / self.contacts['rij']
        
        connmatrix = self.contacts['connmatrix']
        if not including_non_bonds:
            xij_hat[~connmatrix] = nan
            yij_hat[~connmatrix] = nan
        
        delta_x = base['delta_x']
        delta_y = base['delta_y']
        
        delta_xij = delta_x - delta_x[:,np.newaxis]
        delta_yij = delta_y - delta_y[:,np.newaxis]
        
        L1def = base["L1def"]
        L2def = base["L2def"]        
        
        delta_xij = delta_xij - (self.contacts['nx'] * connmatrix) * L1def[0] - (self.contacts['ny'] * connmatrix) * L2def[0]
        delta_yij = delta_yij - (self.contacts['nx'] * connmatrix) * L1def[1] - (self.contacts['ny'] * connmatrix) * L2def[1]        
        
        u_parr = (xij_hat * delta_xij) + (yij_hat * delta_yij)
        u_perp = sqrt(delta_xij**2 + delta_yij**2 - u_parr**2)

        return u_parr, u_perp
        
    def get_scaled_uparperps(self, deformation=1):
        assert(self.loaded)
        u_parr, u_perp = self.get_uparrs(deformation=deformation)

        delta_ij_ellenbroek = 1 / (1 + self.contacts['rij'] / self.contacts['dij'])
        
        u_parr_scaled = (u_parr / delta_ij_ellenbroek ** 0.25)
        u_perp_scaled = (u_perp * delta_ij_ellenbroek ** 0.25)
        
        u_parr_scaled = u_parr_scaled[isfinite(u_parr_scaled)]
        u_perp_scaled = u_perp_scaled[isfinite(u_perp_scaled)]
        
        return u_parr, u_parr_scaled, u_perp, u_perp_scaled

    def find_first_ccs(self):
        """Returns gamma_min for mk and bk (in that order)
           nan if not found..."""
        u_par, u_perp = self.get_uparrs(2, True)
        delta = self.contacts["dijfull"]
        Rij = self.contacts["rij"]
        gammas = -(delta / u_par) * (1 - (delta / Rij) * (u_perp / u_par)**2)
        
        try:
            gmk = amin(gammas[(gammas>0) * (delta < 0)])
        except ValueError:
            gmk = np.nan
        
        try:
            gbk = amin(gammas[(gammas>0) * (delta > 0)])
        except ValueError:
            gbk = np.nan
            
        return gmk*2, gbk*2  # epsilon, gamma, you know the drill...

    def get_el_con(self):
        el_con = V_harm.El_Con(self.Hess, self.rattlers, self.packing)
        c = el_con.pop("c")
        for i, element in enumerate(c):
            el_con["c{num}".format(num=i+1)] = element
        el_con["Galpha"] = self.deform_packing(deformation=2)['energy_cost']/self.packing["L"]**2/4
        return el_con