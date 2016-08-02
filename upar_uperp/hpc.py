# -*- coding: utf-8 -*-
"""
Code to calculate linear response to a given deformation.

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
        """
        Create a new HessianPackingCalculator. group can either be 
        a pyTables group corresponding to a packing or a packing in
        dictionary format
        """
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
        """
        Deform a packing and calculate linear response.
        deformation=1 is pure shear, 2 is simple shear.
        gamma=1 is the magnitude of the applied deformation.

        Returns a dict
        {
            'L1def': L1 (= (Lxx, Lxy) ) deformation,
            'L2def': L2 (= (Lyx, Lyy) ) deformation, 
            'delta_x': Δx for each particle,
            'delta_y': Δy for each particle,
            'forces_ext': forces on boundaries,
            'energy_cost': energy cost
        }

        some internal calculation results are also exposed.
        """
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
        """
        Plot linear response of particles to an external deformation,
        as calculated by deform_packing.
        delta_x = deform_packing(...)['delta_x']
        delta_y = deform_packing(...)['delta_y']
        L1def = deform_packing(...)['L1def']
        L2def = deform_packing(...)['L2def']
        """
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
        """
        Determine uperp and upar for each particle pair for a given deformation.
        * deformation=1 is pure shear, =2 is simple shear.
        * including_non_bonds=True includes uperp/upar for particles not in contact;
          otherwise these entries are np.nan. 
        
        returns upar, uperp (both NxN arrays containing upar_{i,j} and uperp_{i,j}).
        """
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

        delta_xij = delta_xij - self.contacts['nx'] * L1def[0] - self.contacts['ny'] * L2def[0]
        delta_yij = delta_yij - self.contacts['nx'] * L1def[1] - self.contacts['ny'] * L2def[1]
        
        u_parr = (xij_hat * delta_xij) + (yij_hat * delta_yij)
        u_perp = sqrt(delta_xij**2 + delta_yij**2 - u_parr**2)

        return u_parr, u_perp
        
    def get_scaled_uparperps(self, deformation=1):
        """
        Return scaled upar/uperp as in
            Ellenbroek et al., Phys. Rev. E 80, 061307, doi:10.1103/PhysRevE.80.061307
            
        """
        assert(self.loaded)
        u_parr, u_perp = self.get_uparrs(deformation=deformation)

        delta_ij_ellenbroek = 1 / (1 + self.contacts['rij'] / self.contacts['dij'])
        
        u_parr_scaled = (u_parr / delta_ij_ellenbroek ** 0.25)
        u_perp_scaled = (u_perp * delta_ij_ellenbroek ** 0.25)
        
        u_parr_scaled = u_parr_scaled[isfinite(u_parr_scaled)]
        u_perp_scaled = u_perp_scaled[isfinite(u_perp_scaled)]
        
        return u_parr, u_parr_scaled, u_perp, u_perp_scaled

    def find_first_ccs(self):
        """
        Calculate the first contact change under simple shear from linear response.
        Two methods are used:
          - simple linear: γ = δ/u∥
          - full quadratic: solve (Ri + Rj)² = (δ + u∥ * γ)² + (u⊥ * γ)²

        Returns {'gmk_FQ': ...,  # smallest making strain, full quadratic
                 'gbk_FQ': ...,  # smallest breaking strain, full quadratic
                 'gmk_SL': ...,  # smallest making strain, simple linear
                 'gmk_FQ': ...   # smallest breaking strain, simple linear
                }
        Returns gamma_min for mk and bk (in that order)
           nan if not found..."""
        u_par, u_perp = self.get_uparrs(2, True)
        u_par = -u_par
        delta = self.contacts["dijfull"]
        Rij = self.contacts["rij"]

        # we use 2 methods to determine gamma
        # full quadratic
        a = u_par**2 + u_perp**2
        b = 2 * Rij * u_par
        c = -2 * Rij * delta - delta**2

        # + for bk
        gammas = (-b + sqrt(b**2 - 4*a*c)) / (2*a)
        null, gbk_FQ = self.find_ccs_from_gammas(gammas)

        # - for mk
        gammas = (-b - sqrt(b**2 - 4*a*c)) / (2*a)
        gmk_FQ, null = self.find_ccs_from_gammas(gammas)

        #gammas = (delta / u_par) * (1 - (delta / Rij) * (u_perp / u_par)**2)
        #gmk_SQ, gbk_SQ = self.find_ccs_from_gammas(gammas)
        #
        #gammas = (delta / u_par) * (1 - (delta / Rij) * (u_perp**2 / (u_perp**2 + u_par**2)))
        #gmk_SQ2, gbk_SQ2 = self.find_ccs_from_gammas(gammas)

        # simple linear
        gammas = (delta / u_par)
        gmk_SL, gbk_SL = self.find_ccs_from_gammas(gammas)

        return {k: v for k,v in locals().items() if k.startswith('gmk') or k.startswith('gbk')}

    def find_ccs_from_gammas(self,gammas):
        rattlers = self.contacts["rattlers"]
        connmatrix = self.contacts["connmatrix"] 

        gmk = np.nan
        gbk = np.nan

        gammas_f = gammas[gammas>0]
        order = argsort(gammas_f)

        for i in range(0, len(gammas_f), 2):
            idx = where(gammas == gammas_f[order[i]])
            if idx[0][0] not in rattlers and idx[0][1] not in rattlers:
                if connmatrix[idx][0] and not isfinite(gbk):
                    gbk = 2*gammas_f[order[i]]
                elif (not connmatrix[idx][0]) and not isfinite(gmk):
                    gmk = 2*gammas_f[order[i]]
                if isfinite(gmk) and isfinite(gbk):
                    break

        return gmk, gbk

    def get_el_con(self):
        """
        Calculate elastic constants via V_harm.El_Con
        (prevents re-calculating K)
        """
        el_con = V_harm.El_Con(self.Hess, self.rattlers, self.packing)
        c = el_con.pop("c")
        for i, element in enumerate(c):
            el_con["c{num}".format(num=i+1)] = element
        el_con["Galpha"] = self.deform_packing(deformation=2)['energy_cost']/self.packing["L"]**2/4
        return el_con
