# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 12:07:12 2013

@author: Merlijn van Deen
"""

import make_shear_graphs
import load_packing
import V_harm

print "Loading..."
pack = load_packing.loadPacking('d:/h5/xconf/Xconf_2011-07-31_5.txt')

#print "Plotting..."
#make_shear_graphs.doPackingStuff(pack)
#
#print "Calculating G et al..."
#conts = V_harm.get_contacts(pack)
#K = V_harm.Hess(conts, pack)
#elc = V_harm.El_Con(K, conts['rattlers'], pack)
#V_harm.Pmod(elc)

from hpc import HessianPackingCalculator
HPC = HessianPackingCalculator(None, pack)
u_parr, u_parr_scaled, u_perp, u_perp_scaled = HPC.get_scaled_uparperps()

print "saving... runtime HPC: ", time.time() - start, " s"
#np.savez(npfn,
#         u_parr = u_parr, u_parr_scaled = u_parr_scaled,
#         u_perp = u_perp, u_perp_scaled = u_perp_scaled,
#         K0_inv = HPC.K0_inv,
#         rij = HPC.contacts['rij'], dij = HPC.contacts['dij'])