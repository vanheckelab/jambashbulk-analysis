from pylab import *
import calculate_contacts, load_packing

packing = load_packing.loadPackings('N32~P1e-3~9000.txt')[0]
conts = calculate_contacts.get_contacts(packing)

print "Rattlers: ", where(sum(conts['connmatrix'], axis=0) == 0)

dij = conts['dij']

print "Maximum overlap: ", nanmax(dij), " (between particles %i and %i) " % unravel_index(nanargmax(dij), dij.shape)
