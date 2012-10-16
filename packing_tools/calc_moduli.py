import tables
import V_harm

import sys

if len(sys.argv) < 2:
    print "Usage: calc_moduli.py file.h5"

t = tables.File(sys.argv[1])
for group in iter(t.root._v_groups.items()[0][1]).next():
    packing = dict((x, group._v_attrs[x]) for x in group._v_attrs._v_attrnames)
    packing['particles'] = group.particles.read()
    
    contacts = V_harm.get_contacts(packing)

    info = {'N': packing['N'], 'P0': packing['P0'], 'PackingNumber': int(packing['PackingNumber'])}
    info.update(V_harm.V_harm(contacts, packing))
    info.update(V_harm.grad(contacts, packing))
    K = V_harm.Hess(contacts, packing)

    info.update(V_harm.El_Con(K, contacts['rattlers'], packing))
    print repr(info).replace("\n", "\t")
