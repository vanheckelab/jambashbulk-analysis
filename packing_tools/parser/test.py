fn = "particlesN256~P1e-3~SR010~step021~9212.txt"

import sys; sys.path.append('/home/valhallasw/src/phd-library/'); from packing_tools.load_packing import loadPackings
import time
import parser

print time.time()
starttime = time.time()
loadPackings(fn)
print time.time(), time.time()-starttime
starttime = time.time()
[x for x in parser.read_packings(fn)]
print time.time(), time.time()-starttime
