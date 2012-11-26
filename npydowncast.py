import sys
import numpy as np
from numpy import dtype, float64, float128

m = {dtype(float128): dtype(float64)}
data = np.load(sys.argv[1])
newdata = data.astype([(name, m.get(dt,dt)) for name, (dt, alignment) in sorted(data.dtype.fields.items(), key=lambda x: x[1][1])])
np.save(sys.argv[1] + ".reduced.npy", newdata)
