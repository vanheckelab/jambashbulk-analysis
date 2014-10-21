import os
from ctypes import *

if os.name == "nt":
    from numpy import float64, float64 as longdouble, array
    libc = cdll.msvcrt 
    cparser = CDLL(os.path.join(os.path.split(__file__)[0], "_parser.dll"))
    c_longdouble = c_double
    def raise_error():
        raise WinError()
        
else:
    from numpy import float64, float128 as longdouble, array
    libc = cdll.LoadLibrary("libc.so.6")
    cparser = CDLL(os.path.join(os.path.split(__file__)[0], "_parser.so"))
    
    def raise_error():
        errno = c_int.in_dll(libc, "errno").value
        raise IOError(errno, libc.strerror(errno))       

import numpy as np

import os



fopen = libc.fopen
fclose = libc.fclose

libc.strerror.restype = c_char_p

class Header(Structure):
    _fields_ = [("N", c_int),
                ("L", c_longdouble),
                ("L1x", c_longdouble),
                ("L1y", c_longdouble),
                ("L2x", c_longdouble),
                ("L2y", c_longdouble),
                ("P", c_longdouble),
                ("P0", c_longdouble)
               ]

    def __repr__(self):
        jetzers = []
        for f,c in self._fields_:
            jetzers.append("%s=%s" % (f, getattr(self, f)))
        return "; ".join(jetzers)

def create_packing(H, particles):
    x = particles[:,0]
    y = particles[:,1]
    r = particles[:,2]

    x_major = float64(x); x_minor = float64(x-longdouble(x_major))
    y_major = float64(y); y_minor = float64(y-longdouble(y_major))
    r = float64(r)

    particles = array(zip(x_major, x_minor, y_major, y_minor, r), dtype=[('x', float64), ('x_err', float64),
                                                                         ('y', float64), ('y_err', float64),
                                                                         ('r', float64)])

    return {'P0': H.P0,
     'L': H.L,
     'N': H.N,
     'L1': array([H.L1x, H.L1y]),
     'L2': array([H.L2x, H.L2y]),
     'P': H.P,
     'P0': H.P0,
     'particles': particles}

def read_packings(fn):
    fptr = fopen(fn, "r")
    if not fptr:
        raise_error()

    try:
        h = Header()
        while(cparser.read_header(fptr, byref(h)) == 8):
#            print ".",
            import sys; sys.stdout.flush()
            particles = np.zeros([h.N, 3], dtype=longdouble)
            retval = cparser.read_particles(fptr, particles.ctypes.data_as(POINTER(c_longdouble)))
            if retval == 0:
                yield create_packing(h, particles)
            elif retval == c_int.in_dll(cparser, "_eof").value:
                raise EOFError()
            else:
                raise Exception("Unknown error in read_particles; return value was %i" % retval)
    finally:
        fclose(fptr)

if __name__ == "__main__":
    it = read_packings('N32~P1e-1~9000.txt')
    for i in it:
        print i
        break
