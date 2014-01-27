import os

import numpy as np

import ctypes
from ctypes import Structure, c_int, c_longdouble, cdll, byref

class PACKINGPARAMS(Structure):
    _fields_ = [('P', c_longdouble),
                ('phi', c_longdouble),
                ('Z', c_longdouble),
                ('Ncorrected', c_int),
                ('sxx', c_longdouble),
                ('sxy', c_longdouble),
                ('syy', c_longdouble),
                ('U', c_longdouble),
                ('H', c_longdouble),
                ('gg', c_longdouble)]

dll_path = os.path.join(
        os.path.split(os.path.abspath(__file__))[0],
        'bin',
        'jamBashbulk.so')
dll = cdll.LoadLibrary(dll_path)

__version__ = ctypes.c_char_p.in_dll(dll, "FILE_HEADER_C").value

c_get_packing_data = dll.get_packing_data

def get_packing_data(N, P0, x, y, r, alpha, delta, L):
    packingdata = PACKINGPARAMS()

    N = c_int(N)
    P0 = c_longdouble(P0)

    x = np.longdouble(x)
    x_ptr = x.ctypes.data_as(ctypes.POINTER(c_longdouble))
    y = np.longdouble(y)
    y_ptr = y.ctypes.data_as(ctypes.POINTER(c_longdouble))
    r = np.longdouble(r)
    r_ptr = r.ctypes.data_as(ctypes.POINTER(c_longdouble))

    alpha = c_longdouble(alpha)
    delta = c_longdouble(delta)
    L = c_longdouble(L)

    c_get_packing_data(N, P0,
                       x_ptr, y_ptr, r_ptr,
                       alpha, delta, L,
                       byref(packingdata))

    return dict((field[0], getattr(packingdata, field[0])) for field in packingdata._fields_)

def convertLvectors(L1, L2):
    L = np.sqrt(L1[0] * L2[1])
    alpha = L2[0] / L
    delta = np.sqrt(L2[1] / L1[0]) - 1.0;

    return {'alpha': alpha, 'delta': delta, 'L': L}

def updatePackingParams(packing):
    x = np.longdouble([row["x"] for row in packing["particles"]]) + np.longdouble([row["x_err"] for row in packing["particles"]])
    y = np.longdouble([row["y"] for row in packing["particles"]]) + np.longdouble([row["y_err"] for row in packing["particles"]])
    r = np.longdouble([row["r"] for row in packing["particles"]])

    pd = get_packing_data(
        packing["N"],
        packing["P0"],
        x, y, r,
        **convertLvectors(packing["L1"], packing["L2"]))
    # {'phi': 0.8011895053312124,
    #   'Ncorrected': 21,
    #   'gg': 1.818953279284722e-15,
    #   'P': 0.00010000000000000351,
    #   'U': 1.4889779929429936e-06,
    #   'H': 0.012376671564328283,
    #   'sxx': -0.00010000000000007477,
    #   'sxy': 5.570205614118088e-18,
    #   'syy': -9.999999999993227e-05,
    #   'Z': 4.095238095238095
    #}
    packing.update(pd)
    packing.update(convertLvectors(packing["L1"], packing["L2"]))
    packing["N - Ncorrected"] = packing["N"] - packing["Ncorrected"]
    packing["runtime (s)"] = np.nan
    packing["current time"] = 'unknown'
    packing["FIRE iteration count"] = -1
    packing["CG iteration count"] = -1
    packing["dH"] = np.nan
    packing["dU"] = np.nan
    packing["maxGrad"] = packing["gg"]
    packing["Uhelper"] = packing["U"]
    packing["Z_calc"] = packing["Z"]

    return packing

if __name__ == '__main__':
    N = 16
    P0 = 0.001
    data = [
    8.3954120736007511 ,    5.9713347462851643 ,    1.0000000000000000 ,
    7.7732133087976493 ,    7.8703131844924932 ,    1.0000000000000000 ,
    6.3809728475044440 ,    9.3027880393547451 ,    1.0000000000000000 ,
    2.8009430583128803 ,    8.6790180704999307 ,    1.0000000000000000 ,
    7.8966028959990910 ,    0.6292000854999064 ,    1.0000000000000000 ,
    4.7338725981321079 ,    8.1703371898420650 ,    1.0000000000000000 ,
    2.3103976904041286 ,    2.7410961329297243 ,    1.0000000000000000 ,
    6.2834832183083400 ,    1.8104815211023599 ,    1.0000000000000000 ,
    0.7857890932593265 ,    0.8893224628442692 ,    1.3999999999999999 ,
    7.8349418788527788 ,    3.6388967215266238 ,    1.3999999999999999 ,
    4.4582447016779932 ,    3.8105760378512087 ,    1.3999999999999999 ,
    0.9302650466673738 ,    4.6991730092283172 ,    1.3999999999999999 ,
    6.0060829446461095 ,    6.1396900822624010 ,    1.3999999999999999 ,
    0.6616435186167540 ,    7.5983195833701074 ,    1.3999999999999999 ,
    3.2124952245017365 ,    6.3164065410892209 ,    1.3999999999999999 ,
    4.0087550307426671 ,    1.0510396876630096 ,    1.3999999999999999 ,
    ]

    L1 = np.array([9.4949297755650788 , 0.0000000000000000])
    L2 = np.array([0.3131037141205079 , 9.4800157565608026])

    x = np.longdouble(data[::3])
    y = np.longdouble(data[1::3])
    r = np.longdouble(data[2::3])
    print "%s, using:" % os.path.split(__file__)[1]
    print __version__
    print get_packing_data(N, P0, x, y, r, **convertLvectors(L1, L2))
