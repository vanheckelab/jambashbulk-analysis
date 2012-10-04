import os

import numpy as np

import ctypes
from ctypes import Structure, c_int, c_longdouble, cdll, byref

def to_longdouble(value, cutoff=4):
    import re
    part1 = value[:cutoff]; part2 = re.sub('[0-9]', '0', part1) + value[cutoff:]
    return np.longdouble(part1)+np.longdouble(part2)

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
    print type(N), type(P0), type(alpha), type(delta), type(L)
    packingdata = PACKINGPARAMS()

    N = np.int_([N])
    N_ptr = N.ctypes.data_as(ctypes.POINTER(c_int))
    P0 = np.longdouble([P0])
    P0_ptr = P0.ctypes.data_as(ctypes.POINTER(c_longdouble))

    x = np.longdouble(x.copy())
    x_ptr = x.ctypes.data_as(ctypes.POINTER(c_longdouble))
    y = np.longdouble(y.copy())
    y_ptr = y.ctypes.data_as(ctypes.POINTER(c_longdouble))
    r = np.longdouble(r.copy())
    r_ptr = r.ctypes.data_as(ctypes.POINTER(c_longdouble))

    alpha = np.longdouble([alpha])
    alpha_ptr = alpha.ctypes.data_as(ctypes.POINTER(c_longdouble))
    delta = np.longdouble([delta])
    delta_ptr = delta.ctypes.data_as(ctypes.POINTER(c_longdouble))
    L = np.longdouble([L])
    L_ptr = L.ctypes.data_as(ctypes.POINTER(c_longdouble))

    c_get_packing_data(N_ptr.contents, P0_ptr.contents,
                       x_ptr, y_ptr, r_ptr,
                       alpha_ptr.contents, delta_ptr.contents, L_ptr.contents,
                       byref(packingdata))

    return dict((field[0], getattr(packingdata, field[0])) for field in packingdata._fields_)

def convertLvectors(L1, L2):
    L = np.sqrt(L1[0] * L2[1])
    alpha = L2[0] / L
    delta = np.sqrt(L2[1] / L1[0]) - 1.0;

    return {'alpha': alpha, 'delta': delta, 'L': L}

def main():
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
    np.set_printoptions(precision=1)
    for lastdigits in ["75", "76", "77", "78", "79", "80", "81"]:
        L1 = np.longdouble([to_longdouble("9.4949297755650788") , "0.0000000000000000"])
        L2 = np.longdouble([to_longdouble("0.31310371412050"+lastdigits) , to_longdouble("9.4800157565608026")])

        x = np.longdouble(data[::3])
        y = np.longdouble(data[1::3])
        r = np.longdouble(data[2::3])
        #print "%s, using:" % os.path.split(__file__)[1]
        #print __version__
        packdata = get_packing_data(N, P0, x, y, r, **convertLvectors(L1, L2))
        #print packdata
        print lastdigits, packdata["gg"], 0.000109166-packdata["U"], packdata["sxy"]

if __name__=="__main__":
    main()
