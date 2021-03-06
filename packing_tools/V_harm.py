# -*- coding: utf-8 -*-
"""
Tools for working with packings:
 - get_contacts to determine overlaps and the connection matrix,
 - V_harm to calculate energy, enthalpy and internal pressure,
 - grad to calculate the gradient and stresses,
 - Hess to calculate the extended Hessian,
 - LinDef and El_Con to determine elastic moduli.

Set V_harm.VISCOUS=True to disable rattler detection (for use in simulations with 
viscous interactions, where the contept 'rattler' doesn't make sense)

"""
import numpy as np
from numpy import floor, sqrt, sum, tril, where, diag, array, int_

try:
    from numpy import float128
except ImportError:
    print "float64 instead of float128..."
    from numpy import float64 as float128

from numpy import linalg as LA

from copy import copy

VISCOUS = False

def get_contacts(packing):
    """ In: packing={'particles': {'x': [...], 'y': [...], 'r': [...]},
                    {'L1': ..., 'L2': ..., 'L': ...}}

        Out: {xij, yij, rij,
              dij (= overlaps, no underlaps),
              dijfull (= both, underlaps as negative values),
              nx, ny (= periodic boundary crossings for particle pair i,j),
              rattlers (= list of all rattlers),
              connmatrix (= Cij, boolean array whether particle i and j are in contact)
              }
    """
    particles = packing['particles']
    
    x=float128(particles['x'])
    y=float128(particles['y'])
    r=float128(particles['r'])
 
    try:
        x += float128(particles['x_err'])
        y += float128(particles['y_err'])
    except ValueError:
        pass
    except KeyError:
        pass
    
    #xi, xj = np.meshgrid(x,x)
    #yi, yj = np.meshgrid(y,y)
    #Ri, Rj = np.meshgrid(r,r)
    xi = x; xj = x[:,np.newaxis]
    yi = y; yj = y[:,np.newaxis]
    Ri = r; Rj = r[:,np.newaxis]
    Rij=Ri+Rj
    
    lxx, lxy = packing["L1"]
    lyx, lyy = packing["L2"]
    
    L = sqrt(lxx * lyy)

    xij = xj-xi
    yij = yj-yi
    
    ny = -int_(floor((yij + 0.5 * lyy)/lyy))
    
    xij += ny * lyx
    yij += ny * lyy
    
    nx = -int_(floor((xij + 0.5 * lxx)/lxx))

    xij += nx * lxx
    yij += nx * lxy
    
    rij = sqrt(xij**2 + yij**2)
    dij = Rij - rij
    connmatrix = dij > 0
    np.fill_diagonal(connmatrix, False)
    
    if not VISCOUS:
        while True:
            Ncontacts = sum(connmatrix, axis=0)
            rattlers = where((Ncontacts > 0) * (Ncontacts <= 2))[0]                
            if len(rattlers) == 0:
                break
            connmatrix[rattlers,:] = False
            connmatrix[:, rattlers] = False
    
    
    Ncontacts = sum(connmatrix, axis=0)
    
    if VISCOUS:
        rattlers = np.array([], dtype=np.int)
    else:
        rattlers=where((Ncontacts <= 2))[0]
    dijfull = dij.copy()
    dij[~connmatrix] = 0
    return {'xij': xij,'yij': yij,'rij':rij,'dij': dij, 'dijfull': dijfull, 'nx':nx,'ny':ny, 'connmatrix': connmatrix, 'rattlers':rattlers} 

def V_harm(conts,packing):
    """
    Calculates energy, pressure and enthalpy.
    In: conts = get_contacts(packing)['connmatrix'],
        packing = packing
    
    Out: {'U':U,'P':P,'H':H}
    """
    #the elastic constant :
    k_el=1
    
    # packing global data : 
    P0=packing['P0']
    lxx, lxy = packing["L1"]
    lyx, lyy = packing["L2"]
    
    L = sqrt(lxx * lyy)
    
    # distancs matrices
    xij = conts['xij']
    yij = conts['yij']
    rij = conts['rij']
    dij = conts['dij']
    
    # energy :                 
    U=0.5*sum(0.5*k_el*dij**2)
    
    # Pressure : 
    P=0.5*0.5*sum(k_el*dij*rij)/L**2
    
    H=U+P0*L**2

    return {'U':U,'P':P,'H':H}        
        
        
def grad(conts,packing):
    """
    Calculate gradients and stresses.
    """
    #the elastic constant :
    k_el=1
    
    # packing global data : 
    P0=packing['P0']
    particles=packing['particles']
 
    lxx, lxy = packing["L1"]
    lyx, lyy = packing["L2"]
    
    L = sqrt(lxx * lyy)
    r=particles['r']
    
    
    # distancs matrices
    xij = conts['xij']
    yij = conts['yij']
    rij = conts['rij']
    dij = conts['dij']
    
    # contact informations
    con=conts['connmatrix']
    nx=conts['nx']
    ny=conts['ny']
    xhatij=xij/(rij+diag(r))
    yhatij=yij/(rij+diag(r))
    
    # Gradient
    gradx=0.5*(-sum(k_el*dij*xhatij,0))
    
    
    grady=0.5*(-sum(k_el*dij*yhatij,0))
    # Gradient alpha
    g_alpha= -0.5*sum(k_el*dij*ny*xhatij)
    # Gradient beta
    g_delta=np.nan # STUK -0.5*sum(k_el*dij*(ny*yhatij-nx*xhatij/(1+delta)**2))
    # Gradient L
    P=0.5*0.5*sum(k_el*dij*rij)/L**2
    g_L=2*L*(P-P0)
    
    
    
    # Stress and Fabric tensors :
    #        Stress :        
    sxx=0.5 * sum(k_el * dij * xhatij * xij)/L**2
    sxy=0.5 * sum(k_el * dij * xhatij * yij)/L**2
    syy=0.5 * sum(k_el * dij * yhatij * yij)/L**2
    
    Sig=array([sxx,sxy,sxy,syy])
    Sig=Sig.reshape((2,2))
    
    
    #        fabric :        
    rxx=0.5 * sum(xhatij * xhatij)/L**2
    rxy=0.5 * sum(xhatij * yhatij)/L**2
    ryy=0.5 * sum(yhatij*yhatij)/L**2
    
    Rho=array([sxx,sxy,sxy,syy])
    Rho=Rho.reshape((2,2))
    
    
    return {'gx':gradx, 'gy':grady,'g_alpha':g_alpha,'g_delta':g_delta,'g_L':g_L, 'Sig' : Sig, 'Rho':Rho}
        
        
def Hess(conts,packing):
    """
    Calculate extended Hessian.
      
    In: conts = get_contacts(packing)['connmatrix'],
        packing = packing  

    Out: extended Hessian (N+4 by N+4 array)
    """
    #the elastic constant :
    k_el=1
    
    # packing global data : 
    P0=packing['P0']
    particles=packing['particles']

    lxx, lxy = packing["L1"]
    lyx, lyy = packing["L2"]
    
    L = sqrt(lxx * lyy)
    r=particles['r']
    # distancs matrices
    xij = conts['xij']
    yij = conts['yij']
    rij = conts['rij']
    dij = conts['dij']
    
    # contact informations
    con = conts['connmatrix']
    nx = conts['nx']
    ny = conts['ny']
    xhatij = con * xij/(rij+diag(r))
    yhatij = con * yij/(rij+diag(r))
    
    
    d_r = dij/(rij+diag(r))
    
    # Hessian : regular core block, with only the particle position as variables
    
    KXX = d_r * yhatij**2 - xhatij**2
    KXX = k_el*KXX
    
    KYY = d_r * xhatij**2 - yhatij**2
    KYY = k_el*KYY
    
    KXY = -xhatij * yhatij * (1 + d_r)
    KXY = k_el*KXY
    
    # Diagonal needs to be filled, but this will be done later...
    
    
    # interactions with boundaries : 
    
    kb_xa = np.hstack((sum(nx * KXX,0),sum(nx * KXY,0)))
    kb_xb = np.hstack((sum(nx * KXY,0),sum(nx * KYY,0)))
    kb_xc = np.hstack((sum(ny * KXX,0),sum(ny * KXY,0)))
    kb_xd = np.hstack((sum(ny * KXY,0),sum(ny * KYY,0)))
    
    Kbxy = np.vstack((kb_xa,kb_xb,kb_xc,kb_xd))
    
    #print Kbxy
    
    # FIND THE RIGHT SCRIPTURE TO ASSEMPLE ARRAYS
    a = [-sum(KXX * nx**2)     ,                     0,                0,                0]
    b=  [  -sum(KXY * nx**2)   ,-sum(KYY * nx**2)     ,                0,                0]
    c=  [  -sum(KXX * nx*ny)   ,-sum(KXY * nx * ny)-P0,-sum(KXX * ny**2),                0]
    d=  [  -sum(KXY * nx*ny)+P0,-sum(KYY * nx * ny)   ,-sum(KXY * ny**2),-sum(KYY * ny**2)]
    
    
    kb_loc=np.array([a,b,c,d])
    #kb.reshape((4,4))        
    
    kb=(kb_loc+kb_loc.transpose()-diag(diag(kb_loc)))/2
    
    # recomposition : 
    
    KXX=KXX-diag(sum(KXX,0))
    KYY=KYY-diag(sum(KYY,0))
    KXY=KXY-diag(sum(KXY,0))
    
    # K0=[KXX,KXY ; KXY,KYY]
    
    K0=np.vstack((np.hstack((KXX,KXY)),np.hstack((KXY,KYY))))
    K1=np.vstack((K0 , Kbxy))
    K2=np.vstack((Kbxy.transpose() , kb))
                    

    K = np.hstack((K1 , K2))        
    
    return K
    
    
            # TO BE CHECKED : 
    # not sure of the formats of all those length Lxy...
    # 
    # CAREFUL : with K a float128, should use linear algebra from SCIPY (and not numpy)
    # otherwise, the writing is the same :
    # from scipy import linalg as LA
    # [lam, vv]=LA.eigh(K)
    #
    # CAREFUL : not sure the h5 files are in float128 
        


# LinDef : returns the energy difference for a given strain (a,b,d) = ((a,d),(d,b)).
# It works also for a series of strains (save the diagonalization calculation time).
class PackingException(Exception):
    pass

def LinDef(K,rat,strain,packing):
    """
    Calculate linear elastic response to "strain" deformation.
    K = extended hessian = Hess(get_contacts(packing)['connmatrix'], packing),
    rat = list of rattlers = get_contacts(packing)['rattlers'],
    strain = deformation triplet [Δxx, Δyy, Δxy=Δyx]
    packing = packing

    returns: energy change ΔU
    """
    dU=np.array([])
    sz=strain.shape
    
    L1=packing['L1']
    L2=packing['L2']
    
    lxx, lxy = packing["L1"]
    lyx, lyy = packing["L2"]
    
    L = sqrt(lxx * lyy)
    
    if sz[1]!=3 :
        raise Exception("wrong strain input")
    
    N = packing['particles'].shape[0]
    Nr= N-rat.size
    all_rat=np.append(rat,rat+N)

    K= np.delete(K,all_rat,0)
    K= np.delete(K,all_rat,1)
    
    sK=K.shape[0]

    K0=K[0:sK-4,0:sK-4]
    if K0.shape == (0,0):
        raise PackingException("All particles are rattlers!")
    [lam,V]=LA.eigh(K0)
   
    for line in strain:
        Str_Ts=np.array([line[0],line[2],line[2],line[1]]).reshape(2,2)
        defor=np.array([Str_Ts.dot(L1),Str_Ts.dot(L2)]).reshape(4,1)

        def_ket=np.append(np.zeros((1,2*Nr),dtype=float128),defor)
        F_ket=K.dot(def_ket)
        
        F0_ket=-F_ket[0:sK-4]

        # $ K^{-1} = Q \Lambda Q^{-1} $
        #   Q = [v1, v2, v3, v4] (columns), Q^{-1} = Q.T
        #
        # so: u = K^{-1} F = Q \Lambda^{-1} Q^{-1} F
        # where       proj =   \Lambda^{-1} Q^{-1} F
        #           resp_0 = Q \Lambda^{-1} Q^{-1} F
        # with the detail that we set 1/λ = 0 if λ = 0 - a zero eigenvalue
        # means a zero-energy mode (not an infinite-energy mode)
        
        proj = (V.transpose().dot(F0_ket))/lam
        proj[0:2]=0
        
        resp_0=np.append(V.dot(proj),np.zeros((1,4),dtype=float128))

        response = def_ket+resp_0

        Edif=0.5*response.transpose().dot(K).dot(response)
        
        dU=np.append(dU,Edif)

    return dU
        
# El_con : calls lin-def for some given strains and deduces the elastic constants.

def El_Con(K,rattlers,packing):
    """ Calculate linear elastic moduli c1..c6.
    K = extended hessian = Hess(get_contacts(packing)['connmatrix'], packing),
    rat = list of rattlers = get_contacts(packing)['rattlers'],
    packing = packing.clear

    out = {'c': [c1, c2, c3, c4, c5, c6],
           'Gdc':Gdc,'Gac':Gac,
           'Udc':Udc,'Uac':Uac,
           'Dac':Dac}
    """
    lxx, lxy = packing["L1"]
    lyx, lyy = packing["L2"]
    
    L = sqrt(lxx * lyy)
    strain =np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1]])

    dU=LinDef(K,rattlers,strain,packing)
    c1=2*dU[0]/L**2
    c2=2*dU[1]/L**2
    c3=dU[2]/2/L**2
    c4=(dU[3]-dU[1]-dU[0])/L**2
    c5=(dU[4]-dU[0]-dU[2])/2/L**2
    c6=(dU[5]-dU[1]-dU[2])/2/L**2
    
    
    Gdc = 0.5*(c3 + (c1+c2-2*c4)/4)
    
    Gac = (c3 - (c1 + c2 - 2*c4)/4)**2 + (c5-c6)**2
    
    Gac=sqrt(Gac/8)
    #
    #print Gdc-sqrt(2)*Gac
    
    
    B = 0.25*(c1 + c2 + 2*c4)
    
    Udc = B + Gdc
    
    
    A2=sqrt((c1-c2)**2 + 4*(c5+c6)**2 )/2
    A4=sqrt(((c5-c6)**2 + (c3-(c1+c2-2*c4)/4)**2)/4)
    
    Uac = (A2**2+A4**2)/2
    Uac=sqrt(Uac)
    
    Dac=A2**2+4*A4**2
    Dac=sqrt(Dac/8)
    
    
    c=np.array([c1,c2,c3,c4,c5,c6],dtype=float128)
    
    return {'c':c,'Gdc':Gdc,'Gac':Gac,'Udc':Udc,'Uac':Uac,'Dac':Dac}
        
        
def Pmod(el_C):
    """Print elastic moduli as returned by El_Con"""
    c=el_C['c']
    print '\n'+'lxxxx : '+str(c[0])
    print 'lyyyy : '+str(c[1])
    print 'lxyxy : '+str(c[2])
    print 'lxxyy : '+str(c[3])
    print 'lxxxy : '+str(c[4])
    print 'lyyxy : '+str(c[5])+'\n'


    print 'Gdc : '+str(el_C['Gdc'])
    print 'Gac : '+str(el_C['Gac'])
    print 'Udc : '+str(el_C['Udc'])
    print 'Uac : '+str(el_C['Uac'])
    print 'Dac : '+str(el_C['Dac'])+'\n'
    