'''
 egsa87wgs84 - Convert geodetic coordinates from the HGRS87 (EGSA87/ΕΓΣΑ87) Greek geodetic system to WGS84

 egsa87wgs84(x, y, z) 

    x, y        - EGSA87 coordinates (meters)

 Returns
    phi, lambda - latitude & longitude in radians

 see https://en.wikipedia.org/wiki/Hellenic_Geodetic_Reference_System_1987
 Based on icoordstrans -- https://github.com/skozan/icoordstrans/

 Manolis Lourakis 2020
 Institute of Computer Science, Foundation for Research & Technology - Hellas
 Heraklion, Crete, Greece

 Dimitrios "mitselec" Dimopoulos - Python version

 Feb  2020  - Initial version. (v. 1.0)
 Feb  2021  - Minor update. (v. 1.1)
 July 2021  - Python Version
 '''

import math

def egsa87wgs84(x, y):
    GE_WGS84_Alpha=6378137.000
    GE_WGS84_F_INV=298.257223563

    kappa0=0.9996                    # PC_GR87_KAPPA
    lambda0=(24.00*math.pi/180.00)   # PC_GR87_LAMBDA0
    xoffset=500000                   # PC_GR87_XOFFSET
    yoffset=0.00                     # PC_GR87_YOFFSET

    # xy_to_wgs84_philambda

    #1. project_xy_to_philambda
    # utm_project_xy_to_philambda: x, y --> phi2, lambda2, gamma

    x=x-xoffset
    y=y-yoffset
    e2=1 - (1 - 1/GE_WGS84_F_INV)**2
    f=1/1/GE_WGS84_F_INV
    et2=e2/((1-f)*(1-f))

    # f0=fit(y)
    l=y/kappa0
    f0=l/GE_WGS84_Alpha
    f0_old=10*f0
    acount=0

    while(abs(f0-f0_old)>10**(-17) and acount<100):
        acount=acount+1
        f0_old=f0

        #Mi=ellipsoid_arc(f0);
        e4=e2 * e2
        e6=e4 * e2
        e8=e6 * e2
        M0=1 + 0.75*e2 + 0.703125*e4 + 0.68359375*e6 + 0.67291259765625*e8
        M2=0.375*e2 + 0.46875*e4 + 0.5126953125*e6 + 0.538330078125*e8
        M4= 0.05859375*e4 + 0.1025390625*e6 + 0.25*e8
        M6= 0.01139322916666667*e6 + 0.025634765625*e8
        M8= 0.002408551504771226*e8
        Mi=GE_WGS84_Alpha*(1-e2)*(M0*f0-M2*math.sin(2*f0)+M4*math.sin(4*f0)- M6*math.sin(6*f0)+M8*math.sin(8*f0))

        f0=f0 + (l - Mi)/(GE_WGS84_Alpha/math.sqrt(1-e2*math.sin(f0)*math.sin(f0)))     # ellipsoid_main_normal_section_radius(f0)


    N0=GE_WGS84_Alpha/math.sqrt(1-e2*math.sin(f0)*math.sin(f0))                         # ellipsoid_main_normal_section_radius(f0)
    t=math.tan(f0)
    n2=et2*math.cos(f0)*math.cos(f0)
    P=(x/(kappa0*N0))
    P2=P*P

    phi2=(((-(61+90*t*t+45*t*t*t*t)*P2/720 +\
    (5+3*t*t+6*n2-3*n2*n2-6*t*t*n2-9*t*t*n2*n2)/24)*P2-\
    (1+n2)/2)*P2)*t + f0

    lambda2=((((5+6*n2+28*t*t+8*t*t*n2+24*t*t*t*t)*P2/120 -\
    (1+2*t*t+n2)/6)*P2+1)*P)/math.cos(f0) + lambda0

    #2. displace_geodetic_system: phi2, lambda2 --> phi, lambda
    #2.1 philambda_to_xyz()
    rad=GE_WGS84_Alpha/math.sqrt(1-e2*math.sin(phi2)*math.sin(phi2))

    x=rad*math.cos(phi2)*math.cos(lambda2)
    y=rad*math.cos(phi2)*math.sin(lambda2)
    z=rad*(1-e2)*math.sin(phi2)

    #2.2 displace_geocentric_system: x,y,z --> x2, y2, z2
    x2=x - 199.72 # GS_GR87_DX  +199,723 
    y2=y +  74.03 # GS_GR87_DY   -74,030
    z2=z + 246.02 # GS_GR87_DZ  -246,018

    #2.3 ellipsoid_xyz_to_philambda: x2, y2, z2 --> phi, lambda
    aradius=GE_WGS84_Alpha;
    # sphere_xyz_to_philambda
    if(abs(z2)<aradius):
        phi=math.asin(z2/aradius)
    else:
        if(z2>0):
            phi=0.5*math.pi
        else:
            phi=-0.5*math.pi

    if(abs(x2)>0.001):
        lambda1=math.atan(y2/x2)
    else:
        if(y2>0):
            lambda1=0.5*math.pi
        else:
            lambda1=-0.5*math.pi

    if(x2<0):
        lambda1=math.pi - lambda1

    phi=math.atan( z2*(1+et2) / math.sqrt(x2*x2+y2*y2) )
    acount=0
    aradius_old=10**38
    while(abs(aradius-aradius_old)>0.00005 and acount<100):
        acount+=1
        aradius_old=aradius
        aradius=GE_WGS84_Alpha/math.sqrt(1-e2*math.sin(phi)*math.sin(phi)) # ellipsoid_main_normal_section_radius(phi);
        phi=math.atan( (z2+e2*aradius*math.sin(phi)) / math.sqrt(x2*x2+y2*y2) )
        
    return(phi,lambda1)

  
