/////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Convert geodetic coordinates from WGS84 to the HGRS87 (EGSA87/ΕΓΣΑ87) Greek geodetic system
//  and vice versa.
//
//  Copyright (C) 2021  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
/////////////////////////////////////////////////////////////////////////////////////////////////


#include <math.h>
#include <float.h> // DBL_MAX

#include "egsa87.h"

// philam contains the latitude & longitude in radians
// xy contains the EGSA87 coordinates upon return (in meters)
// 
// Based on wgs84egsa87.m from https://github.com/mlourakis/egsa87
void wgs84egsa87(double philam[2], double xy[2])
{
const double GE_WGS84_Alpha=6378137.000;
const double GE_WGS84_F_INV=298.257223563;

double phi, lambda, t, e2, rad, x, y, z, x2, y2, z2, phi2, lambda2, kappa0, lambda0, dl, n2;
double xoffset, yoffset, aradius, aradius_old, f, et2, e4, e6, e8, M0, M2, M4, M6, M8, Mi, Ni, L;
int acount;

  phi=philam[0];  
  lambda=philam[1];  

// wgs84_philambda_to_xy

//1. displace_geodetic_system: phi, lambda --> phi2, lambda2
  //1.1 philambda_to_xyz()
  t=(1.0 - 1.0/GE_WGS84_F_INV);
  e2=1.0 - t*t;
  rad=GE_WGS84_Alpha/sqrt(1.0-e2*sin(phi)*sin(phi));

	x=rad*cos(phi)*cos(lambda);
	y=rad*cos(phi)*sin(lambda);
	z=rad*(1.0-e2)*sin(phi);

  //1.2 displace_geocentric_system: x,y,z --> x2, y2, z2
  x2=x + 199.72; // GS_GR87_DX  +199,723 
  y2=y -  74.03; // GS_GR87_DY   -74,030
  z2=z - 246.02; // GS_GR87_DZ  -246,018

  //1.3 ellipsoid_xyz_to_philambda: x2, y2, z2 --> phi2, lambda2
  aradius=GE_WGS84_Alpha;
  // sphere_xyz_to_philambda
  if(fabs(z2)<aradius)
    phi2=asin(z2/aradius);
  else{
    if(z2>0)
      phi2=0.5*M_PI;
    else
      phi2=-0.5*M_PI;
  }

  if(fabs(x2)>0.001)
    lambda2=atan(y2/x2);
  else{
    if(y2>0)
      lambda2=0.5*M_PI;
    else
      lambda2=-0.5*M_PI;
  }

  if(x2<0)
    lambda2=M_PI - lambda2;

  f=1.0/1.0/GE_WGS84_F_INV;
  et2=e2/((1.0-f)*(1.0-f));
  phi2=atan( z2*(1.0+et2) / sqrt(x2*x2+y2*y2) );
  acount=0;
  aradius_old=DBL_MAX;
  while(fabs(aradius-aradius_old)>0.00005 && acount<100){
    acount=acount+1;
    aradius_old=aradius;
    aradius=GE_WGS84_Alpha/sqrt(1.0-e2*sin(phi2)*sin(phi2));// ellipsoid_main_normal_section_radius(phi2);
    phi2=atan( (z2+e2*aradius*sin(phi2)) / sqrt(x2*x2+y2*y2) );
  }


//2. project_philambda_to_xy: phi2, lambda2  -->  x, y, gamma
  // utm_project_philambda_to_xy
  kappa0=0.9996; // PC_GR87_KAPPA
  lambda0=(24.00*M_PI/180.00); // PC_GR87_LAMBDA0
  xoffset=500000; // PC_GR87_XOFFSET
  yoffset=0.00; // PC_GR87_YOFFSET
  dl=lambda2 - lambda0;
  t=tan(phi2);
  n2=e2*cos(phi2)*cos(phi2)/(1.0-e2);
  L=dl*cos(phi2);
  //Ni=ellipsoid_main_normal_section_radius(phi2);
  Ni=GE_WGS84_Alpha/sqrt(1.0-e2*sin(phi2)*sin(phi2));

  //Mi=ellipsoid_arc(phi2);
  e4 = e2 * e2;
  e6 = e4 * e2;
  e8 = e6 * e2;
  M0 = 1.0 + 0.75*e2 + 0.703125*e4 + 0.68359375*e6 + 0.67291259765625*e8;
  M2 = 0.375*e2 + 0.46875*e4 + 0.5126953125*e6 + 0.538330078125*e8;
  M4= 0.05859375*e4 + 0.1025390625*e6 + 0.25*e8;
  M6= 0.01139322916666667*e6 + 0.025634765625*e8;
  M8= 0.002408551504771226*e8;
  Mi=GE_WGS84_Alpha*(1.0-e2)*(M0*phi2-M2*sin(2*phi2)+M4*sin(4*phi2)-
    M6*sin(6*phi2)+M8*sin(8*phi2));

  x= (((5-18*t*t+t*t*t*t+14*n2-58*t*t*n2)*L*L/120.00 + (1.0-t*t+n2)/6.00)*L*L+1.0)*
       L*kappa0*Ni + xoffset;
  y= Mi + (Ni*t/2)*L*L + (Ni*t/24)*(5-t*t+9*n2+4*n2*n2)*L*L*L*L +
    (Ni*t/720)*(61-58*t*t)*L*L*L*L*L*L;
  y=y*kappa0 + yoffset;
  //gamma=t*L+t*L*L*L*(1.0+3*n2+4*n2*n2)/3+t*(2-t*t)*L*L*L*L*L/15;

  xy[0]=x;
  xy[1]=y;
}


// xy contains the EGSA87 coordinates in meters
// philam contains the latitude & longitude upon return (in radians)
//
// Based on egsa87wgs84.m from https://github.com/mlourakis/egsa87
void egsa87wgs84(double xy[2], double philam[2])
{
const double GE_WGS84_Alpha=6378137.000;
const double GE_WGS84_F_INV=298.257223563;

const double kappa0=0.9996; // PC_GR87_KAPPA
const double lambda0=(24.00*M_PI/180.00); // PC_GR87_LAMBDA0
const double xoffset=500000; // PC_GR87_XOFFSET
const double yoffset=0.00; // PC_GR87_YOFFSET

double x, y, t, phi, lambda, f, et2, l, e2, e4, e6, e8, M0, M2, M4, M6, M8, Mi, N0, P, P2;
double f0, f0_old, n2, phi2, lambda2, rad, z, x2, y2, z2, aradius, aradius_old;
int acount;

  x=xy[0];
  y=xy[1];

// xy_to_wgs84_philambda

//1. project_xy_to_philambda
// utm_project_xy_to_philambda: x, y --> phi2, lambda2, gamma

  x=x-xoffset;
  y=y-yoffset;
  t=(1.0 - 1.0/GE_WGS84_F_INV);
  e2=1.0 - t*t;
  f=1.0/1.0/GE_WGS84_F_INV;
  et2=e2/((1.0-f)*(1.0-f));

  // f0=fit(y)
  l=y/kappa0;
  f0=l/GE_WGS84_Alpha;
  f0_old=10*f0;
  acount=0;
  while(fabs(f0-f0_old)>1e-17 && acount<100){
    acount=acount+1;
    f0_old=f0;

    //Mi=ellipsoid_arc(f0);
    e4=e2 * e2;
    e6=e4 * e2;
    e8=e6 * e2;
    M0=1.0 + 0.75*e2 + 0.703125*e4 + 0.68359375*e6 + 0.67291259765625*e8;
    M2=0.375*e2 + 0.46875*e4 + 0.5126953125*e6 + 0.538330078125*e8;
    M4= 0.05859375*e4 + 0.1025390625*e6 + 0.25*e8;
    M6= 0.01139322916666667*e6 + 0.025634765625*e8;
    M8= 0.002408551504771226*e8;
    Mi=GE_WGS84_Alpha*(1.0-e2)*(M0*f0-M2*sin(2*f0)+M4*sin(4*f0)-
      M6*sin(6*f0)+M8*sin(8*f0));

    f0=f0 + (l - Mi)/(GE_WGS84_Alpha/sqrt(1.0-e2*sin(f0)*sin(f0))); // ellipsoid_main_normal_section_radius(f0);
  }

  N0=GE_WGS84_Alpha/sqrt(1.0-e2*sin(f0)*sin(f0)); // ellipsoid_main_normal_section_radius(f0);
  t=tan(f0);
  n2=et2*cos(f0)*cos(f0);
  P=(x/(kappa0*N0));
  P2=P*P;

  phi2=(((-(61+90*t*t+45*t*t*t*t)*P2/720 +
     (5+3*t*t+6*n2-3*n2*n2-6*t*t*n2-9*t*t*n2*n2)/24)*P2-
     (1.0+n2)/2)*P2)*t + f0;

  lambda2=((((5+6*n2+28*t*t+8*t*t*n2+24*t*t*t*t)*P2/120 -
    (1.0+2*t*t+n2)/6)*P2+1.0)*P)/cos(f0) + lambda0;

  //gamma=x*t/(kappa0*N0)-t*(1.0+t*t-n2-2*n2*n2)*x*x*x/(3*kappa0*kappa0*kappa0*N0*N0*N0)+t*(2+5*t*t+3*t*t*t*t)*pow(x/(kappa0*N0),5)/15;


//2. displace_geodetic_system: phi2, lambda2 --> phi, lambda
  //2.1 philambda_to_xyz()
  rad=GE_WGS84_Alpha/sqrt(1.0-e2*sin(phi2)*sin(phi2));

	x=rad*cos(phi2)*cos(lambda2);
	y=rad*cos(phi2)*sin(lambda2);
	z=rad*(1.0-e2)*sin(phi2);

  //2.2 displace_geocentric_system: x,y,z --> x2, y2, z2
  x2=x - 199.72; // GS_GR87_DX  +199,723 
  y2=y +  74.03; // GS_GR87_DY   -74,030
  z2=z + 246.02; // GS_GR87_DZ  -246,018

  //2.3 ellipsoid_xyz_to_philambda: x2, y2, z2 --> phi, lambda
  aradius=GE_WGS84_Alpha;
  // sphere_xyz_to_philambda
  if(fabs(z2)<aradius)
    phi=asin(z2/aradius);
  else{
    if(z2>0)
      phi=0.5*M_PI;
    else
      phi=-0.5*M_PI;
  }

  if(fabs(x2)>0.001)
    lambda=atan(y2/x2);
  else{
    if(y2>0)
      lambda=0.5*M_PI;
    else
      lambda=-0.5*M_PI;
  }

  if(x2<0)
    lambda=M_PI - lambda;

  phi=atan( z2*(1.0+et2) / sqrt(x2*x2+y2*y2) );
  acount=0;
  aradius_old=DBL_MAX;
  while(fabs(aradius-aradius_old)>0.00005 && acount<100){
    acount=acount+1;
    aradius_old=aradius;
    aradius=GE_WGS84_Alpha/sqrt(1-e2*sin(phi)*sin(phi)); // ellipsoid_main_normal_section_radius(phi);
    phi=atan( (z2+e2*aradius*sin(phi)) / sqrt(x2*x2+y2*y2) );
  }

  philam[0]=phi;
  philam[1]=lambda;
}
