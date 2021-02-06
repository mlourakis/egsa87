% wgs84egsa87 - Convert geodetic coordinates from WGS84 to the HGRS87 (EGSA87/ΕΓΣΑ87) Greek geodetic system
%
%    [x y] = wgs84egsa87(phi, lambda)
%
%    phi        - latitude in radians
%    lambda     - longitude in radians
%
% Returns
%    x, y       - EGSA87 coordinates (meters)

% see https://en.wikipedia.org/wiki/Hellenic_Geodetic_Reference_System_1987
% Based on icoordstrans -- https://github.com/skozan/icoordstrans/
%
% Manolis Lourakis 2020
% Institute of Computer Science, Foundation for Research & Technology - Hellas
% Heraklion, Crete, Greece

% Feb  2020  - Initial version. (v. 1.0)
% Feb  2021  - Minor update. (v. 1.1)


function [x y] = wgs84egsa87(phi, lambda)
  GE_WGS84_Alpha=6378137.000;
  GE_WGS84_F_INV=298.257223563;

% wgs84_philambda_to_xy

%1. displace_geodetic_system: phi, lambda --> phi2, lambda2
  %1.1 philambda_to_xyz()
  e2=1 - (1 - 1/GE_WGS84_F_INV)^2;
  rad=GE_WGS84_Alpha/sqrt(1-e2*sin(phi)*sin(phi));

	x=rad*cos(phi)*cos(lambda);
	y=rad*cos(phi)*sin(lambda);
	z=rad*(1-e2)*sin(phi);

  %1.2 displace_geocentric_system: x,y,z --> x2, y2, z2
  x2=x + 199.72; % GS_GR87_DX  +199,723 
  y2=y -  74.03; % GS_GR87_DY   -74,030
  z2=z - 246.02; % GS_GR87_DZ  -246,018

  %1.3 ellipsoid_xyz_to_philambda: x2, y2, z2 --> phi2, lambda2
  aradius=GE_WGS84_Alpha;
  % sphere_xyz_to_philambda
  if(abs(z2)<aradius)
    phi2=asin(z2/aradius);
  else
    if(z2>0)
      phi2=0.5*pi;
    else
      phi2=-0.5*pi;
    end
  end

  if(abs(x2)>0.001)
    lambda2=atan(y2/x2);
  else
    if(y2>0)
      lambda2=0.5*pi;
    else
      lambda2=-0.5*pi;
    end
  end

  if(x2<0)
    lambda2=pi - lambda2;
  end

  f=1/1/GE_WGS84_F_INV;
  et2=e2/((1-f)*(1-f));
  phi2=atan( z2*(1+et2) / sqrt(x2*x2+y2*y2) );
  acount=0;
  aradius_old=Inf;
  while(abs(aradius-aradius_old)>0.00005 && acount<100)
    acount=acount+1;
    aradius_old=aradius;
    aradius=GE_WGS84_Alpha/sqrt(1-e2*sin(phi2)*sin(phi2)); % ellipsoid_main_normal_section_radius(phi2);
    phi2=atan( (z2+e2*aradius*sin(phi2)) / sqrt(x2*x2+y2*y2) );
  end


%2. project_philambda_to_xy: phi2, lambda2  -->  x, y, gamma
  % utm_project_philambda_to_xy
  kappa0=0.9996; % PC_GR87_KAPPA
  lambda0=(24.00*pi/180.00); % PC_GR87_LAMBDA0
  xoffset=500000; % PC_GR87_XOFFSET
  yoffset=0.00; % PC_GR87_YOFFSET
  dl=lambda2 - lambda0;
  t=tan(phi2);
  n2=e2*cos(phi2)*cos(phi2)/(1-e2);
  L=dl*cos(phi2);
  %Ni=ellipsoid_main_normal_section_radius(phi2);
  Ni=GE_WGS84_Alpha/sqrt(1-e2*sin(phi2)*sin(phi2));

  %Mi=ellipsoid_arc(phi2);
  e4 = e2 * e2;
  e6 = e4 * e2;
  e8 = e6 * e2;
  M0 = 1 + 0.75*e2 + 0.703125*e4 + 0.68359375*e6 + 0.67291259765625*e8;
  M2 = 0.375*e2 + 0.46875*e4 + 0.5126953125*e6 + 0.538330078125*e8;
  M4= 0.05859375*e4 + 0.1025390625*e6 + 0.25*e8;
  M6= 0.01139322916666667*e6 + 0.025634765625*e8;
  M8= 0.002408551504771226*e8;
  Mi=GE_WGS84_Alpha*(1-e2)*(M0*phi2-M2*sin(2*phi2)+M4*sin(4*phi2)-...
    M6*sin(6*phi2)+M8*sin(8*phi2));

  x= (((5-18*t*t+t*t*t*t+14*n2-58*t*t*n2)*L*L/120.00 + (1-t*t+n2)/6.00)*L*L+1)*...
       L*kappa0*Ni + xoffset;
  y= Mi + (Ni*t/2)*L*L + (Ni*t/24)*(5-t*t+9*n2+4*n2*n2)*L*L*L*L +...
    (Ni*t/720)*(61-58*t*t)*L*L*L*L*L*L;
  y=y*kappa0 + yoffset;
  %gamma=t*L+t*L*L*L*(1+3*n2+4*n2*n2)/3+t*(2-t*t)*L*L*L*L*L/15;

end
