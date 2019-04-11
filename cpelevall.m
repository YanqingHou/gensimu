function [elev] = cpelevall(xsat,ysat,zsat,plhstation,xyzstation,lp);
%CPAZIELE: Compute elevation/azimuth of satellites/pseudolites
%
% This routine computes satellite positions in an earth-fixed WGS'84
% frame. If station coordinates are specified (in WGS'84), elevations and
% azimuths are computed as well.
%
% Syntax:
%    [elev] = cpaziele1 (xsat,ysat,zsat,plhstation,xyzstation)
%
% Input arguments:
%    xsat       - xsat(i) = X-coordinate (WGS'84) for satellite (i)
%    ysat       - ysat(i) = Y-coordinate (WGS'84) for satellite (i)
%    zsat       - zsat(i) = Z-coordinate (WGS'84) for satellite (i)
%    plhstation - WGS'84 coordinates for station (phi,lambda,height)
%    xyzstation - WGS'84 coordinates for station (x,y,z)
%
% Output arguments:
%    elev       - elev(i,j) = elevation for satellite (i)
% 

% ----------------------------------------------------------------------
% File.....: cpelev.m
% Date.....: 25-MAY-1999
% Version..: 1.0sv
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% Edited by: Sandra Verhagen
%            28-jun-2000

% -------------------------------------------------
% --- Compute azimuths/elevations, if requested ---
% -------------------------------------------------
for p = 1:lp
   
xsat = xsat - xyzstation(p,1);
ysat = ysat - xyzstation(p,2);
zsat = zsat - xyzstation(p,3);
 
rotmat(1,1) = - sin(plhstation(p,1)) * cos(plhstation(p,2));
rotmat(1,2) = - sin(plhstation(p,1)) * sin(plhstation(p,2));
rotmat(1,3) =   cos(plhstation(p,1))                     ;
rotmat(2,1) = -                      sin(plhstation(p,2));
rotmat(2,2) =                        cos(plhstation(p,2));
rotmat(2,3) =   0d0;
rotmat(3,1) = - cos(plhstation(p,1)) * cos(plhstation(p,2));
rotmat(3,2) = - cos(plhstation(p,1)) * sin(plhstation(p,2));
rotmat(3,3) = - sin(plhstation(p,1))                     ;

  
  neh  = rotmat * [xsat ysat zsat]';
  dist = sqrt(neh(1,:).^2 + neh(2,:).^2 + neh(3,:).^2);
  elev(p,:) = (asin(-neh(3,:)./dist))*180/pi;
  
end;

% --------------------------------
% --- End of function cpelev1  ---
% --------------------------------
