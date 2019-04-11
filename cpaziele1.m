function [azim,elev,rotmat] = cpaziele1 (xsat,ysat,zsat,xyzstation,tsat);
%CPAZIELE1: Compute elevation/azimuth of satellites/pseudolites
%
% This routine computes satellite positions in an erath-fixed WGS'84
% frame. If station coordinates are specified (in WGS'84), elevations and
% azimuths are computed as well.
%
% Syntax:
%    [azim,elev,rotmat] = cpaziele1 (xsat,ysat,zsat,xyzstation,tsat)
%
% Input arguments:
%    xsat       - xsat(i,j) = X-coordinate (WGS'84) for satellite (i) on epoch (j)
%    ysat       - ysat(i,j) = Y-coordinate (WGS'84) for satellite (i) on epoch (j)
%    zsat       - zsat(i,j) = Z-coordinate (WGS'84) for satellite (i) on epoch (j)
%    plhstation - WGS'84 coordinates for station (phi,lambda,height)
%    xyzstation - WGS'84 coordinates for station (x,y,z)
%    tsat       - Array with times for which positions are to be computed
%
% Output arguments:
%    azim       - azim(i,j) = elevation for satellite (i) on epoch (j)
%    elev       - elev(i,j) = elevation for satellite (i) on epoch (j)
%    rotmat     - Rotation-matrix from WGS'84/XYZ to WGS'84/PLH
% 

% ----------------------------------------------------------------------
% File.....: cpaziele.m
% Date.....: 25-MAY-1999
% Version..: 1.0sv
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% Edited by: Sandra Verhagen
%            11-DEC-2000
% ----------------------------------------------------------------------

% ---------------------------
% --- Initialize matrices ---
% ---------------------------

azim = zeros(size(xsat,1),length(tsat));
elev = zeros(size(xsat,1),length(tsat));

plhstation = xyz2plh (xyzstation);
% -------------------------------------------------
% --- Compute azimuths/elevations, if requested ---
% -------------------------------------------------
xsat = xsat - xyzstation(1);
ysat = ysat - xyzstation(2);
zsat = zsat - xyzstation(3);
 
rotmat(1,1) = - sin(plhstation(1)) * cos(plhstation(2));
rotmat(1,2) = - sin(plhstation(1)) * sin(plhstation(2));
rotmat(1,3) =   cos(plhstation(1))                     ;
rotmat(2,1) = -                      sin(plhstation(2));
rotmat(2,2) =                        cos(plhstation(2));
rotmat(2,3) =   0d0;
rotmat(3,1) = - cos(plhstation(1)) * cos(plhstation(2));
rotmat(3,2) = - cos(plhstation(1)) * sin(plhstation(2));
rotmat(3,3) = - sin(plhstation(1))                     ;

for i = 1:size(xsat,1);
  
  neh  = rotmat * [xsat(i,:)' ysat(i,:)' zsat(i,:)']';
  dist = sqrt(neh(1,:).*neh(1,:) + neh(2,:).*neh(2,:) + neh(3,:).*neh(3,:));
  azim(i,:) = mod(rad2deg(atan2 (neh(2,:),neh(1,:))),360);
  elev(i,:) = rad2deg(asin(-neh(3,:)./dist));
end;

% --------------------------------
% --- End of function cpaziele ---
% --------------------------------
