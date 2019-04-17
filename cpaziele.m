function [xsat,ysat,zsat,azim,elev,rotmat] = cpaziele (curweek,tsat,eph,station,sys,mode)
%CPAZIELE: Compute satellite positions and elevation/azimuth
%
% This routine computes satellite positions in an erath-fixed WGS'84
% frame. If station coordinates are specified (in WGS'84), elevations and
% azimuths are computed as well.
%
% Syntax:
%    [xsat,ysat,zsat,azim,elev,rotmat] = cpaziele (tsat,eph,station);
%
% Input arguments:
%    tsat    - Array with times for which positions are to be computed
%    eph     - Satellite ephemeris data
%    station - WGS'84 coordinates for station (optional)
%
% Output arguments:
%    xsat    - xsat(i,j) = X-coordinate (WGS'84) for satellite (i) on epoch (j)
%    ysat    - ysat(i,j) = Y-coordinate (WGS'84) for satellite (i) on epoch (j)
%    zsat    - zsat(i,j) = Z-coordinate (WGS'84) for satellite (i) on epoch (j)
%    azim    - azim(i,j) = elevation for satellite (i) on epoch (j)
%    elev    - elev(i,j) = elevation for satellite (i) on epoch (j)
%    rotmat  - Rotation-matrix from WGS'84/XYZ to WGS'84/PLH
% 

% ----------------------------------------------------------------------
% File.....: cpaziele.m
% Date.....: 25-MAY-1999
% Version..: 1.0
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

% ---------------------------
% --- Initialize matrices ---
% ---------------------------
prns=unique(eph(:,1));
satnum=length(prns);
xsat = zeros(satnum,length(tsat));
ysat = zeros(satnum,length(tsat));
zsat = zeros(satnum,length(tsat));
azim = zeros(satnum,length(tsat));
elev = zeros(satnum,length(tsat));
% -------------------------------------------------------
% --- Compute satellite positions (XYZ WGS84 & AZ/EL) ---
% -------------------------------------------------------

for i = 1:satnum%size(eph,1)

  [tmp] = satposef (prns(i),curweek,tsat,eph,sys,mode);
  xsat(i,:) = tmp(:,1)';
  ysat(i,:) = tmp(:,2)';
  zsat(i,:) = tmp(:,3)';
  
end

% -------------------------------------------------
% --- Compute azimuths/elevations, if requested ---
% -------------------------------------------------
% 
% if nargin > 2;
%   
%   xsat = xsat - station(1);
%   ysat = ysat - station(2);
%   zsat = zsat - station(3);
%   
%   plh = xyz2plh (station,'WGS-84');
% 
%   rotmat(1,1) = - sin(plh(1)) * cos(plh(2));
%   rotmat(1,2) = - sin(plh(1)) * sin(plh(2));
%   rotmat(1,3) =   cos(plh(1))              ;
%   rotmat(2,1) = -               sin(plh(2));
%   rotmat(2,2) =                 cos(plh(2));
%   rotmat(2,3) =   0d0;
%   rotmat(3,1) = - cos(plh(1)) * cos(plh(2));
%   rotmat(3,2) = - cos(plh(1)) * sin(plh(2));
%   rotmat(3,3) = - sin(plh(1))              ;
% 
%   for i = 1:size(eph,1);
%     
%     neh  = rotmat * [xsat(i,:)' ysat(i,:)' zsat(i,:)']';
%     dist = sqrt(neh(1,:).*neh(1,:) + neh(2,:).*neh(2,:) + neh(3,:).*neh(3,:));
%     azim(i,:) = mod(rad2deg(atan2 (neh(2,:),neh(1,:))),360);
%     elev(i,:) = rad2deg(asin(-neh(3,:)./dist));
%     
%   end;
% 
%   xsat = xsat + station(1);
%   ysat = ysat + station(2);
%   zsat = zsat + station(3);
% 
% end;

% --------------------------------
% --- End of function cpaziele ---
% --------------------------------
