function [idx1, amat] = cpdesign (xsat,ysat,zsat,station,elev,cutoff);
%CPDESIGN: Create the design-matrix
%
% The function creates the design-matrix, which is necessary for the
% computation of DOP-values and reliability measures. Note that this
% function computes the design matrix for a single epoch, therefore
% the input-arguments (xsat,ysat,zsat,station,elev) should be given
% for that single epoch as well.
%
% Syntax:
%    [idx1, amat] = cpdesign (xsat,ysat,zsat,station,elev,cutoff);
%
% Input arguments:
%    xsat    - X-positions of satellites (WGS'84, Earth-fixed)
%    ysat    - Y-positions of satellites (WGS'84, Earth-fixed)
%    zsat    - Z-positions of satellites (WGS'84, Earth-fixed)
%    station - Station coordinates (WGS'84, Earth-fixed)
%    elev    - Elevation of each satellite on eah time in tsat
%    cutoff  - Cutoff-elevation
%
% Output arguments:
%    idx1    - Satellites in view (above cutoff-elevation)
%    amat    - Design-matrix

% ----------------------------------------------------------------------
% File.....: cpdesign.m
% Date.....: 12-OCT-2000
% Version..: 2.0
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

xsat = xsat - station(1);
ysat = ysat - station(2);
zsat = zsat - station(3);
dist = sqrt(xsat.*xsat + ysat.*ysat + zsat.*zsat);

idx1 = find(elev>cutoff);
amat = zeros (length(idx1),4);
  
for j = 1:length(idx1);
    
  amat(j,1) = xsat(idx1(j))/dist(idx1(j));
  amat(j,2) = ysat(idx1(j))/dist(idx1(j));
  amat(j,3) = zsat(idx1(j))/dist(idx1(j));
  amat(j,4) = -1;
 
end;
