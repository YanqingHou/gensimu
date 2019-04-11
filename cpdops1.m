function [gdop,pdop,tdop] = ...
	cpdops1 (xsat,ysat,zsat,station);
%CPDOPS: Compute the different existing DOP's
%
% The function computes the different existing DOP'values, including
% GDOP, PDOP, TDOP, HDOP, FDOP, LDOP and VDOP as a function of time.
%
% Syntax:
% [gdop,pdop,tdop,hdop,fdop,ldop,vdop] = cpdops ...
%    (tsat,xsat,ysat,zsat,station,elev,cutoff,rotmat);
%
% Input arguments:
%    tsat    - Times (in seconds into the current GPSWEEK)
%    xsat    - X-positions of satellites (WGS'84, Earth-fixed)
%    ysat    - Y-positions of satellites (WGS'84, Earth-fixed)
%    zsat    - Z-positions of satellites (WGS'84, Earth-fixed)
%    station - Station coordinates (WGS'84, Earth-fixed)
%    elev    - Elevation of each satellite on eah time in tsat
%    cutoff  - Cutoff-elevation
%    rotmat  - Rotation-matrix from WGS'84/XYZ to WGS'84/PLH
%
% Output arguments:
%    gdop    - Geometrical Dilution of Precision
%    pdop    - Positional  Dilution of Precision
%    tdop    - Time Dilution of Precision
%    hdop    - Horizontal Dilution of Precision
%    fdop    - Lattitude Dilution of Precision
%    ldop    - Longitude Dilution of Precision
%    vdop    - Vertical Dilution of Precision
%

% ----------------------------------------------------------------------
% File.....: pldops.m
% Date.....: 12-OCT-2000
% Version..: 1.0
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

% ---------------------
% --- Compute DOP's ---
% ---------------------

for i = 1:length(tsat);
  
  [amat] = Amatrix (xsat,ysat,zsat,station,m);
  
  if m >= 4;
    
    qx      = inv(amat' * amat);
    
    gdop(i) = sqrt(trace(qx));
    pdop(i) = sqrt(trace(qx(1:3,1:3)));
    tdop(i) = sqrt(qx(4,4));
    
  end;
  
end;