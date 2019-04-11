function [gdop,pdop,tdop,hdop,fdop,ldop,vdop] = ...
	cpdops (tsat,xsat,ysat,zsat,station,elev,cutoff,rotmat);
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

% -----------------------------
% --- initialize the arrays ---
% -----------------------------

gdop = NaN * zeros (size(tsat));
pdop = NaN * zeros (size(tsat));
tdop = NaN * zeros (size(tsat));
vdop = NaN * zeros (size(tsat));
hdop = NaN * zeros (size(tsat));
fdop = NaN * zeros (size(tsat));
ldop = NaN * zeros (size(tsat));

% ---------------------
% --- Compute DOP's ---
% ---------------------

for i = 1:length(tsat);
  
  [idx1,amat] = cpdesign (xsat(:,i),ysat(:,i),zsat(:,i),station,elev(:,i),cutoff);
  
  if length (idx1) >= 4;
    
    qx      = inv(amat' * amat);
    
    gdop(i) = sqrt(trace(qx));
    pdop(i) = sqrt(trace(qx(1:3,1:3)));
    tdop(i) = sqrt(qx(4,4));

    qx = rotmat * qx(1:3,1:3) * rotmat';
    
    hdop(i) = sqrt(qx(1,1)+qx(2,2));
    fdop(i) = sqrt(qx(1,1));
    ldop(i) = sqrt(qx(2,2));
    vdop(i) = sqrt(qx(3,3));
    
  end;
  
end;