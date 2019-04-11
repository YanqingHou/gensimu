function [mdbs, snrx, nabx, naby, nabz, nabt] = ...
	cprels (tsat,xsat,ysat,zsat,station,elev,cutoff);
%CPRELS: Compute some reliability measures
%
% The function computes some of the existing reliability measures
% as a function of time.
%
% Syntax:
%    [mdbs, snrx, nabx, naby, nabz, nabt] = cprels ...
%       (tsat,xsat,ysat,zsat,station,elev,cutoff);
%
% Input arguments:
%    tsat    - Epochs (in seconds into the current GPSWEEK)
%    xsat    - X-positions of satellites (WGS'84, Earth-fixed)
%    ysat    - Y-positions of satellites (WGS'84, Earth-fixed)
%    zsat    - Z-positions of satellites (WGS'84, Earth-fixed)
%    station - Station coordinates (WGS'84, Earth-fixed)
%    elev    - Elevation of each satellite on eah time in tsat
%    cutoff  - Cutoff-elevation
%
% Output arguments:
%    mdbs    - Internal reliability, minimal detectable bias
%    snrx    - External reliability, norm of the influence of an
%              observational error the size of the mdb on all unknown 
%              parameters (sqare root)
%    nabx    - External reliability, effect of an error the size of the
%              mdb on the X-coordinate
%    naby    - External reliability, effect of an error the size of the
%              mdb on the Y-coordinate
%    nabz    - External reliability, effect of an error the size of the
%              mdb on the Z-coordinate
%    nabt    - External reliability, effect of an error the size of the
%              mdb on the time-parameter

% ----------------------------------------------------------------------
% File.....: cprels.m
% Date.....: 12-OCT-2000
% Version..: 2.0
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

% -----------------------------
% --- Initialize the arrays ---
% -----------------------------

% ---------------------------
% --- Initialize matrices ---
% ---------------------------

mdbs = NaN * zeros(size(elev,1),length(tsat));
snrx = NaN * zeros(size(elev,1),length(tsat));
nabx = NaN * zeros(size(elev,1),length(tsat));
naby = NaN * zeros(size(elev,1),length(tsat));
nabz = NaN * zeros(size(elev,1),length(tsat));
nabt = NaN * zeros(size(elev,1),length(tsat));

% -----------------------------
% --- Initialize parameters ---
% -----------------------------

lam0     = 17.0747;

% ------------------------------------
% --- Compute reliability measures ---
% ------------------------------------

for i = 1:length(tsat);
 
  [idx1,amat] = cpdesign (xsat(:,i),ysat(:,i),zsat(:,i),station,elev(:,i),cutoff);
  
  if length (idx1) >= 4;

    qxhat = inv(amat' * amat);
    qyhat = amat * qxhat * amat';
    qehat = eye(size(qyhat))-qyhat;

    mdbs(idx1,i) = sqrt (lam0*ones(size(amat,1),1) ./ diag(qehat));
    snrx(idx1,i) = sqrt (mdbs(idx1,i) .^ 2-lam0*ones(size(amat,1),1));
    tnab         = abs(qxhat*amat'*diag(mdbs(idx1,i)));
    nabx(idx1,i) = tnab(1,:)';
    naby(idx1,i) = tnab(2,:)';
    nabz(idx1,i) = tnab(3,:)';
    nabt(idx1,i) = tnab(4,:)';
    
  end;
  
end;
