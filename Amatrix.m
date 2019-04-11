function [amat] = Amatrix (xsat,ysat,zsat,station,m);
%CPDESIGN1: Create the design-matrix
%
% The function creates the design-matrix, which is necessary for the
% computation of DOP-values and reliability measures. Note that this
% function computes the design matrix for a single epoch, therefore
% the input-arguments (xsat,ysat,zsat,station,elev) should be given
% for that single epoch as well.
%
% Syntax:
%    [amat] = cpdesign1 (xsat,ysat,zsat,station,m);
%
% Input arguments:
%    xsat    - X-positions of satellites (WGS'84, Earth-fixed)
%    ysat    - Y-positions of satellites (WGS'84, Earth-fixed)
%    zsat    - Z-positions of satellites (WGS'84, Earth-fixed)
%    station - Station coordinates (WGS'84, Earth-fixed)
%    m       - number of satellites
%
% Output arguments:
%    amat    - Design-matrix

% ----------------------------------------------------------------------
% File.....: cpdesign1.m
% Date.....: 12-OCT-2000
% Version..: 2.0sv
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% Edited by: Sandra Verhagen
%            11-DEC-2000
% ----------------------------------------------------------------------



xsat = xsat - station(1);
ysat = ysat - station(2);
zsat = zsat - station(3);
dist = sqrt(xsat.*xsat + ysat.*ysat + zsat.*zsat);
amat = zeros (m,4);
amat(:,1) = -xsat./dist;
amat(:,2) = -ysat./dist;
amat(:,3) = -zsat./dist;

amat(:,4) = 1;

