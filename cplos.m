function [los] = cplos (xyz,station);
%CPLOS: Create the design-matrix for pseudolite-case
%
% Syntax:
%    [los] = cplos (xyz,station);
%
% Input arguments:
%    xyz     - XYZ-positions of pseudo- or satellites (WGS'84, Earth-fixed)
%    station - Station coordinates (WGS'84, Earth-fixed)
%
% Output arguments:
%    los     - line-of-sight vectors from station to pseudo- or satellites
% ________________________________________________________________________
% File.....: cplos.m
% Date.....: 22-DEC-2000
% Version..: 1.0sv
% Author...: Sandra Verhagen
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% _________________________________________________________________________

npl = size(xyz,1);
for i = 1:npl
   xyzt = xyz(i,:) - station;
   dist = sqrt(sum(xyzt.*xyzt)); 
   los(i,:) = -xyzt./dist;
end

 
