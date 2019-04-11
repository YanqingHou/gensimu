function xyz = plh2xyz(plh)
%PLH2XYZ Ellipsoidal coordinates to Cartesian Coordinates 
%        Converts ellipsoidal coordinates Phi, Lambda and h into
%        cartesian coordinates X, Y and Z:
%
%           xyz = plh2xyz(plh,ellips)
%
%        Ellips is a text string with the name of the ellipsoid or
%        a vector with the semi-major axis a and flattening 1/f.
%        Default for ellips is 'WGS-84'.

%        H. van der Marel, LGR, 07-05-95
%        (c) Geodetic Computing Centre, TU Delft


if nargin<2, ellips='WGS-84';, end
if isstr(ellips)
  [a,f] = inqell(ellips);
else
  a=ellips(1);
  f=1/ellips(2);
end
% excentricity e (squared) 
e2 = 2*f - f^2;

[m,n]=size(plh);
if n~=3 & m==3, plh=plh';, end

N = a ./ sqrt(1 - e2 .* sin(plh(:,1)).^2);
xyz = [ (N+plh(:,3)).*cos(plh(:,1)).*cos(plh(:,2)) ...
        (N+plh(:,3)).*cos(plh(:,1)).*sin(plh(:,2)) ...
        (N-e2.*N+plh(:,3)).*sin(plh(:,1))          ];

if n~=3 & m==3, xyz=xyz';, end

