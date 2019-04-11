function plh = xyz2plh(xyz,ellips,method)
%XYZ2PLH Cartesian Coordinates to Ellipsoidal coordinates.
%        Converts cartesian coordinates X, Y and Z into ellipsoidal 
%        coordinates Phi, Lambda and h:
% 
%             plh = xyz2plh(xyz,ellips)
%
%        Ellips is a text string with the name of the ellipsoid or
%        a vector with the semi-major axis a and flattening 1/f.
%        Default for ellips is 'WGS-84'.
%        This subroutine uses Bowring's method by default. The more
%        conventional iterative method can be also be used
% 
%             plh = xyz2plh(xyz,ellips,1)
%
%        This method is less precise on the surface of the earth, and should
%        only be used above 10-20 km of height.

%        H. van der Marel, LGR, 07-05-95
%        (c) Geodetic Computing Centre, TU Delft


if nargin<2, ellips='WGS-84';, end
if nargin<3, method=0;, end
if isstr(ellips)
  [a,f] = inqell(ellips);
else
  a=ellips(1);
  f=1/ellips(2);
end
% excentricity e (squared) and semi-minor axis
e2 = 2*f - f^2;
b=(1-f)*a;

[m,n]=size(xyz);
if n~=3 & m==3, xyz=xyz';, end

r  = sqrt(xyz(:,1).^2+xyz(:,2).^2);

if method==1
% compute phi via iteration
  Np = xyz(:,3);
  for i=1:4
    phi = atan( (xyz(:,3) + e2.*Np) ./ r );
    N = a ./ sqrt(1 - e2 .* sin(phi).^2);
    Np = N .* sin(phi);
  end
else
% compute phi using B.R. Bowring's equation (default method)
  u    = atan2 ( xyz(:,3).*a , r.*b ); 
  phi  = atan2 ( xyz(:,3) + (e2/(1-e2)*b) .* sin(u).^3, r - (e2*a) .* cos(u).^3 );
  N = a ./ sqrt(1 - e2 .* sin(phi).^2);
end

plh  = [ phi                       ...
         atan2(xyz(:,2),xyz(:,1))  ...
         r ./ cos(phi) - N          ];

if n~=3 & m==3, plh=plh';, end
    


