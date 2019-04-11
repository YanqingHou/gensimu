function xyz = plh2xyzwgs(plh)
%PLH2XYZWGS Ellipsoidal coordinates to Cartesian Coordinates 
%           Converts ellipsoidal coordinates Phi, Lambda and h into
%           cartesian coordinates X, Y and Z:
%
%           xyz = plh2xyzwgs(plh)
%
%        Ellips is WGS-84.

%        Based on plh2xyz.m written by
%        H. van der Marel, LGR, 07-05-95
%        (c) Geodetic Computing Centre, TU Delft


a=6378137.;
f=1/298.257223563;


% excentricity e (squared) 
e2 = 2*f - f^2;

[m,n]=size(plh);
if n~=3 & m==3, plh=plh';, end

N = a ./ sqrt(1 - e2 .* sin(plh(:,1)).^2);
xyz = [ (N+plh(:,3)).*cos(plh(:,1)).*cos(plh(:,2)) ...
        (N+plh(:,3)).*cos(plh(:,1)).*sin(plh(:,2)) ...
        (N-e2.*N+plh(:,3)).*sin(plh(:,1))          ];

if n~=3 & m==3, xyz=xyz';, end
    


