function [a,f,GM] = inqell(ellips)
%INQELL  Semi-major axis, flattening and GM for various ellipsoids.
%        Returns the semi-major axis a, flattening f and GM for a 
%        particular ellipsoid.
%        Default for ellips is 'WGS-84'.
%        A '?' returns all known ellipsoids.

%        H. van der Marel, LGR, 07-05-95
%        (c) Geodetic Computing Centre, TU Delft

ell= ['AIRY         ';
      'BESSEL       ';
      'CLARKE       ';
      'INTERNATIONAL';
      'HAYFORD      ';
      'GRS80        ';
      'WGS-84       '];
par= [6377563.396 , 299.324964    , NaN        ;
      6377397.155 , 299.1528128   , NaN        ;
      6378249.145 , 293.465       , NaN        ;
      6378388.      297.00        , NaN        ;
      6378388.      297.00        , 3.986329e14;      
      6378137.    , 298.257222101 , 3.986005e14;
      6378137.    , 298.257223563 , 3.986005e14];

if nargin==0,ellips='unknown';, end
ellips=deblank(upper(ellips));
if strcmp(ellips,'?')
  disp('Ellipsoids:');
  disp(' ');
  disp(ell);
  return
end

i=0;
for j=1:size(par,1)
  if strcmp(deblank(ell(j,:)),ellips), i=j;, end
end
if i==0 
  i=size(par,1); 
  disp(['Warning: Ellipsoid ',ellips,' not found, ',deblank(ell(i,:)), ...
  ' selected instead']);
end

a=par(i,1);
f=1/par(i,2);
GM=par(i,3);
