function [PDOP,M] = dopworld(file,date,phi,lam)

% PDOP = dopworld (file, date,phi,lam) compute worldwide pdop values
%
% INPUT:
% file : YUMA ephemerides file 'yumaXX.txt' (XX is weeknumber)
% date : date and time 'dd-MONTH-year hh:mm:ss.sss' (MONTH: 3 first letters of month)
% phi  : latitude, vector in radians (default pi/2 : -pi/2)
% lam  : longitude, vector in radians (default -pi : pi)
%
% OUTPUT:
% PDOP : matrix with pdops for each gridcell
%_____________________________________________________________________________________

% define grid

if nargin < 3
   grid = pi/72;
   phi = pi/2 : -grid : -pi/2;
   lam = -pi : grid : pi;
end

pmax = size(phi,2);
lmax = size(lam,2);

% start with ephemerides data, epoch time, positions of satellites

eph = rdyuma(file);

tsat = mktsat(date,date,300);

[xs,ys,zs] = cpaziele(tsat,eph);

% give cutoff elevation (degrees)

cutoff = 10; 

% initialize PDOP matrix

PDOP = zeros(pmax,lmax);

% start computation

for i = 1:pmax
   
   for j = 1:lmax
      
      plh = [phi(i) lam(j) 0];
      xyz = plh2xyz(plh);
      [az,el] = cpaziele1 (xs,ys,zs,xyz,tsat);
      [m,gdop,pdop] = cpdops1(tsat,xs,ys,zs,xyz,el,cutoff);
      PDOP(i,j) = pdop;
      M(i,j) = m;
      
   end
   
   a=i+1
   
end

p = find(PDOP>7);
PDOP(p) = 7*ones(size(p));

pldopworld(eph,tsat,xs,ys,zs,PDOP);

% end of computation