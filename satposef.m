function [xsat, vsat, asat]=satpos(prn,t,eph)
%SATPOS  Compute GPS satellite Position, Velocity and Acceleration.
%        SATPOS returns the satellite position, velocity and acceleration 
%        of the GPS satellite('s) with pseudo-random noise number prn at 
%        time t in a inertial "WGS-84" frame.
%        Syntax:
%                [xsat, vsat, asat]=satpos(prn,t,eph)
%        with
%                prn      array   column vector with satellite prn numbers 
%                t        array   column vector with the time in seconds into
%                                 the GPS week
%                eph      matrix  matrix with ephemeris parameters from RXNAV.
%                -------
%                xsat     matrix  satellite position [ x, y, z ]
%                vsat     matrix  satellite velocity [ vx, vy, vz ]
%                asat     matrix  satellite acceleration [ ax, ay ,az ] 
%
%        In case T is a single number, and PRN a vector, the same time is 
%        used for all computations. In case T is an array, and PRN a
%        single number, the same satellite is used. 
%        The computed position, velocity and acceleration is stored
%        in the rows of xsat, vsat, asat. First for all PRN at T(1), next
%        all PRN at T(2), etc. 
%
%        See also SATCLK, RXNAV and RXOBSE.

%        The satellite ephemeris data is stored in eph, idx contains
%        the row number(s) with ephemeris corresponding to the satellite(s)
%        and time(s) given in prn and t.
%        The index vector idx can be constructed with the function
%
%             idx=seleph(prn,t,eph)
%
%        See also SELEPH.

%        H. van der Marel, LGR, 29-04-95
%        (c) Geodetic Computing Centre, TU Delft

%
% Declare global GPSWEEK
%

global GPSWEEK;

% Some options
%
inert=0;  %  -->  0  Earth fixed WGS-84 
	  %  -->  1  rotating WGS-84
%
% Some constants
%
omgedo = 7292115.1467e-11;     % angular velocity of Earth
gm     = 398.60044e12;         % gravitational constant
c      = 299792458;            % speed of light
epsec  = 1e-12;                % stop criterion for iterations
maxit  = 10;                   % maximum iterations for E

%
% Broadcast Satellite Ephemeris data stored in Eph(_,:)
%
%            name      description
%
% Eph(_, 1)  prn       Satellite prn number
% Eph(_, 2)  toc       Time of clock (seconds into GPS week)
% Eph(_, 3)  svbi      SV clock bias (seconds)
% Eph(_, 4)  svdr      SV clock drift (sec/sec) 
% Eph(_, 5)  svdrr     SV clock drift rate (sec/sec2)
% Eph(_, 6)  aode      AODE (age of data ephemeris) (sec)
% Eph(_, 7)  crs       C_rs (meters)  
% Eph(_, 8)  deltan    delta mean motion (radians/sec)
% Eph(_, 9)  m0        M_0 (radians)
% Eph(_,10)  cuc       C_uc (radians)
% Eph(_,11)  e         eccentricity
% Eph(_,12)  cus       C_us (radians)
% Eph(_,13)  a         semi-major axis (meters)
% Eph(_,14)  toe       Time of ephemeris (seconds into GPS week)
% Eph(_,15)  cic       C_ic (radians)
% Eph(_,16)  omega_0   Omega_0 (radians)
% Eph(_,17)  cis       C_is (radians)
% Eph(_,18)  inc_0     i_0 (radians)
% Eph(_,19)  crc       C_rc (meters)
% Eph(_,20)  s_omega   omega (radians)
% Eph(_,21)  omega_dot omega dot (radians/sec)
% Eph(_,22)  inc_dt    IDOT (radians/sec)
% Eph(_,23)  codl2     codes on L2 channel
% Eph(_,24)  week      GPS week # (for TOE and TOC)
% Eph(_,25)  fl2p      L2 P data flag
% Eph(_,26)  svac      SV accuracy
% Eph(_,27)  svhe      SV health (MSB only)
% Eph(_,28)  tgd       TGD (seconds)
% Eph(_,29)  aodc      AODC (seconds)
% Eph(_,30)  how       transmission time of message (seconds into GPS week)


%
% Compute all requested satellite positions, velocities and accelerations
%

lprn = length(prn);
ltim = length(t);

if lprn==ltim, l=lprn;, else, l=lprn*ltim;, end

for count=1:l
  if lprn==ltim
    itim=count;
    iprn=count;
  else
    itim = fix((count-1)/lprn) + 1; 
    iprn= count - (itim-1)*lprn; 
  end

  time     = t(itim);

  idx      = seleph(prn(iprn),time,eph);
  if idx==-1, disp('No ephemeris data found');break, end

  prni     = eph(idx, 1);
  toc      = eph(idx, 2);
  svbi     = eph(idx, 3);
  svdr     = eph(idx, 4);
  svdrr    = eph(idx, 5);
  aode     = eph(idx, 6);
  crs      = eph(idx, 7);
  deltan   = eph(idx, 8);
  m0       = eph(idx, 9);
  cuc      = eph(idx,10);
  e        = eph(idx,11);
  cus      = eph(idx,12);
  a        = eph(idx,13).^2;
  toe      = eph(idx,14);
  cic      = eph(idx,15);
  omega_0  = eph(idx,16);
  cis      = eph(idx,17);
  inc_0    = eph(idx,18);
  crc      = eph(idx,19);
  s_omega  = eph(idx,20);
  omega_dot= eph(idx,21);
  inc_dt   = eph(idx,22);
  codl2    = eph(idx,23);
  week     = eph(idx,24);
  fl2p     = eph(idx,25);
  svac     = eph(idx,26);
  svhe     = eph(idx,27);
  tgd      = eph(idx,28);
  aodc     = eph(idx,29);
% how      = eph(idx,30);    % not used (only defined in version 2)

% Delta time (in sec)

  dt = time - toe + 604800*(GPSWEEK-week);

% Eccentric anomaly 

%     n      mean motion (mean angular satellite velocity)
%     eccan  eccentic anomaly

  n = sqrt(gm/(a^3)) + deltan;
  m = m0 + n*dt;
  eccan = m;

  for iter=1:maxit
    oldan = eccan;
    eccan = m + e*sin(oldan);
    if abs(oldan-eccan)<epsec, break, end
  end 

% Satellite clock error (including relativistic correction and Tgd)

  dtc        =  time - toc + 604800*(GPSWEEK-week);
  dtrel      = -2*sqrt(gm)*e*sqrt(a)*sin(eccan)/(c^2);
  dts(count,1) = svbi + svdr*dtc + svdrr*(dtc^2) + dtrel - tgd;

% Position in the osculating plane

%     nu     true anomaly
%     phi    argument of perigee + true anomaly
%     u      argument of latitude
%     r      radius
%     xo     position in the orbital plane

  nus   = sqrt(1-e^2) * sin(eccan);
  nuc   = cos(eccan)-e;
  nue   = 1-e*cos(eccan);

  nu    = atan2(nus,nuc);
  phi   = nu + s_omega;
  h2    = [cos(2*phi) ; sin(2*phi)];

  u     = phi     + [cuc cus] * h2;
  r     = a * nue + [crc crs] * h2;

  h     = [cos(u), sin(u) ];

  xo    = r * h;

% Position in WGS-84

%     inc    inclination of the orbital plane
%     omeg   longitude of ascending node
%     p      1th row of rot. matrix orbital -> geocentric
%     q      2nd row of rot. matrix orbital -> geocentric
%     qq     3rd row of rot. matrix orbital -> geocentric
%     xgc    cartesian coordinates of the satellite

  inc   = inc_0 + inc_dt*dt  + [cic cis] * h2;
  if inert==1
    omeg = omega_0 + omega_dot * dt;
  else
    omeg = omega_0 + omega_dot * dt - omgedo * (toe+dt);
  end
  p  = [   cos(omeg)          , sin(omeg)          ,   0      ];
  q  = [  -cos(inc)*sin(omeg) , cos(inc)*cos(omeg) , sin(inc) ];
  qq = [   sin(inc)*sin(omeg) ,-sin(inc)*cos(omeg) , cos(inc) ];

  xsat(count,:) = xo*[p;q];

  if nargout>1

% Velocity
      
%     eccand first derivative of eccentric anomaly
%     nud    first derivative of true anomaly
%     ud     first derivative of argument of latitude
%     rd     first derivative of radius
%     xod    first derivatives of position in orbital plane
%     incd   first derivative of inclination of orbital plane
%     omegd  first derivative of longitude of ascending node
%     pd     first derivative of p
%     qd     first derivative of q
%     Vsat   first der. of cartesian coordinates of satellite

  eccand = n / nue;
  nud    = sqrt(1-e^2) / nue * eccand;
  h2d    = 2 * nud * [ -h(2) ; h(1)] ;

  ud     = nud                   + [cuc cus ] * h2d;
  rd     = a*e*sin(eccan)*eccand + [crc crs ] * h2d;

  hd     = [-h(2) h(1)] ;

  xod    = rd*h + (r*ud)*hd;

  incd   = inc_dt + [cic cis] * h2d;
  if inert==1
    omegd = omega_dot;
  else
    omegd = omega_dot - omgedo;
  end

  pd = [ -p(2) , p(1),  0 ] * omegd;
  qd = [ -q(2) , q(1),  0 ] * omegd + qq *incd;

  vsat(count,:) = xo*[pd;qd] + xod*[p;q];

  end
  
  if nargout>2
      
% Accelleration

%     nudd   second derivative of true anomaly
%     udd    second derivative of argument of latitude
%     rdd    second derivative of radius
%     xodd   second derivative of position in orbital plane
%     incdd  second derivative of inclination of orbital plane
%     pdd    second derivative of p
%     qdd    second derivative of q
%     asat   second der. of cartesian coordinates of satellite

  nudd   = -2 * e * nus / nue^2 * eccand^2;
  h2dd  = 2 * nudd * [ -h(2); h(1) ] - 4*nud^2 * h2 ;

  udd    = nudd + [cuc cus ] * h2dd;
  rdd    = a * e * nuc / nue * eccand^2 + [crc crs ] * h2dd;

  xodd   = (rdd-r*ud^2) * h  + (2*rd*ud+r*udd) * hd;
      
  incdd  =  [cic cis ] * h2dd;

  pdd = - p * omegd^2;
  qdd = - [q(1),q(2),0]*omegd^2 - q*incd^2 + [-qq(2),qq(1),0]*2*incd*omegd ...
    + qq*incdd;

  asat(count,:) = xo*[pdd;qdd] + 2 * xod * [pd;qd] + xodd * [p;q];

  end
  
end
