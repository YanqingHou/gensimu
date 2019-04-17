function [xsat]=satposef(prn,curweek,t,eph,sys,mode)
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

% global GPSWEEK;

% Some options
%
% inert=0;  %  -->  0  Earth fixed WGS-84
% 	  %  -->  1  rotating WGS-84
% %
% % Some constants
% %
% omgedo = 7292115.1467e-11;     % angular velocity of Earth
% gm     = 398.60044e12;         % gravitational constant
% c      = 299792458;            % speed of light
% epsec  = 1e-12;                % stop criterion for iterations
% maxit  = 10;                   % maximum iterations for E

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
% cordform='Inertial';
cordform='ECEF';

lprn = length(prn);
ltim = length(t);
xsat=nan(ltim,3);
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
    
    idx      = seleph(prn(iprn),curweek,time,eph);
    if idx==-1, disp('No ephemeris data found');break, end
    cursec=struct('sec',0,'dec',0);
    cursec.sec=time;
    ephi=eph(idx,:);
    if strcmp(mode,'alm')
        [~,rsat,~]=Alm2Pos(curweek,cursec,ephi,sys,prn(iprn),cordform);
    elseif strcmp(mode,'eph')
        [~,rsat,~]=Eph2Pos(curweek,cursec,ephi,sys,prn(iprn),cordform);
    else
        error('unkonwn mode!');
    end
    xsat(count,:)=rsat';
end
