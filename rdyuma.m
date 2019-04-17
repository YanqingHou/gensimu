function [eph] = rdyuma (file,addid,ymaperiod)
%RDYUMA: Read satellite ephemeris from a YUMA-file
%
% This routine reads satellite data from a YUMA almanac. It can be used for
% either GPS satellites or GLONASS satellites. Data is stored in an array
% "eph", as defined by Hans van der Marel, and explained below. Note that
% not all items in the "eph" definition are present in YUMA-almanacs. Those
% items will be set to 0 (zero) by this routine.
%
% GLONASS almanacs sometimes store the frequency-number as well as the
% slot-number in the "ID:" field. In that case, the slot number is assumed
% to be in the first two positions, and is used as the ID of the satellite.
%
% Syntax:
%    function [eph] = rdyuma (file,addid,yumaperiod);
%
% Input arguments:
%    file  : Filename with YUMA-almanac
%    addid : Number to be added to the ID of a satellite. This can be used
%            for example to distinguish between GPS and GLONASS satellites
%            (for GLONASS, set to "32"). If left out, nothing is added and
%            the ID remains as it was read from the YUMA file.
%   yumaperiod: period of 1024 weeks of the ymaalmanak.
%   period 1= 6jan 1980 - 21 aug 1999
%   period 2= 22 aug 1999 - 7 apr 2019
%   period 3= 8 apr 2019 - 21 nov 2038
%   Default (if not specified) is period from commando clock (now)
%
% Output arguments:
%    eph   : Ephemeris as found in YUMA-almanac
%
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

% ----------------------------------------------------------------------
% File.....: rdyuma.m
% Date.....: 30-NOV-1999
% Version..: 2.0
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

% ----------------------------------------------------------
% --- Declare global variables for the GPSWEEK / GPSECHO ---
% ----------------------------------------------------------

% global GPSWEEK;
global GPSECHO;

% ------------------------------------
% --- Value to add to satellite ID ---
% ------------------------------------
if ~exist('addid')
    addid = 0;
end
if isempty(addid)
    addid = 0;
end

if ~exist('ymaperiod')
    time        = clock;
    ymaperiod   = ceil(((ymd2mjd(time(1:3))-ymd2mjd([1980 1 6]))/7)/1024);
end

% ----------------------
% --- Open YUMA-file ---
% ----------------------

fid=fopen(file);
if fid==-1
    disp(['Error opening YUMA ephemeris file: ' file]);
    return;
end

% -------------------------------------
% --- Read ephemeris and close file ---
% -------------------------------------

count = 0;

while feof(fid) == 0
    
    line=upper(fgetl(fid));
    
    if ischar(line) & length(line) >= 3
        
        if ~isempty(GPSECHO)
            disp(line);
        end
        
        switch (line(1:3))
            
            case 'ID:'
                count = count + 1;
                eph(count,1:30) = 0;
                eph(count, 1) = str2num (line(28:length(line)));
                if eph(count,1) > 100, eph(count,1) = floor(eph(count,1)/100); end
                eph(count,1) = eph(count,1) + addid;
            case 'HEA'
                eph(count,27) = str2num (line(28:length(line)));
            case 'ECC'
                eph(count,11) = str2num (line(28:length(line)));
            case 'TIM'
                eph(count,14) = str2num (line(28:length(line)));
            case 'ORB'
                eph(count,18) = str2num (line(28:length(line)));
            case 'RAT'
                eph(count,21) = str2num (line(28:length(line)));
            case 'SQR'
                eph(count,13) = str2num (line(28:length(line)));
            case 'RIG'
                eph(count,16) = str2num (line(28:length(line)));
            case 'ARG'
                eph(count,20) = str2num (line(28:length(line)));
            case 'MEA'
                eph(count, 9) = str2num (line(28:length(line)));
            case 'Af0'
                eph(count, 3) = str2num (line(28:length(line)));
            case 'Af1'
                eph(count, 4) = str2num (line(28:length(line)));
            case 'WEE'
                eph(count,24) = str2num (line(28:length(line)));
                eph(count,24) = eph(count,24) + (ymaperiod-1)*1024;
            otherwise;
        end
        
    end
    
end

% --------------------------------
% --- Set GPSWEEK if necessary ---
% --------------------------------

% if length(GPSWEEK)==0
%   GPSWEEK=eph(1,24);
% end
% 北斗时间转为GPS时间
if addid==260
    eph(:,14)=eph(:,14)+14;
    eph(:,24)=eph(:,24)+1356;
    %如果秒数超过一周，则调整周数
    extrawk=eph(:,14)>=604800;
    eph(extrawk,24)=eph(extrawk,24)+1;
    eph(extrawk,14)=eph(extrawk,14)-604800;
end
% 在历书中，toa=toe
eph(:,2)=eph(:,14);
% --------------------------------------
% --- Echo some results if necessary ---
% --------------------------------------

if length(GPSECHO)~=0
    disp(['Ready reading Navigation data from YUMA-dataset: ' file]);
    disp ([num2str(count) ' ephemeris sets found']);
end
