function [eph] = rd_alma (file,addid)
%RD_ALMANAC: Read satellite ephemeris from Trimble almanac file
%
% This routine reads satellite data from a Trimble almanac file. It can be used for
% ALL GNSS satellites. Data is stored in an array
% "eph", as defined by Hans van der Marel, and explained below. Note that
% not all items in the "eph" definition are present in almanacs. Those
% items will be set to 0 (zero) by this routine.
%
% GLONASS almanacs sometimes store the frequency-number as well as the
% slot-number in the "ID:" field. In that case, the slot number is assumed
% to be in the first two positions, and is used as the ID of the satellite.
%
% Syntax:
%    function [eph] = rd_alma (file,addid);
%
% Input arguments:
%    file  : Filename with almanac
%    addid : Number to be added to the ID of a satellite. This can be used
%            for example to distinguish between GPS and GLONASS satellites
%            (for GLONASS, set to "32"). If left out, nothing is added and
%            the ID remains as it was read from the ALMANAC file.
%
% Output arguments:
%    eph   : Ephemeris as found in almanac
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
% File.....: rd_alma.m
% Date.....: 30-July-2013
% Version..: 2.0
% Author...: Jingyu Zhang
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

global GPSWEEK;

% ------------------------------------
% --- Value to add to satellite ID ---
% ------------------------------------
if ~exist('addid')
   addid = 0;
end
if isempty(addid)
   addid = 0;
end

% ----------------------
% --- Open ALMA-file ---
% ----------------------

fid=fopen(file);
if fid==-1;
  disp(['Error opening Almanac file: ' file]);
  return;
end;

% -------------------------------------
% --- Read ephemeris and close file ---
% -------------------------------------

count = 0;
eph=zeros(1,30);
num=1;           % To indicate the number of line
while feof(fid) == 0;
    
    line=upper(fgetl(fid));
    count=count+1;
    row=mod(count,14);
    if row~=0
        switch row
            
            case 1
                n=0; 
                while length(line)/n~=10
                    eph(num,1)=str2num(line(1+10*n:10*(n+1)));
                    n=n+1;
                    num=num+1;
                end
                
            case 2
                mm=n-1;
                col=0;
                while mm~=-1
                    leph=size(eph,1);
                    if line(1,1+10*col:10*(col+1))=='         1'
                        eph(leph-mm,27)=str2num(line(1,1+10*col:10*(col+1)));
                    else
                        eph(leph-mm,27)=0;
                    end
                    mm=mm-1;
                    col=col+1;
                end
                
            case 3
                mm=n-1;
                col=0;
                while mm~=-1
                    leph=size(eph,1);
                    eph(leph-mm,11)=str2num(line(1,1+10*col:10*(col+1)));
                    mm=mm-1;
                    col=col+1;
                end
                
            case 4
                mm=n-1;
                col=0;
                while mm~=-1
                    leph=size(eph,1);
                    eph(leph-mm,13)=str2num(line(1,1+10*col:10*(col+1)));
                    mm=mm-1;
                    col=col+1;
                end
                
            case 5
                mm=n-1;
                col=0;
                while mm~=-1
                    leph=size(eph,1);
                    eph(leph-mm,16)=str2num(line(1,1+10*col:10*(col+1)))*pi/180;
                    mm=mm-1;
                    col=col+1;
                end
                
            case 6
                mm=n-1;
                col=0;
                while mm~=-1
                    leph=size(eph,1);
                    eph(leph-mm,20)=str2num(line(1,1+10*col:10*(col+1)))*pi/180;
                    mm=mm-1;
                    col=col+1;
                end
                
            case 7
                mm=n-1;
                col=0;
                while mm~=-1
                    leph=size(eph,1);
                    eph(leph-mm,9)=str2num(line(1,1+10*col:10*(col+1)))*pi/180;
                    mm=mm-1;
                    col=col+1;
                end
                
            case 8
                mm=n-1;
                col=0;
                while mm~=-1
                    leph=size(eph,1);
                    eph(leph-mm,14)=str2num(line(1,1+10*col:10*(col+1)));
                    mm=mm-1;
                    col=col+1;
                end
                
            case 9
                mm=n-1;
                col=0;
                while mm~=-1
                    leph=size(eph,1);
                    eph(leph-mm,18)=(str2num(line(1,1+10*col:10*(col+1)))+54)*pi/180;
                    mm=mm-1;
                    col=col+1;
                end
                
            case 10
                mm=n-1;
                col=0;
                while mm~=-1
                    leph=size(eph,1);
                    eph(leph-mm,21)=str2num(line(1,1+10*col:10*(col+1)))*pi*0.001/180;
                    mm=mm-1;
                    col=col+1;
                end
                
            case 13
                mm=n-1;
                col=0;
                while mm~=-1
                    leph=size(eph,1);
                    eph(leph-mm,24)=str2num(line(1,1+10*col:10*(col+1)));
                    mm=mm-1;
                    col=col+1;
                end
                
                
            otherwise
        end
        
    end
    
end

% --------------------------------
% --- Set GPSWEEK if necessary ---
% --------------------------------
 
if length(GPSWEEK)==0
  GPSWEEK=eph(1,24);
end  





