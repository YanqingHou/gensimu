function [ eph, alpha, beta, utc, leap ] = rxnav(file)
%RXNAV   Read RINEX Navigation File 
%        RXNAV(file) reads a RINEX Navigation File and return ephemeris data 
%        EPH, ionosphere parameters ALPHA and BETA, UTC conversion parameters 
%        A0, A1, T and W, and leap second LEAP.
%        Syntax:
%                 [eph, alpha, beta, utc, leap ] = rxnav(file)
%        with:
%                  file       string   name of the RINEX meteo file
%                  eph        matrix   table with ephemeris parameters; a set 
%                                      of ephemerides parameters is stored in a
%                                      row of eph.
%                  alpha,beta array    ionosphere parameters 
%                  utc        array    utc conversion parameters ao,a1,t,w 
%                  leap       number   leap second 
%

%        H. van der Marel, LGR, 29-04-95
%        (c) Geodetic Computing Centre, TU Delft

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

% Declare global variable for the GPSWEEK

global GPSWEEK
global GPSECHO

fid=fopen(file);
if fid==-1, disp(['Error opening RINEX Navigation file ' file]);return, end
%
% read navigation file header
%
while feof(fid)~=1
  line=fgetl(fid);
  if ~isstr(line),break,end
  if length(GPSECHO)~=0, disp(line), end
  if length(findstr(line(61:length(line)),'RINEX VERSION / TYPE'))~=0
    version=str2num(line(1:6));
    if strcmp(line(21:21),'N')==0, disp('This is not a RINEX Navigation file');return,end
  end
  if length(findstr(line(61:length(line)),'ION ALPHA'))~=0
     alpha = sscanf(line(1:60),'%g',4);
  end
  if length(findstr(line(61:length(line)),'ION BETA'))~=0
     beta  = sscanf(line(1:60),'%g',4);
  end
  %if length(findstr(line(61:length(line)),'DELTA-UTC: A0,A1,T,W'))~=0
  %   utc(1) = str2num(line( 4:22));
  %   utc(2) = str2num(line(23:41));
  %   utc(3) = str2num(line(42:50));
  %   utc(4) = str2num(line(51:59));
  %end
  if length(findstr(line(61:length(line)),'LEAP SECONDS'))~=0
    leap=str2num(line(1:6));
  end
  if length(findstr(line(61:length(line)),'END OF HEADER'))~=0,break,end
  if length(findstr(line,blanks(80)))~=0,break,end
end
%
% read ephemeris data
%
count=0;
if version==1, NumRec=7;, else, NumRec=8;, end  
while feof(fid)==0
  count=count+1;
  for Rec=1:NumRec
    line=fgetl(fid);
    %disp(line);
    if ~isstr(line),count=count-1;break,end
    if Rec==1
       eph(count,1)=str2num(line(1:2));  
       T=sscanf(line(3:22),'%d %d %d %d %d %f',6);
       ic=2;
       i1=2;
       i2=4;
    elseif Rec==8
       i2=1;
    end
    for i=i1:i2
      ic=ic+1;
      ephi=str2num(line(19*(i-1)+4:19*i+3));
      eph(count,ic)=ephi;
    end
    i1=1;
  end
% convert TOC (yy-mm-dd hh:mm:ss) to seconds into the GPS week
  week=eph(count,24);
  days=ymd2mjd(T(1:3)')-gps2mjd([week 0]);
  eph(count,2)=days*86400+T(4)*3600+T(5)*60+T(6);
% set the global GPSWEEK (if necessary)
  if length(GPSWEEK)==0
     GPSWEEK=week;
  end  
end
if length(GPSECHO)~=0, disp(['Ready reading Navigation data: ' num2str(count) ' ephemeris sets found']);, end  
fclose(fid);


