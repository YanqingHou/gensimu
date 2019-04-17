function [ outputEphemeris] = readRinexNav( filePath )
%readRinexNav Reads a mixed RINEX navigation file *.nav and returns the
%loaded ephemeris for each constellation
%   Reads Keplerian and Cartesian type ephemeris coming from RINEX 3.02
%   Files can be downlaoded from here: ftp://cddis.gsfc.nasa.gov/gnss/data/campaign/mgex/daily/rinex3/2015/
%   Download in *.p format and convert to .nav using rtklib

%%%%%-------Input
%       fileName = File adress

%%%%%------- Output
%       outputEphemeris = Class containing the ephemeris for each
%       constellation

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

'Loading ephemeris...'
endOfHeader = 0;
ionoAlpha=0;ionoBeta=0;
navFile = fopen(filePath);
global GPSWEEK
%Read header
while (~endOfHeader)
    line = fgetl(navFile);
    lineSplit = strsplit(line);
    
    if strfind(line,'RINEX VERSION')
        Version = lineSplit(2);
        vernum=str2double(Version);
        %         if ~strcmp(Version,'3.02')
        if vernum<3.0
            error 'Not the correct version, should be 3.02'
        end
        
        
    elseif strfind(line,'DATE')
        date = lineSplit(3);
        year = str2double(date{1,1}(1:4));
        month = str2double(date{1,1}(5:6));
        day = str2double(date{1,1}(7:8));
        DOY=date2doy([real(year),real(month),real(day)]);
        %         DOY=Date2DayOfYear(real(year),real(month),real(day));
        
    elseif strfind(line,'IONOSPHERIC CORR')
        if strcmp(lineSplit(1), 'GPSA')
            ionoAlpha = str2double(lineSplit(2:5));
        elseif strcmp(lineSplit(1), 'GPSB')
            ionoBeta = str2double(lineSplit(2:5));
        end
    elseif strfind (line,'LEAP SECONDS')
        leapSeconds = str2double(lineSplit(2));
    elseif strfind(line,'END OF HEADER')
        endOfHeader = 1;
    end
end

%Pointer line set at the end of the header.
ionosphericParameters = [ionoAlpha; ionoBeta];


%read body

gpsEphemeris =  [];
glonassEphemeris = [];
beidouEphemeris = [];

keplerArray = zeros(22,1); %Vector containing Keplerian elements type ephemeris (GPS, Beidou, Galileo)
cartesianArray = zeros(19,1); %Vector containing Cartesian type ephemeris (GLONASS, SBAS)
while ~feof(navFile)
    line = fgetl(navFile);
    %     lineSplit = strsplit(line);
    
    constellation = line(1);
    if ischar(constellation) %New Ephemeris
        switch constellation
            case {'G', 'C'}                %If the ephemeris is ether for GPS or Beidou, store Keplerian elements
                
                %%Read All of the ephemeris
                Time=str2double(strsplit(line(5:23)));
                svprn = str2double([line(2), line(3)]);
                af2=str2double(line(24:42));
                af1=str2double(line(43:61));
                af0=str2double(line(62:end));
                
                line=fgetl(navFile);
                TODE=str2double(line(5:23));
                crs=str2double(line(24:42));
                deltan=str2double(line(43:61));
                M0=str2double(line(62:end));
                
                
                
                line=fgetl(navFile);
                cuc=str2double(line(5:23));
                ecc=str2double(line(24:42));
                cus=str2double(line(43:61));
                roota=str2double(line(62:end));
                
                
                
                line=fgetl(navFile);
                toe=str2double(line(5:23));
                cic=str2double(line(24:42));
                Omega0=str2double(line(43:61));
                cis=str2double(line(62:end));
                
                
                
                line=fgetl(navFile);
                i0=str2double(line(5:23));
                crc=str2double(line(24:42));
                omega=str2double(line(43:61));
                Omegadot=str2double(line(62:end));
                
                
                line=fgetl(navFile);
                idot=str2double(line(5:23));
                CodesOnL2=str2double(line(24:42));
                week_toe=str2double(line(43:61));
                L2Pflag=str2double(line(62:end));
                
                
                line=fgetl(navFile);
                SVaccuracy=str2double(line(5:23));
                SVhealth=str2double(line(24:42));
                tgd=str2double(line(43:61));
                IODC=str2double(line(62:end));
                
                
                line=fgetl(navFile);
                transmissionTime=str2double(line(5:23));
                fitInterval=str2double(line(24:42));
                
                %Conversion to the format required by function
                %sat_coordinates_XYZ
                keplerArray(1)  = svprn;
                keplerArray(2)  = 0;%toc, assign it later
                keplerArray(3)  = af2;
                keplerArray(4)  = af1;
                keplerArray(5)  = af0;
                keplerArray(6)  = TODE;
                keplerArray(7)  = crs;
                keplerArray(8)  = deltan;
                keplerArray(9)  = M0;
                keplerArray(10)  = cuc;
                keplerArray(11)  = ecc;
                keplerArray(12)  = cus;
                keplerArray(13)  = roota;
                keplerArray(14)  = toe;
                keplerArray(15)  = cic;
                keplerArray(16)  = Omega0;
                keplerArray(17)  = cis;
                keplerArray(18)  = i0;
                keplerArray(19)  = crc;
                keplerArray(20)  = omega;
                keplerArray(21)  = Omegadot;
                keplerArray(22)  = idot;
                keplerArray(23)  = CodesOnL2;
                keplerArray(24)  = week_toe;
                keplerArray(25)  = L2Pflag;
                keplerArray(26)  = SVaccuracy;
                keplerArray(27)  = SVhealth;
                keplerArray(28)  = tgd;
                keplerArray(29)  = IODC;
                keplerArray(30)  = transmissionTime;
                
                if constellation == 'G'
                    % convert TOC (yy-mm-dd hh:mm:ss) to seconds into the GPS week
                    %   week=eph(count,24);
                    days=ymd2mjd(Time(1:3))-gps2mjd([week_toe 0]);
                    sec_in_week=days*86400+Time(4)*3600+Time(5)*60+Time(6);
                    keplerArray(2) = sec_in_week;
                    % set the global GPSWEEK (if necessary)
                    if length(GPSWEEK)==0
                        GPSWEEK=week_toe;
                    end
                    gpsEphemeris =  [gpsEphemeris keplerArray];
                elseif constellation == 'C'
                    
                    % convert TOC (yy-mm-dd hh:mm:ss) to seconds into the GPS week
                    %   week=eph(count,24);
                    deltat=14;
                    if week_toe<1356
                       week_toe=week_toe+1356;%北斗系统与GPS系统之间有1356周的时间差
%                        deltat=14;
                    end
                    keplerArray(24)  = week_toe;
                    days=ymd2mjd(Time(1:3))-gps2mjd([week_toe 0]);
                    sec_in_week=days*86400+Time(4)*3600+Time(5)*60+Time(6)+deltat;
                    keplerArray(2) = sec_in_week;
%                     keplerArray(14) = sec_in_week;
                    keplerArray(30) = sec_in_week;
                    % set the global GPSWEEK (if necessary)
                    if length(GPSWEEK)==0
                        GPSWEEK=week_toe;
                    end
                    beidouEphemeris =  [beidouEphemeris keplerArray];
                else
                    error 'Unknown constellation'
                    %Should never reach this point, as there is a case
                    %above.
                end
                
                
            otherwise
                %error 'Unknown constellation'
                
                
        end
        
    else
        error ('Wrong counting. New ephemeris expected.')
    end
    
end

% Construct output
% outputEphemeris.glonassEphemeris        = real(glonassEphemeris);
outputEphemeris.gpsEphemeris            = real(gpsEphemeris);
outputEphemeris.beidouEphemeris         = real(beidouEphemeris);
outputEphemeris.ionosphericParameters   = real(ionosphericParameters);
outputEphemeris.DOY                     = real(DOY);
outputEphemeris.leapSeconds             = real(leapSeconds);


fclose(navFile);
'Ephemeris loaded correctly'


end
