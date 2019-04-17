function [bdsweek,secs]=utc2bdt(str)
% utc0=[2006,1,1,0,0,0];


% function [gpsweek,secs] = utc2gpst (str);
%UTC2GPST: Convert time/date string to GPSWEEK / Seconds
%
% This routine converts a time and date, given as a string in the form of
% 'dd-mon-yyyy hh:mm:ss.sss' into the gpsweek and seconds into this
% week. If the global variable GPSWEEK (uppercase) is defined, the seconds
% will be relative to this week.
%
% Syntax:
%    [gpsweek,secs] = utc2gpst (str);
%
% Input arguments:
%    str    : Date/time string in the form: 'dd-mon-yyyy hh:mm:ss.sss'
%
% Output arguments:
%    gpsweek: gpsweek number
%    secs   : Seconds relative to this week
%

% ----------------------------------------------------------------------
% File.....: utc2gpst.m
% Date.....: 27-MAY-1999
% Version..: 1.0
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

% ---------------------------------------
% --- Declare global variable GPSWEEK ---
% ---------------------------------------

global BDSWEEK;

% ----------------------
% --- Convert string ---
% ----------------------

[y, m, d, hh, mm, ss] = datevec (str);
bdst    = mjd2bdt(ymd2mjd([y m d]));
bdst(2) = bdst(2) + 3600*hh + 60*mm + ss;

bdsweek = bdst(1);
secs    = bdst(2);

% -------------------------------------------------------------
% --- if GPSWEEK was defined, adapt the results accordingly ---
% -------------------------------------------------------------

if length(BDSWEEK) ~= 0

  secs    = secs + (bdsweek - BDSWEEK) * 7*24*3600;
  bdsweek = GPSWEEK;

end


end
function bdt=mjd2bdt(mjd)
%MJD2BDS Modified Julian Date to BDS Week and Second past Sunday 0hrs. 
%        MJD2YMD(mjd) returns the two-element row vector with the decimal
%        BDS week number and the number of seconds past midnight of last 
%        saturday/sunday
%
%           bdt = [ BDSweek BDSsecond]
%
%        Input is the Modified Julian Date MJD (JD-2400000.5). If the 
%        Modified Julian Date is not OK the function returns a row vector 
%        of NaN's.

%        H. van der Marel, LGR, 29-04-95
%        (c) Geodetic Computing Centre, TU Delft

% Number of days since starting epoch of BDS weeks (01-Jan-2006)
mjds = mjd - 53736;
      
% Week number, day and second in week
nweek  = floor(mjds/7);
day    = mjds - nweek*7;
second = day*86400;
bdt = [ nweek second ];
end