function [gpsweek,secs] = utc2gpst (str);
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

global GPSWEEK;

% ----------------------
% --- Convert string ---
% ----------------------

[y m d hh mm ss] = datevec (str);
gpst    = mjd2gps(ymd2mjd([y m d]));
gpst(2) = gpst(2) + 3600*hh + 60*mm + ss;

gpsweek = gpst(1);
secs    = gpst(2);

% -------------------------------------------------------------
% --- if GPSWEEK was defined, adapt the results accordingly ---
% -------------------------------------------------------------

if length(GPSWEEK) ~= 0;

  secs    = secs + (gpsweek - GPSWEEK) * 7*24*3600;
  gpsweek = GPSWEEK;

end ;

% --------------------------------
% --- End of function utc2gpst ---
% --------------------------------



