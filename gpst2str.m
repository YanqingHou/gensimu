function [str] = gpst2str (secs,gpsweek);
%GPST2STR: Convert gpsweek/seconds into a string
%
% This routine converts a gpsweek/seconds combination into a string 
% of the form 'dd-mon-yyyy hh:mm:ss.sss'. The time itself is not
% transformed, the string will still represent GPS-time.
%
% Syntax:
%    [str] = gpst2str (secs,gpsweek);
%
% Input arguments:
%    secs   : Seconds relative to this week
%    gpsweek: gpsweek number (optional, if not specified, global
%             variable GPSWEEK will be used);
%
% Output arguments:
%    str    : Date/time string in the form: 'dd-mon-yyyy hh:mm:ss.sss'
%

% ----------------------------------------------------------------------
% File.....: gpst2utc.m
% Date.....: 27-MAY-1999
% Version..: 1.0
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

% -----------------------------------------------
% --- Declare and use global variable GPSWEEK ---
% -----------------------------------------------

global GPSWEEK;

if nargin < 2; 

  if length(GPSWEEK) == 0;
    error ('GPSWEEK was not specified');
    return;
  end;
  
  gpsweek = GPSWEEK;

end;

% -----------------------------------
% --- Convert date/time to string ---
% -----------------------------------

months = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];

gpsweek = gpsweek + floor (secs/604800);
secs    = mod (secs,604800);
ymd     = mjd2ymd (gpsweek*7 + secs/86400 + 44244);
secs    = rem(secs,86400);
hh      = floor (secs/3600);
mm      = floor (secs/60);
mm      = mm - hh*60;
ss      = rem (secs,60);

str = sprintf ('%2d-xxx-%4d %2.2d:%2.2d:%2.2d.%3.3d', ...
	       floor(ymd(3)),ymd(1),hh,mm,floor(ss),round(1000*(rem(ss,1))));

str(4:6) = months(ymd(2),1:3);

% --------------------------------
% --- End of function utc2gpst ---
% --------------------------------
