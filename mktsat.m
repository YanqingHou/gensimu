function [strweek,tsat] = mktsat (startt,endt,interval)
%MKTSAT: Create an array of epochs for which positions are to be computed
%
% This routine creates an array of epochs for which satellite positions
% (and azimuths, elevations, DOPS, etc) are to be computed. An epochs is
% specified as a number of seconds into the GPSWEEK (as given in the global
% variable GPSWEEK). Warnings are issued for "unlikely situations"
% 
% Syntax:
%    [tsat] = mktsat (start,end,interval);
%
% Input arguments:
%    start  : Start time, in the form 'dd-mon-yyyy hh:mm:ss.sss'
%    end    : End time, in the form 'dd-mon-yyyy hh:mm:ss.sss'
%             Note that the month should be specified as the first three
%             characters of the name of the month, so september 4th 1999
%             will be "04-SEP-1999"
%    interval: Interval between epochs, in seconds
%
% Output arguments:
%    tsat: Array with epochs (seconds into GPSWEEK)

% ----------------------------------------------------------------------
% File.....: mktsat.m
% Date.....: 25-MAY-1999
% Version..: 1.0
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

% -------------------------------------------------
% --- Declare and check global variable GPSWEEK ---
% -------------------------------------------------

% global GPSWEEK;
% 
% if length(GPSWEEK) == 0
%   error ('GPSWEEK not set, exiting!');
%   return;
% end;

% --------------------------------------------------
% --- Convert UTC times/dates to gpsweek/seconds ---
% --------------------------------------------------

[strweek,strsecs] = utc2gpst (startt);
[endweek,endsecs] = utc2gpst (endt);

% -----------------------------------------------------------------
% --- Check if weeks are not too different from current GPSWEEK ---
% -----------------------------------------------------------------

% if abs (strsecs) > 31*24*3600;
%   disp ('WARNING: Starting time and ephemeris differ at least a month');
% end;
% 
% if abs (endsecs) > 31*24*3600;
%   disp ('WARNING: End time and ephemeris differ at least a month');
% end;

% --------------------------------------
% --- Convert to the current strweek ---
% --------------------------------------

% strsecs = strsecs + (strweek - GPSWEEK) * 7*24*3600;
endsecs = endsecs + (endweek - strweek) * 7*24*3600;

% ----------------------------
% --- Reasonable time-span ---
% ----------------------------

if strsecs > endsecs
  disp ('Warning: End time is before start time, empty tsat returned');
  tsat = [];
  return;
end

if endsecs-strsecs > 2*24*3600
  disp ('Warning: Time-span covers more than 2 days!');
end

% --------------------------
% --- Create tsat-vector ---
% --------------------------

tsat = strsecs:interval:endsecs;

% --------------------------------------
% --- Warning if result is very long ---
% --------------------------------------

if length (tsat) > 1000
  disp ('Warning: over 1000 epochs selected, might take a while ...');
end
 
% ------------------------------
% --- End of function mktsat ---
% ------------------------------
