function gps=mjd2gps(mjd)
%MJD2GPS Modified Julian Date to GPS Week and Second past Sunday 0hrs. 
%        MJD2YMD(mjd) returns the two-element row vector with the decimal
%        GPS week number and the number of seconds past midnight of last 
%        saturday/sunday
%
%           gps = [ GPSweek GPSsecond]
%
%        Input is the Modified Julian Date MJD (JD-2400000.5). If the 
%        Modified Julian Date is not OK the function returns a row vector 
%        of NaN's.

%        H. van der Marel, LGR, 29-04-95
%        (c) Geodetic Computing Centre, TU Delft

% Number of days since starting epoch of GPS weeks (Sunday 06-Jan-1980)
mjds = mjd - 44244;
      
% Week number, day and second in week
nweek  = floor(mjds/7);
day    = mjds - nweek*7;
second = day*86400;
gps = [ nweek second ];
