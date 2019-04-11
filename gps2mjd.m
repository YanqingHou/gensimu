function mjd=gps2mjd(gps)
%GPS2MJD GPS Week and Second past Sunday 0hrs to Modified Julian Date. 
%        GPS2MJD(gps) returns the Modified Julian Date MJD (JD-2400000.5) 
%        given the two-element row vector with the decimal GPS week number 
%        and the number of seconds past midnight of last saterday/sunday
%
%           gps = [ GPSweek GPSsecond]
%
%        If the GPSweek and Second are not OK the function returns the 
%        a NaN.

%        H. van der Marel, LGR, 29-04-95
%        (c) Geodetic Computing Centre, TU Delft

mjd = gps(:,1)*7 + gps(:,2)/86400 + 44244;
