function ymd=mjd2ymd(mjd)
%MJD2YMD Modified Julian Date to Gregorian Calendar.
%        MJD2YMD(mjd) returns the three-element row vector with the decimal
%        year, month and day in the Gregorian Calendar
%
%           ymd = [ year month day ]
%
%        Input is the Modified Julian Date MJD (JD-2400000.5) for 0 hrs. 
%        If the Modified Julian Date is not OK the function returns the 
%        a row vector of NaN's.
%
%        See also YMD2MJD, YMD2STR and STR2YMD.

%        H. van der Marel, LGR, 29-04-95
%        (c) Geodetic Computing Centre, TU Delft

% Check if Modified Julian Date is acceptable, else replace by a NaN
i=find(mjd<=-2395522 | mjd>=1000000000);mjd(i)=NaN*i;

% Express day in Gregorian calendar
jd=fix(mjd)+2400001;

n4  = 4 * ( jd + fix((fix((6*fix((4*jd-17918)/146097))/4)+1)/2) - 37 );
nd10=10 * fix(rem(n4-237,1461)/4) + 5;

year=fix(n4/1461)-4712;
month=rem(fix(nd10/306)+2,12)+1;
day=fix(rem(nd10,306)/10)+1;

day=day+mjd-fix(mjd);

ymd=[year month day];
