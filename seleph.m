function idx = seleph(prn,curweek,t,eph)
%SELEPH  Select satellite ephemeris
%        SELEPH selects the correct satellite ephemeris from eph and returns
%        the row number idx.
%        Syntax:
%                idx = sateph(prn,t,eph)
%        with
%                prn      satellite prn numbers 
%                t        time in seconds into refweek 
%                eph      matrix with ephemeris parameters from RXNAV
%                -------
%                idx      row number of the most recent ephemeris
% 
%        The satellite ephemeris data is stored in eph, idx contains
%        the row number(s) with ephemeris corresponding to the satellite(s)
%        and time(s) given in prn and t.

%        H. van der Marel, LGR, 29-04-95
%        (c) Geodetic Computing Centre, TU Delft


%
% Broadcast Satellite Ephemeris data stored in Eph(_,:)
%
%            name      description
%
% Eph(_, 1)  prn       Satellite prn number
% Eph(_, 2)  toc       Time of clock (seconds into GPS week)
% Eph(_,24)  week      GPS week # (for TOE and TOC)

%
% Declare global GPSWEEK
%

% global GPSWEEK;

%
% Loop through all ephemeris sets and select most recent one 
%

% [n,m] = size(eph);
% dtmin=inf;
% idx=-1;
% for count=1:n
%   if prn==eph(count, 1) 
%      dt=abs( t - eph(count,2) + 604800*(curweek-eph(count,24)));
%      if dt<dtmin
%         idx=count;
%         dtmin=dt;
%      end
%   end
% end
% 条件是：eph中只有prn卫星的星历。即eph第一列数值都一样
if length(unique(eph(:, 1)))==1 && prn==eph(1, 1) 
dt=abs( t - eph(:,2) + 604800*(curweek-eph(:,24)));
[dtmin,idx]=min(dt);
else
    error('eph contains too much prns');
end
