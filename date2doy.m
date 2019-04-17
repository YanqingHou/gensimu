function doy=date2doy(date)
%%convert date to doy
% epoch=time2epoch(date);
% epoch(2)=1; epoch(3)=1;
% epoch(4)=0; epoch(5)=0; epoch(6)=0;
% 
% doy=timediff(date,epoch2time(epoch))/86400.0+1;

%% method 2
%   lineSplit = strsplit(line);
% yr=date(1);
doy = datenum(date)-datenum([date(1),1,1])+1;
% doy = today-datenum(['1-Jan-' year(today)])+1;