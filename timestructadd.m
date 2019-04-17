function [timestruct]=timestructadd(time1,time2)
timestruct=time1;
timestruct.sec=time1.sec+time2.sec;
timestruct.dec=time1.dec+time2.dec;
timestruct=dbltimeadd(timestruct,0);