function [cursec,xsat,dts]=Eph2Pos(curweek,cursec,eph,sys,prn,cordform)
%暂时先不考虑使用伪距进行计算
% global GPSWEEK;
% omgedo = 7292115.1467e-11;     % angular velocity of Earth
gm     = 3.986004418E14;        % gravitational constant
c      = 299792458;            % speed of light
epsec  = 1e-12;                % stop criterion for iterations
maxit  = 10;                   % maximum iterations for E
SIN_5=-0.0871557427476582;
COS_5=0.9961946980917456;

toc=eph(2);  af0=eph(3); af1=eph(4);af2=eph(5); TODE=eph(6); crs=eph(7);
deltan=eph(8); M0=eph(9); cuc=eph(10); ecc=eph(11); cus=eph(12); roota=eph(13);
toe=eph(14); cic=eph(15); Omega0=eph(16); cis=eph(17); i0=eph(18); crc=eph(19);
omega=eph(20); Omegadot=eph(21); idot=eph(22); week=eph(24); svh=eph(26);
tgd=eph(28); IODC=eph(29); tor=eph(30);
a=roota^2;

if strcmp(sys,'GPS')
    gm = 3.9860050E14;
    omge = 7.2921151467E-5;
elseif strcmp(sys,'BDS')
    gm = 3.986004418E14;
    omge = 7.292115E-5;
else
    %     gm
    %     omge
end
% 计算卫星钟差参数
% t为当前测量历元信号发射时刻的时间（周内秒），周数为weekcur
% Dtsvstruct=cursec; Dtsvstruct.dec=0;
extrasec=604800*(curweek-week);%当前历元距离时间播发的参考历元是否经历了数周
cursec.sec=cursec.sec+extrasec;%补偿到当前时刻秒计数上

%计算当前历元卫星发射时刻时间距离参考历元toc的时间长度。
% Dtsvstruct.sec=604800*(weekcur-week) - toc;
% dtstruct=timestructadd(t,Dtsvstruct);
dt=cursec.sec+cursec.dec-toc;
% dt = t - toc + 604800*(GPSWEEK-week);
for i=1:2
    dt=dt-(af0+af1*dt+af2*dt^2);
end
dt=af0+af1*dt+af2*dt^2;
cursec=dbltimeadd(cursec,-dt);
% t=t-dt;

% tkstruct=timestructadd(t,Dtsvstruct);
tkstruct=dbltimeadd(cursec,-toe);
tk=tkstruct.sec+tkstruct.dec;
% tk= t - toc + 604800*(GPSWEEK-week);

n = sqrt(gm/(a^3)) + deltan;
m = M0 + n*tk;
eccan = m;

for iter=1:maxit
    oldan = eccan;
    eccan = m + ecc*sin(oldan);
    if abs(oldan-eccan)<epsec, break, end
end
sinE=sin(eccan);  cosE=cos(eccan);

u=atan2(sqrt(1-ecc^2)*sinE,cosE-ecc)+omega;
r=a*(1-ecc*cosE);
i=i0+idot*tk;
sin2u=sin(2*u); cos2u=cos(2*u);
u=u+cus*sin2u+cuc*cos2u;
r=r+crs*sin2u+crc*cos2u;
i=i+cis*sin2u+cic*cos2u;
x=r*cos(u); y=r*sin(u); cosi=cos(i);
%beidou
if(strcmp(sys,'BDS') && prn<=5)
    %     omge=7.292115E-5;
    
    O=Omega0+Omegadot*tk-omge*toe;
    
    sinO=sin(O); cosO=cos(O);
    xg=x*cosO-y*cosi*sinO;
    yg=x*sinO+y*cosi*cosO;
    zg=y*sin(i);
    sino=sin(omge*tk); coso=cos(omge*tk);
    rs(1)= xg*coso+yg*sino*COS_5+zg*sino*SIN_5;
    rs(2)=-xg*sino+yg*coso*COS_5+zg*coso*SIN_5;
    rs(3)=-yg*SIN_5+zg*COS_5;
else
    %     omge=7.2921151467E-5 ;
    %     omge=7.292115E-5;
    
    O=Omega0+(Omegadot-omge)*tk-omge*toe;
    sinO=sin(O); cosO=cos(O);
    rs(1)=x*cosO-y*cosi*sinO;
    rs(2)=x*sinO+y*cosi*cosO;
    rs(3)=y*sin(i);
end
% tkstruct=timestructadd(t,Dtsvstruct);
% tk=tkstruct.sec+tkstruct.dec;
% tk = t - toc + 604800*(GPSWEEK-week);
dts=af0+af1*tk+af2*tk^2;
dts=dts-2*sqrt(gm*a)*ecc*sinE/c^2;
% xsat=rs;
% 以curweek周0时0分0秒为时间基准，建立惯性直角坐标系。
duration=cursec.sec+cursec.dec;
Omega=-omge*duration;
sinOmega=sin(Omega);cosOmega=cos(Omega);

Tr=[cosOmega, sinOmega, 0;
    -sinOmega, cosOmega, 0;
    0,          0,        1];

if strcmp(cordform,'ECEF')
    Trx=eye(3);
elseif strcmp(cordform,'Inertial')
    Trx=Tr;
else
    error('unknown coordinate form');
end
xsat=Trx*rs';

% 	// relativity correction
% 	*dts=*dts-2.0*sqrt(mu*ptEph->dSqrtA)*ptEph->dE*sinE/CLIGHT/CLIGHT;
%
% 	// position and clock error variance
% 	*var=var_uraeph(ptEph->iSva);
% [var]=var_uraeph(svh);