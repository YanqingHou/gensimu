function [cursec,xsat,dts]=Alm2Pos(curweek,cursec,eph,sys,prn)
%��ʱ�Ȳ�����ʹ��α����м���
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
toe=eph(14); cic=eph(15); Omega0=eph(16); cis=eph(17); deltai=eph(18); crc=eph(19);
omega=eph(20); Omegadot=eph(21); idot=eph(22); week=eph(24); svh=eph(26);
tgd=eph(28); IODC=eph(29); tor=eph(30);
a=roota^2;

if strcmp(sys,'GPS')
    gm = 3.9860050E14;
    omge = 7.2921151467E-5;
    i0=0;%GPS��������eph��18������ľ���i0�����������������������������i=deltai+i0;��Ȼdeltai=eph��18������ô���ǽ�i0����Ϊ0���ɡ�

elseif strcmp(sys,'BDS')
    gm = 3.986004418E14;
    omge = 7.292115E-5;
    if prn<=265
        i0=0;
    else
        i0=0.3*pi;
    end
else
    %     gm
    %     omge
end
% ���������Ӳ����
% tΪ��ǰ������Ԫ�źŷ���ʱ�̵�ʱ�䣨�����룩������Ϊweekcur
% Dtsvstruct=cursec; Dtsvstruct.dec=0;
extrasec=604800*(curweek-week);%��ǰ��Ԫ����ʱ�䲥���Ĳο���Ԫ�Ƿ���������
cursec.sec=cursec.sec+extrasec;%��������ǰʱ���������

%���㵱ǰ��Ԫ���Ƿ���ʱ��ʱ�����ο���Ԫtoc��ʱ�䳤�ȡ�
% Dtsvstruct.sec=604800*(weekcur-week) - toc;
% dtstruct=timestructadd(t,Dtsvstruct);
dt=cursec.sec+cursec.dec-toc;
% dt = t - toc + 604800*(GPSWEEK-week);
for i=1:2
    dt=dt-(af0+af1*dt);
end
dt=af0+af1*dt;
cursec=dbltimeadd(cursec,-dt);
% t=t-dt;

% tkstruct=timestructadd(t,Dtsvstruct);
tkstruct=dbltimeadd(cursec,-toe);
tk=tkstruct.sec+tkstruct.dec;
% tk= t - toc + 604800*(GPSWEEK-week);

n = sqrt(gm/(a^3)) ;
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
i=deltai+i0;

x=r*cos(u); y=r*sin(u); cosi=cos(i);
%beidou
if(strcmp(sys,'BDS') && prn<=265)
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
dts=af0+af1*tk;
dts=dts-2*sqrt(gm*a)*ecc*sinE/c^2;
xsat=rs;
% 	// relativity correction
% 	*dts=*dts-2.0*sqrt(mu*ptEph->dSqrtA)*ptEph->dE*sinE/CLIGHT/CLIGHT;
%
% 	// position and clock error variance
% 	*var=var_uraeph(ptEph->iSva);
% [var]=var_uraeph(svh);