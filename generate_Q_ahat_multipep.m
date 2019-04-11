function res=generate_Q_ahat_multipep(option)
% this routine generates float ambiguity variance matrices for GPS, Galileo
% and integrated GPS+Galileo, as well as corresponding float ambiguity
% vectors using Monte Carlo simulations.
% The results are stored in the structure res (length equal to the number
% of epochs considered).
%
% res(k).QaS  : float ambiguity variance matrix, size [nxn], for epoch k
% res(k).ahatS: generated float ambiguity vectors for epoch k, this is a
%               matrix of size [n x Nsamp], where n is the number ambiguities,
%               and Nsamp the number of generated vectors
% res(k).PsS  : bootstrapped success rate for epoch k
%
% with S=i for GPS+Galileo, S=1 for GPS, S=2 for Galileo
% load('options.mat');
% option=opts(158);

freqs=option.freqs;
stdcode=option.stdcode;
stdphase=option.stdphase;
% Nsamp=option.Nsamp;
sdion=option.stdion;
tropo=option.tropo;
ldeg  =option.ldeg;                       % longitude of receiver [degrees]
pdeg  =option.pdeg;
% close all;
% clear all;
% clc
% warning off;
% freqs=[1 1;1 1; 1 1];
freqno1=sum(freqs(:,1));
freqno2=sum(freqs(:,2));

freq1 = [1575.42e6;1227.60e6; 1176.45e6];     % GPS frequencies
freq2 = [1561.098e6;1207.14e6; 1268.56e6];    % BDS frequencies

freq1=freq1(freqs(:,1));
freq2=freq2(freqs(:,2));

sigcode=stdcode*ones(freqno1+freqno2,1);
sigphase=stdphase*ones(freqno1+freqno2,1);
% EDIT THE PARAMETERS BELOW
% Nsamp = 1e3;                          % number of samples (simulated float ambiguity vectors)
% if Nsamp = 0: no simulations

% freq1 = [1575.42e6;1227.6e6];       % GPS frequencies
% freq2 = [1575.42e6;1176.45e6];      % Galileo frequencies


% freq1 = [1575.42e6;1227.6e6;1176.45e6];% freq2 = [1575.42e6;1176.45e6];
% freq2 = [1575.42e6;1207.14e6;1176.45e6];


% sigcode   = [0.2;0.2;0.2;0.2 ];       % undifferenced code standard deviations per frequency and per system
%                                       % first for all GPS frequencies, then
%                                       % for Galileo
% sigphase  = [0.002;0.002;0.002;0.002];% undifferenced phase standard deviations per frequency and per system
cfix      = 0;                        % cfix=1: coordinates fixed (both receivers);
% sdion     = 0.025;                    % undifferenced ionospheric standard deviation [meter]
% tropo     = 'Tfloat';                 % 'Tfixed': ZTD not estimated, 'Tfloat': ZTD estimated

% ldeg  = 115.35;                       % longitude of receiver [degrees]
% pdeg  = -33.3;                        % latitude of receiver [degrees]
%20 days
starttime  = '22-nov-2013 0:00';      % first epoch for which to generate data
endtime    = '28-nov-2013 23:59';     % last epoch for which to generate data

mixeph =      rd_alma('Almanac.alm');
eph   =      mixeph(mixeph(:,1)<33,:);
eph2   =     mixeph(mixeph(:,1)>260,:);


% starttime  = '21-mar-2012 0:00';      % first epoch for which to generate data
% endtime    = '31-mar-2012 24:00';     % last epoch for which to generate data
%
%
% eph        = rdyuma('yumaGPS20120319.txt');  % GPS almanac (use almanac corresponding to the date considered)
% eph2       = rdyuma('yumaGAL20120319.txt');  % Galileo almanac
int        = 1800;                     % interval [sec] for which model will be generated

% no_epochs     = 1;                    % number of epochs used in resolving float solution
cutoff        = 10;                   % cutoff elevation
% FROM HERE ON NO EDITING REQUIRED

prad          = pdeg*pi/180;
lrad          = ldeg*pi/180;

plh           = [prad lrad 0];  % vector [latitudes longitudes heights]
xyz           = plh2xyzwgs(plh);

tsat          = mktsat ( starttime,endtime,int);

[xs,ys,zs]    = cpaziele (tsat,eph);
[xs2,ys2,zs2] = cpaziele (tsat,eph2);
eps=50;
%
lt = size(tsat,2);
% time = 1:lt;
% time=47;
res = struct('Qa',[],'ahat',[],'Ps',[],'Qb',[],'Qab',[]);
if freqno1==0
    for time=1:lt
        for no_epochs=1:eps
            [Qa,Ps,Qb,Qab] = SingleQgg(freq2,sigcode,sigphase,sdion,tropo,no_epochs,cutoff,xs2(:,time),ys2(:,time),zs2(:,time),xyz,plh,cfix);
            res(time,no_epochs).Qa = Qa; % float ambiguity variance matrix for Galileo
            res(time,no_epochs).Qb = Qb(1:4,1:4);
            res(time,no_epochs).Qab = Qab(:,1:4);
            res(time,no_epochs).Ps= Ps;
            %         if Nsamp > 0
            %             res(k).ahat  = mvnrnd(zeros(1,size(Qa,1)),Qa,Nsamp)';
            %         end
        end
    end
elseif freqno2==0
    for time=1:lt
        for no_epochs=1:eps
            [Qa,Ps,Qb,Qab] = SingleQgg(freq1,sigcode,sigphase,sdion,tropo,no_epochs,cutoff,xs(:,time),ys(:,time),zs(:,time),xyz,plh,cfix);
            res(time,no_epochs).Qa = Qa; % float ambiguity variance matrix for Galileo
            res(time,no_epochs).Qb = Qb(1:4,1:4);
            res(time,no_epochs).Qab = Qab(:,1:4);
            res(time,no_epochs).Ps= Ps;
            %         if Nsamp > 0
            %             res(k).ahat  = mvnrnd(zeros(1,size(Qa,1)),Qa,Nsamp)';
            %         end
        end
    end
else
    for time=1:lt
        for no_epochs=1:eps
            [Qa,Ps,Qb,Qab]  = DualQgg(freq1,freq2,...
                sigcode,sigphase,sdion,tropo,no_epochs,cutoff,...
                xs(:,time),ys(:,time),zs(:,time),xs2(:,time),ys2(:,time),zs2(:,time),xyz,plh,cfix);
            res(time,no_epochs).Qa = Qa; % float ambiguity variance matrix for Galileo
            res(time,no_epochs).Qb = Qb(1:4,1:4);
            res(time,no_epochs).Qab = Qab(:,1:4);
            res(time,no_epochs).Ps= Ps;
            %             if Nsamp > 0
            %                 res(k).ahat  = mvnrnd(zeros(1,size(Qa,1)),Qa,Nsamp)';
            %             end
            
        end
    end
end
