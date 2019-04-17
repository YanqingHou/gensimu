% function simuambs
clear;clc;
option=struct('sys',[],'freqs',[],'stdcode',[],'stdphase',[],'Nsamp',[],'stdion',[],'tropo',[],'ldeg',[],'pdeg',[],'filename',[]);
option.Nsamp=10;

freqindex=cell(9,1);
freqindex{1,1}=[1 0 0; 0 0 0]'; freqindex{2,1}=[1 1 0; 0 0 0]'; freqindex{3,1}=[1 1 1; 0 0 0]';
freqindex{4,1}=[0 0 0; 1 0 0]'; freqindex{5,1}=[0 0 0; 1 1 0]'; freqindex{6,1}=[0 0 0; 1 1 1]';
freqindex{7,1}=[1 0 0; 1 0 0]'; freqindex{8,1}=[1 1 0; 1 1 0]'; freqindex{9,1}=[1 1 1; 1 1 1]';
% 经度-纬度
locs=[115+1/3,18+1/3;%南海
    123+56/60,29+56/60;%东海
    113+15/60,23+8/60;%广州
    108+56/60,34+20/60];%西安
freqindexs=[0 0 0; 1 1 0]';%bds dual-freq
stdiind=0.005;
c2pscale=100;
% freqinx=ones(3,2);[55S, 85S, 0N, 35N, 55N]，
cnt=0;
option.stdion=stdiind;
option.tropo='Tfixed';
for i=4
    option.ldeg=locs(i,1);%经度
    option.pdeg=locs(i,2);%纬度
    option.freqs=freqindexs==1;
    %             for stdcind=0.15:0.05:0.40
    %             option.stdcode=stdcind;
    for stdpind=0.001:0.0005:0.005
        option.stdphase=stdpind;
        option.stdcode=stdpind*c2pscale;
        cnt=cnt+1;
        option.filename=getfilename(option,'Haige');
        res=generate_Q_ahat(option);
        save(strcat('./Haigesimudata/',option.filename),'res');
    end
    %                     end    %             end
end
% end
cnt