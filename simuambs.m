% function simuambs
clear;clc;
option=struct('sys',[],'freqs',[],'stdcode',[],'stdphase',[],'Nsamp',[],'stdion',[],'tropo',[],'ldeg',[],'pdeg',[],'filename',[]);
opts=option;
option.Nsamp=1E5;

freqindex=cell(9,1);
freqindex{1,1}=[1 0 0; 0 0 0]'; freqindex{2,1}=[1 1 0; 0 0 0]'; freqindex{3,1}=[1 1 1; 0 0 0]';
freqindex{4,1}=[0 0 0; 1 0 0]'; freqindex{5,1}=[0 0 0; 1 1 0]'; freqindex{6,1}=[0 0 0; 1 1 1]';
freqindex{7,1}=[1 0 0; 1 0 0]'; freqindex{8,1}=[1 1 0; 1 1 0]'; freqindex{9,1}=[1 1 1; 1 1 1]';
% 
% freqinx=ones(3,2);[55S, 85S, 0N, 35N, 55N]£¬
cnt=0;
for ldeg=[115 140]
    option.ldeg=ldeg;                       % longitude of receiver [degrees]
% option.pdeg=35.0;  
    for pdeg=[30 50]
       option.pdeg=pdeg;   
        for find=1:9
            option.freqs=freqindex{find,1}==1;
            for stdpind=0.002:0.001:0.003
                option.stdphase=stdpind;
                for scale=[100 150]
                        option.stdcode=option.stdphase*scale;
                    for stdiind=[0.005,0.010,0.015,0.020,0.030]
                        option.stdion=stdiind;
                        if stdiind>=0.01
                          option.tropo='Tfloat';
                        else
                          option.tropo='Tfixed';
                        end
                        cnt=cnt+1;
                        option.filename=getfilename(option);
                        opts(cnt)=option;
%                         res=generate_Q_ahat(option);
%                         save(strcat('./simudata/',option.filename),'res');
                    end
                    
                end
                
            end
%             end
        end
    end
end
save('options.mat','opts');
% option.freqs=[1 0 0;1 1 0]'==1;%L1, L2, L5; B1, B2 B5;
% option.stdcode=0.3;
% option.stdphase=0.003;
% option.stdion=0.01;
% option.tropo='Tfloat';
% option.ldeg=115.0;                       % longitude of receiver [degrees]
% option.pdeg=35.0;  


