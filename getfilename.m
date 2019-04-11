function filename=getfilename(option)
txt={'L1','L2','L5';'B1','B2','B3'}';
if sum(option.freqs(:,1))>0 && sum(option.freqs(:,2))==0
   option.sys='GPS';
elseif sum(option.freqs(:,1))==0 && sum(option.freqs(:,2))>0
   option.sys='BDS';
elseif sum(option.freqs(:,1))>0 && sum(option.freqs(:,2))>0
   option.sys='GPSBDS';
else
   option.sys='NULL';
end
text=[]; txt1=txt(option.freqs);
for itxt=1:length(txt1)
    text=strcat(text,txt1(itxt));
end
stdp=option.stdphase/0.001;
stdc=option.stdcode/0.01;
stdi=option.stdion/0.001;
%filename: ÆµÂÊ+µØµã+Tropmodel+stdphase+stdcode+stdion£»
filename=strcat(option.sys,'_',text{1},'_E',num2str(option.ldeg),'N',num2str(option.pdeg),...
    '_',option.tropo,'_sdp',num2str(stdp),'mm_sdc',num2str(stdc),'cm_sdi',num2str(stdi),'mm');