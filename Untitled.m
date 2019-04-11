% read the datafile
load('options.mat');
existfile=cell(864,1);
fp=fopen('files.txt');
cnt=0;
while ~feof(fp)
    cnt=cnt+1;
    existfile{cnt,1}=fgets(fp);
end
fclose(fp);
existfile=deblank(existfile);
% exist=load('files.txt');
noexist=[];
for i=1:length(opts)
    sumflag=0;
    for j=1:cnt
        sumflag=sumflag+strcmp(strcat(opts(i).filename,'.mat'),existfile{j,1});
    end
    if ~sumflag
        noexist=[noexist; i];
    end
end
save('simurest.mat','noexist');
for i=1:length(noexist)
    runonce(noexist(i));
end
