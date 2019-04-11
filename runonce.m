function runonce(num)
load('options.mat');
res=generate_Q_ahat(opts(num));
save(strcat('./simudata/',opts(num).filename),'res','-v7.3');