function runonce_mult(num)
load('options.mat');
res=generate_Q_ahat_multipep(opts(num));
save(strcat('./simudata/Mult_',opts(num).filename),'res','-v7.3');