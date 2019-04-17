function [dbltime]=dbltimeadd(dbltime,decsec)
dbltime.dec=dbltime.dec+decsec;
while dbltime.dec<0
    dbltime.sec=dbltime.sec-1;
    dbltime.dec=dbltime.dec+1;
end
while dbltime.dec>1
    dbltime.sec=dbltime.sec+1;
    dbltime.dec=dbltime.dec-1;
end