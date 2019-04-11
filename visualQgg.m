function [Qa,Qa1,Qa2,Ps,Ps1,Ps2,Qb,Qb1,Qb2,Qab,Qab1,Qab2] = visualQgg(freq1,freq2,...
    sigcode,sigphase,sdion,tropo,no_epochs,cutoff,xs,ys,zs,xs2,ys2,zs2,xyz,plh,cfix)


freq  = [freq1;freq2];

model         = 'baseline';         % baseline
scenario      = 'rr';               % stationary
rectype       = 'non-cc';           % non-crosscorrelating

ncode         = length(freq1);      % number of code observations 
nphase        = ncode;              % number of phase observations
ncode2        = length(freq2);      % number of code observations 
nphase2       = ncode2;             % number of phase observations

nn      = 1;


% wavelengths
lambda = 299792458 ./ freq;
if strcmp(rectype,'cc')
    lambda(2) = 299792458 / (freq(1)-freq(2));
end
LAMBDA = diag(lambda);



% single difference ionospheric vc-matrix
sd    = 2;
sig   = [sigcode(1:ncode);sigphase(1:ncode);sigcode(ncode+1:end);sigphase(ncode+1:end)];
mu    = ([lambda(1:ncode)]./ lambda(1)).^2;
mu    = [mu(1:ncode);-mu(1:nphase)];
s2    = sd*sdion^2;

CI    = mu * mu';


mu2   = (lambda(ncode+1:end)./ lambda(ncode+1)).^2  ;
mu2   = [mu2(1:ncode);-mu2(1:nphase)];
CI2   = mu2 * mu2';


el       = cpelev1(xs,ys,zs,plh,xyz);
idx1     = find (el > cutoff);
no_sv    = length(idx1);
refsat   = find(el(idx1)==max(el(idx1)));
xsv      = xs(idx1);
ysv      = ys(idx1);
zsv      = zs(idx1);
el       = el(idx1);
zen      = deg2rad(90 - el)';

el2       = cpelev1(xs2,ys2,zs2,plh,xyz);
idx2      = find (el2 > cutoff);
no_sv2    = length(idx2);
refsat2   = find(el2(idx2)==max(el2(idx2)));
xsv2      = xs2(idx2);
ysv2      = ys2(idx2);
zsv2      = zs2(idx2);
el2       = el2(idx2);
zen2      = deg2rad(90 - el2)';



if ~cfix
    LOS  = cplos([xsv ysv zsv],xyz);   % LOS-vectors
    LOS2 = cplos([xsv2 ysv2 zsv2],xyz);
else 
    LOS  = [];
    LOS2 = [];
end
m2  = no_sv2 - 1;
m  = no_sv - 1;


Im   = eye(m);
em   = ones(m,1);
if strcmp(model,'single')
    D = Im;
else
    D = [Im(:,1:refsat-1) -em Im(:,refsat:m)];
end

if strcmp(tropo,'Tfixed')
    Mz = [];
else
    Mz   = tropomap(zen,plh(3),'s');    % tropo. wet delay mapping function
end

Mu = [LOS Mz];
if ~strcmp(scenario,'gf') && ~isempty(Mu)
    M1 = D * Mu;
elseif strcmp(scenario,'gf')
    M1 = Im;
elseif isempty(Mu)
    M1 = [];
end
M = []; % see eq.(19)
n1 = nphase + ncode;
n2 = 1;
Qy = [];
while (n1 > 0 )
    M = [M ; M1];
    n1 = n1 - 1;
    %elevation dependent weighting as in eq.(2.13) and (2.15)
    Wi = sd * (sig(n2) * diag(1 + 10 * exp(-el/10))).^2;
    Qy = blkdiag(Qy,D*Wi*D'); %first term after =-sign in eq.(2.20)/(2.21)
    n2 = n2 + 1;
end
if sdion > 0
    Wi = sd * (sdion * diag(1 + 10 * exp(-el/10))).^2;
    QI = kron(CI,D*Wi*D');       % 2nd term after =-sign in eq.(2.20)/(2.21)
else
    QI = zeros(size(Qy));
end
iQy = inv(Qy+QI);               % inverse of matrix in eq.(2.20)/(2.21)
N =[];  % see eq.(19)
if nphase>0
    N1 = zeros(ncode*m);
    N2 = kron(LAMBDA(1:nphase,1:nphase),Im);
    N = [N1;N2];
end
NQN = N'*iQy * N;
% MQN = M' * iQy * N;
% MQM = M' * iQy * M;
% iMQM = inv(MQM);
% also for GALILEO if system=GPS + GALILEO

Im2   = eye(m2);
em2   = ones(m2,1);
if strcmp(model,'single')
    D2 = Im2;
else
    D2 = [Im2(:,1:refsat2-1) -em2 Im2(:,refsat2:m2)];
end

if strcmp(tropo,'Tfixed')
    Mz2 = [];
else
    Mz2   = tropomap(zen2,plh(3),'I');    % tropo. wet delay mapping function
end


Mu2 = [LOS2 Mz2];
if ~strcmp(scenario,'gf') && ~isempty(Mu)
    M21 = D2 * Mu2;
elseif strcmp(scenario,'gf')
    M21 = Im2;
elseif isempty(Mu)
    M21 = [];
end


M2 = [];   % see eq.(19)
%%%
n1 = nphase + ncode;
Qy2 = [];
while (n1 > 0 )
    M2 = [M2 ; M21];
    n1 = n1 - 1;
    %elevation dependent weighting as in eq.(2.13) and (2.15)
    Wi = sd * (sig(n2) * diag(1 + 10 * exp(-el2/10))).^2;
    Qy2 = blkdiag(Qy2,D2*Wi*D2'); %first term after =-sign in eq.(2.20)/(2.21)
    n2 = n2 + 1;
end
if sdion > 0
    Wi = sd * (sdion * diag(1 + 10 * exp(-el2/10))).^2;
    QI2 = kron(CI2,D2*Wi*D2');       % 2nd term after =-sign in eq.(2.20)/(2.21)
else
    QI2 = zeros(size(Qy2));
end
iQy2 = inv(Qy2+QI2);               % inverse of matrix in eq.(2.20)/(2.21)

%%%

N2 =[];   % see eq.(19)
if nphase2>0
    N21 = zeros(ncode2*m2);
    N22 = kron(LAMBDA(nphase+1:end,nphase+1:end),Im2);
    N2 = [N21;N22];
end
NQN2 = N2' * iQy2 * N2;
% MQN2 = M2' * iQy2 * N2;
% MQM2 = M2' * iQy2 * M2;

if isempty(M)
    
    ins1 = NQN;
    ins2 = NQN2;
    ins  = blkdiag(NQN,NQN2);
    
    MQM = []; MQN = []; MQM2 = []; MQN2 = [];
else
    MQN2 = M2' * iQy2 * N2;
    MQM2 = M2' * iQy2 * M2;  
    MQN = M' * iQy * N;
    MQM = M' * iQy * M;
    iMQM = inv(MQM);
    
    ins1   = NQN - MQN' * iMQM * MQN;
    ins2   = NQN2 - MQN2' * inv(MQM2) * MQN2;
    
    iMQM = inv(MQM + MQM2);
    
    ins   = NQN - MQN' * iMQM * MQN;
    
    ins12 = -MQN' * iMQM * MQN2;
    ins22  = NQN2 - MQN2' * iMQM * MQN2;
    ins = [ins ins12;ins12' ins22];
end

Qa1 = (1/no_epochs)*inv(ins1); % eq.(58)
Qa2 = (1/no_epochs)*inv(ins2); % eq.(58)
Qa = (1/no_epochs)*inv(ins); % eq.(58)

Qa1   =(tril(Qa1,0)+tril(Qa1,-1)');
Qa2   =(tril(Qa2,0)+tril(Qa2,-1)');
Qa   =(tril(Qa,0)+tril(Qa,-1)');

Ps1 = cpsucrate(Qa1,'LS');
Ps2 = cpsucrate(Qa2,'LS');
Ps  = cpsucrate(Qa,'LS');
% ---------------------
% --- success-rates ---
% ---------------------
% 

    
if strcmp(scenario,'rr')
    Ik  = eye(no_epochs);
    ek  = ones(no_epochs,1);
    
    Qab1 = -Qa1 * kron(ek', MQN'/MQM );      % eq.(3.24)
    MQMk = kron(Ik,MQM);
    Qb1  = inv(MQMk) - MQMk \ kron(ek, MQN) * Qab1; %eq.(3.25)
    
    Qab2 = -Qa2 * kron(ek', MQN2'/MQM2 );      % eq.(3.24)
    MQM2k = kron(Ik,MQM2);
    Qb2  = inv(MQM2k) - MQM2k \ kron(ek, MQN2) * Qab2; %eq.(3.25)
else
    
    Qab1 = -Qa1 * MQN'/MQM ;        % eq.(60)
    Qb1  = inv(MQM) - MQM\MQN * Qab1;  % eq.(61)
    Qab2 = -Qa2 * MQN2'/MQM2 ;        % eq.(60)
    Qb2  = inv(MQM2) - MQM2\ MQN2 * Qab2;  % eq.(61)
    
end

MQM = MQM + MQM2;
MQN = [MQN MQN2];
if strcmp(scenario,'rr')
    Qab = -Qa * kron(ek', MQN'/MQM );      % eq.(3.24)
    MQMk = kron(Ik,MQM);
    Qb  = inv(MQMk) - MQMk \ kron(ek, MQN) * Qab; %eq.(3.25)
else
    
    Qab = -Qa * MQN'/MQM ;        % eq.(60)
    Qb  = inv(MQM) - MQM \ MQN * Qab;  % eq.(61)
end


