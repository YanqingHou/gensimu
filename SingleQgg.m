function [Qa1,Ps1,Qb1,Qab1,DOPs] = SingleQgg(freq,...
    sigcode,sigphase,sdion,tropo,no_epochs,cutoff,xs,ys,zs,xyz,plh,cfix)


% freq  = [freq1;freq2];

model         = 'baseline';         % baseline
scenario      = 'kin';               % rr-stationary
rectype       = 'non-cc';           % non-crosscorrelating

ncode         = length(freq);      % number of code observations 
nphase        = ncode;              % number of phase observations

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
% s2    = sd*sdion^2;

CI    = mu * mu';

el       = cpelev1(xs,ys,zs,plh,xyz);
idx1     = find (el > cutoff);
no_sv    = length(idx1);
refsat   = find(el(idx1)==max(el(idx1)));
xsv      = xs(idx1);
ysv      = ys(idx1);
zsv      = zs(idx1);
el       = el(idx1);
zen      = deg2rad(90 - el)';




if ~cfix
    LOS  = cplos([xsv ysv zsv],xyz);   % LOS-vectors
else 
    LOS  = [];
end
A=[LOS,ones(no_sv,1)];
Q=inv(A'*A);
TDOP=sqrt(Q(4,4));
PDOP=sqrt(Q(1,1)+Q(2,2)+Q(3,3));
GDOP=sqrt(trace(Q));


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

if isempty(M)
    
    ins1 = NQN;
    
    MQM = []; MQN = []; 
else
 
    MQN = M' * iQy * N;
    MQM = M' * iQy * M;
    iMQM = inv(MQM);
    
    ins1   = NQN - MQN' * iMQM * MQN;

    

end

Qa1 = (1/no_epochs)*inv(ins1); % eq.(58)

Qa1   =(tril(Qa1,0)+tril(Qa1,-1)');


Ps1 = cpsucrate(Qa1,'LS');
ADOP=det(Qa1);
ADOP=ADOP^(1/m);
DOPs=[TDOP,PDOP,GDOP,ADOP];

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
else    
    Qab1 = -Qa1 * MQN'/MQM ;        % eq.(60)
    Qb1  = inv(MQM) - MQM\MQN * Qab1;  % eq.(61)  
end