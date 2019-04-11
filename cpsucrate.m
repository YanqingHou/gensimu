function sucrate = cpsucrate (Qa,estimator)

% sucrate = cpsucrate(Qa,estimator)
% computes success rate for the given variance-covariance matrix and integer estimator
%
% Qa       : variance-covariance matrix of the float ambiguities
% estimator: integer estimator to be used:
%            'R' : integer rounding (computed success rate is a lower bound)
%            'B' : integer bootstrapping (no decorrelation)
%            'LS': integer least-squares (bootstrapped success rate is used as a lower bound)
%                  [default]
%
% Note that 'B' and 'LS' will not give the same results, since with 'LS'
% first a decorrelation is carried out.
%
% Sandra Verhagen, August 2002

if nargin<2, estimator = 'LS'; end

if strcmp(estimator,'R')
    D = diag(Qa);
    sucrate = prod ( 2 * normcdf(1./(2*sqrt(D))) -1 );    
elseif strcmp(estimator,'B')
    [L,D] = ldldecom(Qa);
    sucrate = prod ( 2 * normcdf(1./(2*sqrt(D))) -1 );    
elseif strcmp(estimator,'LS')
    [Qzhat,Z,L,D] = decorrel (Qa);
    sucrate = prod ( 2 * normcdf(1./(2*sqrt(D))) -1 );
else
    error('unknown estimator. See help')
end