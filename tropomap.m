function [Mz] = tropomap (zen,h,MF);
% [Mz] = tropomap (zen,h,MF); 
% compute the mapping function for the tropospheric wet delay
%
% IN:
%      zen   zenith angle(s) from receiver to satellite [rad]
%      h     station height above mean sea level
%      MF    mapping function to be used:
%            'N' : Niell        [NOT YET IMPLEMENTED]
%            'I' : Ifadis
%            's' : 1 / cos(zen) [DEFAULT}
%
%
%
% Sandra Buur-Verhagen, January 2006

% START:
% -------------------------------------------------------------------------

if strcmp(MF,'I')

% SAASTAMOINEN model with Ifadis mapping function
T0 = 288.15 - 6.5*h/1000; % Temperature of standard atmosphere [Kelvin]
P0 = 1013.25*(288.15/T0)^(-5.255877); % Pressure of standard atmosphere [mbar]
e = 0.5*exp(24.3702 - 6162.3496/T0); % Partial pressure of water vapour [mbar]
                                     % with the assumption that the relative humidity is 50%

% Parameters of the wet mapping function
a = 0.5236*10^(-3) + 0.2471*10^(-6)*(P0 - 1000) - 0.1328*10^(-4)*sqrt(e) + 0.1724*10^(-6)*(T0 - 288.15);
b = 0.1705*10^(-2) + 0.7384*10^(-6)*(P0 - 1000) + 0.2147*10^(-4)*sqrt(e) + 0.3767*10^(-6)*(T0 - 288.15);
c = 0.05917;
% Mapping function by Ifadis
Mz = (1 + a/(1 + b/(1 +c)))./(cos(zen) + a./(cos(zen) + b./(cos(zen) + c)));

elseif strcmp(MF,'s')
    
    Mz = 1./cos(zen);

end