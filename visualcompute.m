% visualcompute
% perform the actual computations given the GUI settings
% called by VISUAL.M
%
% See also VISUAL

figNumber = watchon;

lc  =  get(findobj(gcf,'Tag','locchoice'),'Value');
if lc == 2 && get(findobj(gcf,'Tag','point'),'Value')

    % select area by pointing in map
    watchoff(figNumber);
    uiwait(msgbox ('Click and drag your mouse to select the desired area','','modal'));
    figure(figNumber); zoom off;
    waitforbuttonpress;
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
    x      = [point1(1,1);point2(1,1)];
    y      = [point1(1,2);point2(1,2)];
    figNumber = watchon;

    xmin = round(min(x));
    xmax = round(max(x));
    ymin = round(min(y));
    ymax = round(max(y));

    reso = str2double(get (findobj(gcf,'Tag','resol'),'String'));
    grint = reso;

    pdeg = ymax:-grint:ymin;
    ldeg = xmin:grint:xmax;
    prad = pdeg*pi/180;
    lrad = ldeg*pi/180;

    set (findobj(gcf,'Tag','long1'),'String',num2str(xmin));
    set (findobj(gcf,'Tag','long2'),'String',num2str(xmax));
    set (findobj(gcf,'Tag','lat1'),'String',num2str(ymin));
    set (findobj(gcf,'Tag','lat2'),'String',num2str(ymax));

elseif lc == 1 && get(findobj(gcf,'Tag','point'),'Value')

    % select one point by pointing in map
    watchoff(figNumber);
    uiwait(msgbox ('Click to select the desired location','','modal'));
    figure(figNumber); zoom off;
    [xmin,ymin]= ginput(1); figNumber = watchon;
    prad = ymin * pi /180;
    lrad = xmin * pi /180;

    set (findobj(gcf,'Tag','long1'),'String',num2str(xmin));
    set (findobj(gcf,'Tag','long2'),'String','');
    set (findobj(gcf,'Tag','lat1'),'String',num2str(ymin));
    set (findobj(gcf,'Tag','lat2'),'String','');

elseif get(findobj(gcf,'Tag','type'),'Value') && lc ==1

    % give desired latitude(s) and longitude(s) by typing

    ymin = str2double(get (findobj(gcf,'Tag','lat1'),'String'));
    ymax = str2double(get (findobj(gcf,'Tag','lat2'),'String'));
    xmin = str2double(get (findobj(gcf,'Tag','long1'),'String'));
    xmax = str2double(get (findobj(gcf,'Tag','long2'),'String'));
    prad = ymin*pi/180;
    lrad = xmin*pi/180;


elseif get(findobj(gcf,'Tag','type'),'Value')&& lc == 2

    % area defined typing

    y=[str2double(get (findobj(gcf,'Tag','lat1'),'String'));...
        str2double(get (findobj(gcf,'Tag','lat2'),'String'))];
    x=[str2double(get (findobj(gcf,'Tag','long1'),'String'));...
        str2double(get (findobj(gcf,'Tag','long2'),'String'))];
    ymin = min(y);
    ymax = max(y);
    xmin = min(x);
    xmax = max(x);

    reso = str2double(get (findobj(gcf,'Tag','resol'),'String'));
    grint = reso;
    pdeg  = ymax:-grint:ymin;
    ldeg  = xmin:grint:xmax;
    prad  = pdeg*pi/180;
    lrad  = ldeg*pi/180;


elseif lc == 3

    % computations for entire world with defined resolution

    reso = str2double(get (findobj(gcf,'Tag','resol'),'String'));
    xmin = -180; xmax = 180;
    ymin = -90; ymax = 90;
    grint = reso * pi/180;
    prad = pi/2 : -grint : -pi/2;
    lrad = -pi : grint : pi;

end

height = str2double(get(findobj(gcf,'Tag','height'),'String'));

% find out which output parameter must be computed;
% if satellite tracks are chosen --> plot
figure(figNumber);
out = get(findobj(gcf,'Tag','output'),'Value');
if out == 10
    visual('initialFigure','begin')
    watchoff(figNumber);
    return
end

% frequencies:
% frequencies depend on system choice
SYSTEM = get(findobj(gcf,'Tag','system'),'Value');
rectype = 'non-cc';


if SYSTEM == 2 % GLONASS

    frq = [1602e6;1246e6];

else     % GPS / Galileo

    frq = [1575.42e6;1227.60e6;1176.45e6;1227.60e6;...
        1575.42e6; 1278.75e6;1207.14e6;1176.45e6];

    % receiver type: cross-correlating when cross-correlated code
    % observation is chosen
    if (get(findobj(gcf,'Tag','pL2L1'),'Value'))
        rectype = 'cc';
    end

end

% measurement scenario:
% - single baseline, stationary or roving receiver (sr of rr) (user defined)
% - geometry free
% - single point

if (get(findobj(gcf,'Tag','receiver'),'Value'))==1
    scenario = 'sr';
else
    scenario = 'rr';
end
if (get(findobj(gcf,'Tag','scenario'),'Value'))==3
    scenario = 'gf';
    model = 'baseline';
    sd    = 2;
elseif (get(findobj(gcf,'Tag','scenario'),'Value'))==1
    model = 'single';
    sd    = 1;
else
    model = 'baseline';
    sd    = 2;
end

% standard deviations of observations (user defined)
freqsel  = [get(findobj(gcf,'Tag','pL1'),'Value');...
    get(findobj(gcf,'Tag','pL2'),'Value');...
    get(findobj(gcf,'Tag','pL5'),'Value');...
    get(findobj(gcf,'Tag','pL2L1'),'Value');...
    get(findobj(gcf,'Tag','pE1'),'Value');...
    get(findobj(gcf,'Tag','pE6'),'Value')
    get(findobj(gcf,'Tag','pE5b'),'Value');...
    get(findobj(gcf,'Tag','pE5a'),'Value')];
sigcode  = [str2double(get(findobj(gcf,'Tag','sigc1'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigc2'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigc3'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigc4'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigc1b'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigc4b'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigc3b'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigc2b'),'String'))];
sigphase = [str2double(get(findobj(gcf,'Tag','sigp1'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigp2'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigp3'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigp4'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigp1b'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigp4b'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigp3b'),'String'));...
    str2double(get(findobj(gcf,'Tag','sigp2b'),'String'))];
acode  = [str2double(get(findobj(gcf,'Tag','ac1'),'String'));...
    str2double(get(findobj(gcf,'Tag','ac2'),'String'));...
    str2double(get(findobj(gcf,'Tag','ac3'),'String'));...
    str2double(get(findobj(gcf,'Tag','ac4'),'String'));...
    str2double(get(findobj(gcf,'Tag','ac1b'),'String'));...
    str2double(get(findobj(gcf,'Tag','ac4b'),'String'));...
    str2double(get(findobj(gcf,'Tag','ac3b'),'String'));...
    str2double(get(findobj(gcf,'Tag','ac2b'),'String'))];
aphase = [str2double(get(findobj(gcf,'Tag','ap1'),'String'));...
    str2double(get(findobj(gcf,'Tag','ap2'),'String'));...
    str2double(get(findobj(gcf,'Tag','ap3'),'String'));...
    str2double(get(findobj(gcf,'Tag','ap4'),'String'));...
    str2double(get(findobj(gcf,'Tag','ap1b'),'String'));...
    str2double(get(findobj(gcf,'Tag','ap4b'),'String'));...
    str2double(get(findobj(gcf,'Tag','ap3b'),'String'));...
    str2double(get(findobj(gcf,'Tag','ap2b'),'String'))];
e0code  = [str2double(get(findobj(gcf,'Tag','ec1'),'String'));...
    str2double(get(findobj(gcf,'Tag','ec2'),'String'));...
    str2double(get(findobj(gcf,'Tag','ec3'),'String'));...
    str2double(get(findobj(gcf,'Tag','ec4'),'String'));...
    str2double(get(findobj(gcf,'Tag','ec1b'),'String'));...
    str2double(get(findobj(gcf,'Tag','ec4b'),'String'));...
    str2double(get(findobj(gcf,'Tag','ec3b'),'String'));...
    str2double(get(findobj(gcf,'Tag','ec2b'),'String'))];
e0phase = [str2double(get(findobj(gcf,'Tag','ep1'),'String'));...
    str2double(get(findobj(gcf,'Tag','ep2'),'String'));...
    str2double(get(findobj(gcf,'Tag','ep3'),'String'));...
    str2double(get(findobj(gcf,'Tag','ep4'),'String'));...
    str2double(get(findobj(gcf,'Tag','ep1b'),'String'));...
    str2double(get(findobj(gcf,'Tag','ep4b'),'String'));...
    str2double(get(findobj(gcf,'Tag','ep3b'),'String'));...
    str2double(get(findobj(gcf,'Tag','ep2b'),'String'))];

ifreq    = find(freqsel==1);

freq  = frq(ifreq);
if SYSTEM == 4
    nfreq = length(find(ifreq<5)); % nfreq = no. of GPS frequencies
else
    nfreq = length(ifreq);
end

if get(findobj(gcf,'Tag','codesel'),'Value')
    ncode  = nfreq;
    ncode2 = length(ifreq)-ncode;
else
    ncode  = 0;
    ncode2 = 0;
end
if get(findobj(gcf,'Tag','phasesel'),'Value')
    nphase  = nfreq;
    nphase2 = length(ifreq)-nphase;
else
    nphase = 0;
    nphase2= 0;
end

% wavelengths
lambda = 299792458 ./ freq;
if strcmp(rectype,'cc')
    lambda(2) = 299792458 / (freq(1)-freq(2));
end
LAMBDA = diag(lambda);

sigcode = sigcode(ifreq);
sigphase= sigphase(ifreq);
acode   = acode(ifreq);
aphase  = aphase(ifreq);
e0code  = e0code(ifreq);
e0phase = e0phase(ifreq);

sig  = [sigcode(1:nfreq);sigphase(1:nfreq);sigcode(nfreq+1:end);sigphase(nfreq+1:end)];
a    = [acode(1:nfreq);aphase(1:nfreq);acode(nfreq+1:end);aphase(nfreq+1:end)];
e0   = [e0code(1:nfreq);e0phase(1:nfreq);e0code(nfreq+1:end);e0phase(nfreq+1:end)];

% get ionospheric standard deviation
if out<7
    ic = get(findobj(gcf,'Tag','ioncase'),'Value');
    if ic == 3 % ionosphere weighted

        if get(findobj(gcf,'Tag','sigion'),'Value')
            sdion = str2double(get(findobj(gcf,'Tag','sdionweight'),'String'));
        elseif get(findobj(gcf,'Tag','baseline'),'Value')
            sdion = 0.68 *(frq(1)/freq(1))* ...
                str2double(get(findobj(gcf,'Tag','blionweight'),'String')) / 1000;
        end

    elseif ic == 1 % ionosphere fixed

        sdion = 0;

    else % ionosphere float

        sdion = 99999;

    end

    % single difference ionospheric vc-matrix

    mu    = ([lambda(1:nfreq)]./ lambda(1)).^2;
    mu    = [mu(1:ncode);-mu(1:nphase)];
    s2    = sd*sdion^2;
    
    CI    = s2 * mu * mu';

    if SYSTEM ==4
        mu2   = (lambda(nfreq+1:end)./ lambda(nfreq+1)).^2  ;
        mu2   = [mu2(1:ncode2);-mu2(1:nphase2)];        
        CI2   = s2 * mu2 * mu2';        
    end


    % troposphere model: float or fixed
    % if tropo = 'float' : mapping function = mapfun

    if (get(findobj(gcf,'Tag','tropocase'),'Value'))==1
        tropo = 'Tfixed';
    else
        tropo = 'Tfloat';
        if get(findobj(gcf,'Tag','mapfun'),'Value')==1
            mapfun = 's';  % 1/cos(z)
        else
            mapfun = 'I';  % Ifadis
        end
    end

end
% get ephemerides and compute satellite positions

almanac       = get(findobj(gcf,'Tag','almanac'),'String');
starttime  = [get(findobj(gcf,'Tag','startdate'),'String') ' ' ...
   get(findobj(gcf,'Tag','starttime'),'String')] ;
endtime    = [get(findobj(gcf,'Tag','enddate'),'String') ' ' ...
   get(findobj(gcf,'Tag','endtime'),'String')] ;
no_epochs     = str2double(get(findobj(gcf,'Tag','nepochs'),'String'));

if length(prad) ==1
    interval   = str2double(get(findobj(gcf,'Tag','intv'),'String'));
    tsat       = mktsat ( starttime,endtime,interval);
else
    tsat       = mktsat ( starttime,starttime,1);
    tsat       = tsat(1);
end
if SYSTEM ==4
    eph        = rdyuma(almanac);
    almanac2   = get(findobj(gcf,'Tag','almanac2'),'String');
    eph2       = rdyuma(almanac2);
    i          = find(eph(:,27)==0); % find healthy satellites
    eph        = eph(i,:);
    i          = find(eph2(:,27)==0); % find healthy satellites GALILEO
    eph2       = eph2(i,:);
    [xs,ys,zs] = cpaziele (tsat,eph);
    [xs2,ys2,zs2] = cpaziele (tsat,eph2);
else
    if strcmp(almanac(1:4),'yuma')
        eph = rdyuma(almanac);
    else
        ephtmp =rxnav(almanac);
        ephtmp = sortrows(ephtmp);
        eph(1,:)=ephtmp(1,:);
        j=1;
        for i=2:length(ephtmp)
            if ephtmp(i,1)==ephtmp(i-1,1),
                continue
            else
                j=j+1;
                eph(j,:)=ephtmp(i,:);
            end
        end
    end
    i          = find(eph(:,27)==0); % find healthy satellites
    eph        = eph(i,:);
    [xs,ys,zs] = cpaziele (tsat,eph);
end

% cutoff elevation and number of epochs
cutoff     = str2double(get(findobj(gcf,'Tag','cutoff'),'String'));

% testing parameter (default)
lam0 = 17.0747;


% MDB/BNR/MDE
if ismember(out,[1 2 3 4 5 6])
    % observation for which to compute output:
    % carrier slip or code outlier
    if (get(findobj(gcf,'Tag','errtype'),'Value'))==1
        output = 'slip';
    else
        output = 'outlier';
    end
    outfreq = (get(findobj(gcf,'Tag','outputfreq'),'Value'));
    if (SYSTEM==4) && (outfreq > 3)
        outputfreq = find(freq==frq(outfreq+1));
        if length(outputfreq)>1, outputfreq=outputfreq(2); end
    elseif SYSTEM==3
        outputfreq = find(freq==frq(outfreq+4));
    else
        outputfreq = find(freq==frq(outfreq));
        outputfreq = outputfreq(1);
    end

    if (outputfreq > nfreq) && (SYSTEM ~= 4)
        h2 = msgbox ('Error: output cannot be computed with the number of frequencies as given');
        watchoff(figNumber);
    elseif isempty(outputfreq)
        h2 = msgbox ('Error: wrong frequency choice for outlier / cycle slip');
        watchoff(figNumber);
    elseif (outputfreq > nfreq) && (SYSTEM == 4)
        sys        = 2;
        ml         = lambda(outputfreq);
        outputfreq = outputfreq - nfreq;
        nc         = ncode2;
        np         = nphase2;
    else
        sys        = 1;
        ml         = lambda(outputfreq);
        nc         = ncode;
        np         = nphase;
    end

    startepoch = str2double(get(findobj(gcf,'Tag','startepoch'),'String'));
    if startepoch > no_epochs
        h2 = msgbox ('Error: epoch specified in output frame > number of epochs');
        watchoff(figNumber);
        return
    end
    % definition of canonical vectors that select observation for
    % which to compute output c_{i} or c_{i+f} in eq.(3.11),(3.12)
    if strcmp(output,'slip')
        if np == 0
            h2 = msgbox ('Error: carrier slip not possible with code-only');
            watchoff(figNumber);
            return
        end
        Ifreq    = eye(np);
        c1       = Ifreq(:,outputfreq);
        c2       = zeros(nc,1);
        c        = [c2 ; c1];
        slip     = no_epochs - startepoch + 1;
        cl       = [zeros(no_epochs-slip,1);ones(slip,1)];  % vector sl in eq.(3.12)
    elseif strcmp(output,'outlier')
        if nc == 0
            h2 = msgbox ('Error: code outlier not possible with phase-only');
            watchoff(figNumber);
            return
        end
        Ifreq    = eye(nc);
        c1       = Ifreq(:,outputfreq);
        c2       = zeros(np,1);
        if strcmp(rectype,'cc')&&(outputfreq==1), c1 = [1;1]; end
        c        =  [c1 ; c2];
        slip     = 1;
        cl       = eye(no_epochs);
        cl       = cl(:,slip);  
    end

end

% definition of vector with epochs and size of slip window
if length(prad)==1
    time = 1:size(tsat,2);
else
    time = 1;
end

% initialize output matrix
OUT = zeros(length(prad),length(lrad));
OUT2= OUT;
% re-order vectors with latitudes and longitudes and compute corresponding XYZ-coordinates
lp = length(prad);
ll = length(lrad);
i=1;
j = ll;
for p = 1: lp
    plh(i:j,:) = [prad(p)*ones(ll,1) lrad' height*ones(ll,1)];  % vector [latitudes longitudes heights]
    i = j+1;
    j = i + ll -1;
end
xyz      = plh2xyzwgs(plh);

% set parameters for and create waitbar
wend     = size(tsat,2) + lp * ll - 1;
widx     = 0;
watchoff(figNumber);
wh       = waitbar(0,'Busy...');
widx2    = 1;
% loop over latitudes and longitudes, when only one location:
% loop over all epochs
for t = time             % loop over all start times for option 'one point'

    for p = 1:lp*ll       % loop over all grid points for option 'area' or 'world'

        widx = widx + 1;

        el       = cpelev1(xs(:,t),ys(:,t),zs(:,t),plh(p,:),xyz(p,:));
        idx1     = find (el > cutoff);
        no_sv    = length(idx1);
        el       = el(idx1);
        refsat   = find(el==max(el));
        xsv      = xs(idx1,t);
        ysv      = ys(idx1,t);
        zsv      = zs(idx1,t);
        zen      = deg2rad(90 - el)';
        no_sv2= 0;
        if SYSTEM == 4
            el2       = cpelev1(xs2(:,t),ys2(:,t),zs2(:,t),plh(p,:),xyz(p,:));
            idx2      = find (el2 > cutoff);
            el2       = el2(idx2);
            no_sv2    = length(idx2);
            refsat2   = find(el2==max(el2));
            xsv2      = xs2(idx2,t);
            ysv2      = ys2(idx2,t);
            zsv2      = zs2(idx2,t);
            zen2      = deg2rad(90 - el2)';
        end
        if out   == 9

            % --------------------
            % number of satellites
            % --------------------

            OUT(p) = no_sv + no_sv2;
            titel    = ' number of satellites ';

            % ____________________

        else
            OUT2(p) = no_sv + no_sv2;
            LOS = zeros(no_sv,3);
            if ~strcmp(scenario,'gf')
                LOS = cplos([xsv ysv zsv],xyz(p,:));   % LOS-vectors
            end
            LOS2 = [];
            if SYSTEM == 4
                LOS2 = cplos([xsv2 ysv2 zsv2],xyz(p,:));
                m2  = no_sv2 - 1;
                if strcmp(model,'single')
                    m2 = no_sv2;
                end
            end
            m  = no_sv - 1;
            if strcmp(model,'single')
                m = no_sv;
            end

            if out < 7  % for internal or external reliability or success rates

                if no_sv<4, OUT(p)=NaN; continue; end

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
                    Mz   = tropomap(zen,plh(3),mapfun);    % tropo. wet delay mapping function
                end
                if ~strcmp(scenario,'gf')

                    if ~strcmp(model,'single')
                        M1 = D * [LOS Mz];           % [\bar{G} \Psi] as in eq.(2.22)
                    else
                        M1 = [LOS em Mz];            % [G em \Psi] as in eq.(2.5)
                    end
                else
                    M1 = Im;
                end
                M = []; % see eq.(19)
                if strcmp(rectype,'non-cc')
                    n1 = nphase + ncode;
                    n2 = 1;
                    Qy = [];
                    while (n1 > 0 )
                        M = [M ; M1];
                        n1 = n1 - 1;
                        %elevation dependent weighting as in eq.(2.13) and (2.15)
                        Wi = sd * (sig(n2) * diag(1 + a(n2) * exp(-el/e0(n2)))).^2;
                        Qy = blkdiag(Qy,D*Wi*D'); %first term after =-sign in eq.(2.20)/(2.21)
                        n2 = n2 + 1;
                    end
                else
                    M = [M1;zeros(size(M1))];
                    W1 = sd * (sig(1) * diag(1 + a(1) * exp(-el/e0(1)))).^2;
                    W2 = sd * (sig(2) * diag(1 + a(2) * exp(-el/e0(2)))).^2;
                    Qy = blkdiag(D*W1*D',D*W2*D');
                    if nphase > 0
                        M = [M;M1;M1];
                        W1 = sd * (sig(3) * diag(1 + a(3) * exp(-el/e0(3)))).^2;
                        W2 = sd * (sig(4) * diag(1 + a(4) * exp(-el/e0(4)))).^2;
                        Qy = blkdiag(Qy,D*W1*D',D*W2*D');
                    end
                end
                if sdion > 0
                    QI   = kron(CI,D*D');       % 2nd term after =-sign in eq.(2.20)/(2.21)
                else
                    QI   = zeros(size(Qy));
                end
                iQy = inv(Qy+QI);               % inverse of matrix in eq.(2.20)/(2.21)
                N =[];  % see eq.(19)
                if nphase>0
                    N1 = zeros(ncode*m);
                    N2 = kron(LAMBDA(1:nphase,1:nphase),Im);
                    N = [N1;N2];
                end
                NQN  = N' * iQy * N;
                MQN  = M' * iQy * N;
                MQM  = M' * iQy * M;
                iMQM = inv(MQM);
                
                % also for GALILEO if system=GPS + GALILEO
                if SYSTEM == 4
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
                        Mz2   = tropomap(zen2,plh(3),mapfun);    % tropo. wet delay mapping function
                    end
                    if ~strcmp(scenario,'gf')
                        if ~strcmp(model,'single')
                            M21 = D2 * [LOS2 Mz2];       % [\bar{G} \Psi] as in eq.(2.22)
                        else
                            M21 = [LOS2 em2 Mz2];        % [G em \Psi] as in eq.(2.5)
                        end
                    else
                        M21 = Im2;
                    end
                    M2 = [];   % see eq.(2.22)
                    n1 = nphase2 + ncode2;
                    Qy2 = [];
                    while (n1 > 0 )
                        M2 = [M2 ; M21];
                        n1 = n1 - 1;
                        %elevation dependent weighting as in eq.(2.13)and (2.15)
                        Wi = sd * (sig(n2) * diag(1 + a(n2) * exp(-el2/e0(n2)))).^2;

                        Qy2 = blkdiag(Qy2,D2*Wi*D2');
                        n2 = n2 + 1;
                    end
                    if sdion > 0
                        QI2   = kron(CI2,D2*D2'); % 2nd term after =-sign in eq.(2.20)/(2.21)
                    else
                        QI2   = zeros(size(Qy2));
                    end
                    iQy2 = inv(Qy2+QI2);
                    N2 =[];   % see eq.(2.22)
                    if nphase2>0
                        N21 = zeros(ncode2*m2);
                        N22 = kron(LAMBDA(nphase+1:end,nphase+1:end),Im2);
                        N2 = [N21;N22];
                    end

                    NQN2 = N2' * iQy2 * N2;
                    MQN2 = M2' * iQy2 * N2;
                    MQM2 = M2' * iQy2 * M2;
                    iMQM2= inv(MQM2);

                end

                if SYSTEM == 4 && ~strcmp(scenario,'gf')
                    iMQM  = inv(MQM + MQM2);
                    ins1  = NQN - MQN' * iMQM * MQN;
                    ins12 = -MQN' * iMQM * MQN2;
                    ins2  = NQN2 - MQN2' * iMQM * MQN2;
                    Qa   = inv([ins1 ins12;ins12' ins2]);
                    Q1   = Qa;
                elseif SYSTEM == 4 && strcmp(scenario,'gf')
                    NQNc  = blkdiag(NQN,NQN2);
                    MQNc  = blkdiag(MQN,MQN2);
                    iMQMc = blkdiag(iMQM,iMQM2);
                    Qa    = inv(NQNc - MQNc' * iMQMc * MQNc);
                    Q1    = inv(NQN - MQN' * iMQM * MQN);
                else
                    Qa    = inv(NQN - MQN' * iMQM * MQN);
                    Q1    = Qa;
                end
                Qahat = (1/no_epochs)*Qa ;% eq.(3.23)
                % ---------------------
                % --- success-rates ---
                % ---------------------
                if ismember(out,[4 5 6])
                    if strcmp(model,'single')
                        h2 = msgbox ('Error: you should choose the baseline scenario');
                        watchoff(figNumber);
                        return
                    end
                    [Qz,Z,LQz,DQz]= decorrel(Qahat);
                    iL = inv(LQz');
                    OUT(p) = prod ( 2 * normcdf(1./(2*sqrt(DQz))) -1 );   % eq.(3.31)

                    if out == 6
                        SR(p) = OUT(p);
                    end

                    titel   = ' success rates ';
                end

                % _________________________

            end   % if-statement output is MDE or success rates

            if ismember(out,[1 2 3 5 6]) % MDBs, MDEs or BNRs

                if sys == 1
                    d  = kron(c,D);
                    dQ = d'*iQy;
                    PN = N*inv(NQN)*N'*iQy;
                    
                    if SYSTEM == 4 && ~strcmp(scenario,'gf')
                        lm  = size(M,1);
                        A   = [[M;M2] blkdiag(N,N2) ];
                        A   = A(1:lm,:);
                        MQN = [MQN MQN2];
                    else
                        A   = [ M N ];
                    end
                    if strcmp(scenario,'sr')
                        PM  = 0;
                    else
                        PM  = M*iMQM*M'*iQy;
                    end
                    im = size(iQy);
                    if strcmp(output,'outlier')
                        Q2  = -iMQM * MQN * Q1;
                        Q3  = iMQM - iMQM * MQN * Q2';
                        PMN = A * [Q3 Q2 ; Q2' Q1] * A' * iQy;
                    else
                        PMN = eye(im);
                    end
                    
                else
                    d  = kron(c,D2);
                    dQ = d'*iQy2;
                    PN = N2*inv(NQN2)*N2'*iQy2;
                    if ~strcmp(scenario,'gf')
                        lm  = size(M,1);
                        A   = [[M;M2] blkdiag(N,N2) ];
                        A   = A(lm+1:end,:);
                        MQN = [MQN MQN2];
                    else
                        A   = [ M2 N2 ];
                    end
                    if strcmp(scenario,'sr')
                        PM  = 0;
                    else
                        PM  = M2*iMQM*M2'*iQy2;
                    end
                    im = size(iQy2);
                    if strcmp(output,'outlier')
                        Q2  = iMQM * MQN * Q1;
                        Q3  = iMQM - iMQM * MQN * Q2';
                        PMN = A * [Q3 Q2 ; Q2' Q1] * A' * iQy2;
                    else
                        PMN = eye(im);
                    end
                    
                    
                end
                % div : denominator in eq. (3.13) and (3.14)
                div = sum( slip * (dQ*(eye(im) - (1-slip/no_epochs) * PM - ...
                    (slip/no_epochs) * PMN))'.*d );

                % -----------------------
                % minimal detectable bias
                % -----------------------

                MDB       = sqrt ( lam0 ./ div );             % eq.(3.13) or (3.14)
                Qhat      = [];
                
                if out == 1
                    if strcmp(output,'slip')
                        MDB = ml * MDB;        % multiplication with wavelength:
                                               % cycle slip MDB in units of range
                    end

                    OUT(p) = max(MDB);
                    titel    = ' minimal detectable biases ';
                    % _________________________

                    % -------------------
                    % bias-to-noise ratio
                    % -------------------

                elseif out == 3
                    % div2: term in parentheses in eq. (3.15) and (3.16)
                    div2 = slip * (dQ*(eye(size(PN)) - (slip/no_epochs) * PN))'.*d;
                    BNR    = (max(MDB)^2) * div2 - lam0 ;      % (3.15) or (3.16)
                    if abs(BNR)<1e-20, BNR = 0; end
                    OUT(p) = sqrt(BNR(1));
                    titel  = ' bias-to-noise ratios (baseline)';

                    % _________________________

                    % -------------------------
                    % minimal detectable effect
                    % -------------------------

                elseif ismember(out,[2 5 6]) % MDE / biased success rates

                    if SYSTEM == 4
                        if ~strcmp(scenario,'gf')
                            M = [M ; M2];
                            N = blkdiag(N,N2);
                            iQy = blkdiag(iQy,iQy2);
                        else 
                            M = blkdiag(M,M2);
                            N = blkdiag(N,N2);
                            iQy = blkdiag(iQy,iQy2);
                        end
                        iMQM = inv(M'* iQy * M);
                        MQN  = M'* iQy * N;
                        if sys == 1
                            d = [d;zeros(length(c)*m2,size(d,2))];
                        else
                            d = [zeros(length(c)*m,size(d,2));d];
                        end
                    end
                    if strcmp(scenario,'gf')||strcmp(scenario,'rr')
                        Ik  = eye(no_epochs);
                        ek  = ones(no_epochs,1);

                        Qab = -Qahat * kron(ek', MQN'*iMQM );      % eq.(3.24)
                        iMQMk = kron(Ik,iMQM);
                        Qb  = iMQMk - iMQMk * kron(ek, MQN) * Qab; %eq.(3.25)
                    else
                        Qab = -Qahat * MQN'*iMQM ;        % eq.(3.24)
                        Qb  = (iMQM - iMQM * MQN * Qab);  % eq.(3.25)
                    end
                    Qhat   = [Qb Qab';Qab Qahat];         % eq.(3.21)
                    nr  = size(N,2);       % number of unknown ambiguities
                    nt  = size(Qhat,1);    % number of unknowns

                    nabx=[];
                    lmdb = length(MDB);
                    mdetmp = zeros(lmdb,1);

                    for jdx = 1:lmdb
                        cnab   = d(:,jdx) * MDB(jdx);

                        % eq.(3.26)
                        if strcmp(scenario,'gf')||strcmp(scenario,'rr')
                            nabx(:,jdx) = Qhat * [ kron(cl , M'*iQy*cnab); slip * N'*iQy*cnab ];
                        else
                            nabx(:,jdx)  = slip*Qhat * [M N]' * iQy * cnab;
                        end
                        mdetmp(jdx) = sqrt (sum(nabx(1:3,jdx).^2));
                    end
                    jdx2 = find(mdetmp==max(mdetmp));
                    if length(jdx2)>1, jdx2=jdx2(1); end

                    if out == 2 % MDE on baseline
                        OUT(p)= mdetmp(jdx);
                        titel    = ' minimal detectable effects (on baseline) ';
                    else
                        bias     = iL*Z'*nabx(nt-nr+1:nt,:);
                        bsr =  [];
                        for si = 1:size(bias,2)
                            bsr(si) = prod ( normcdf((1-2*bias(:,si))./(2*sqrt(DQz'))) ...
                                + normcdf((1+2*bias(:,si))./(2*sqrt(DQz')))-1 );% eq.(3.33)
                        end
                        BSR(p) = min(bsr); % biased success rate
                        if out == 5
                            OUT(p) = BSR(p);
                            titel = ' bias affected success rates ';
                        elseif out == 6
                            MDE(p) = mdetmp(jdx);
                        end
                    end

                end

                % _________________________

            end     % if-statement output is 'MDB', 'MDE' or 'BNR'

            % ---------------------
            % DOPs
            % ---------------------

            if ismember(out,[7 8])
                if no_sv+no_sv2 >= 4;
                    Amat    = [LOS ones(no_sv,1)];
                    if SYSTEM == 4  %eq.(3.36)
                        Amat = [Amat zeros(no_sv,1);...
                            LOS2 zeros(no_sv2,1) ones(no_sv2,1)];
                    end
                    qx      = inv(Amat' * Amat);
                    OUT(p) = sqrt(trace(qx(1:out-4,1:out-4)));  % PDOP: out=7, GDOP: out=8
                    if OUT(p)>10, OUT(p)=10; end                % eq.(3.34) or (3.35)
                    titel    = ' dilution of precision ';
                else
                    OUT(p) = 10; % maximum value for plotting purposes
                end;
            end
        end  % if-statement output is 'number of satellites' or misc.
        % _________________________
        % waitbar
        if widx > widx2*wend/100
            waitbar(widx/wend,wh);
            widx2 = widx2 + 1;
        end
    end   % loop over grid cells

    if lp==1
        if out == 6
            MDEt(t) = MDE(p);
            SRt(t)  = SR(p);
            BSRt(t) = BSR(p);
        else
            OUTt(t) = OUT(p); 
        end

        if out~=5, NO_SV(t)= OUT2(p); else NO_SV(t) = OUT(p); end
    end
end       % loop over time stamps

delete(wh);
if out == 6
    if lp>1
        MDE = reshape (MDE,ll,lp)';
        SR  = reshape (SR,ll,lp)';
        BSR = reshape (BSR,ll,lp)';
    else
        MDE = MDEt;
        SR  = SRt;
        BSR = BSRt;
    end
else
    OUT = reshape (OUT,ll,lp)';
end

% plot OUTPUT
if length(time)==1
    if out~= 9, NO_SV = reshape (OUT2,ll,lp)'; else NO_SV = OUT; end
    if get(findobj(gcf,'Tag','saveout'),'Value')
        filename = get(findobj(gcf,'Tag','filename'),'String');
        if out == 6
            save(filename,'MDE','SR','BSR','NO_SV','xmin','xmax','ymin','ymax', ...
                'starttime')
            return
        else
            save(filename,'OUT','NO_SV','xmin','xmax','ymin','ymax','starttime')
        end
    end


    inter = (ymax-ymin)/6; if inter>1, inter=round(inter); end
    figure(figNumber)
    axes(findobj(gcf,'Tag','ColorBarAxis'));
    cla reset;
    set(gca,'Tag','ColorBarAxis');
    axes(findobj(gcf,'Tag','PictureAxis'));
    cla reset;
    set(gca,'Tag','PictureAxis');

    imagesc([xmin xmax],[ymax ymin],OUT);
    hold on
    grid on
    load coast;
    plot (long,lat,'k-');

    title(titel,'color',[1 1 1]);
    set(gca,'YDir','normal','Color',[1 1 1], 'XColor',[1 1 1], 'YColor',[1 1 1], ...
        'Xgrid','on','Ygrid','on', 'Xtick',xmin:inter:xmax,'Ytick',ymin:inter:ymax, ...
        'Xlim',[xmin xmax],'Ylim',[ymin ymax],...
        'Box','on','Visible','on','Tag','PictureAxis');

    h=colorbar (findobj(gcf,'Tag','ColorBarAxis'));
    set(h,'Tag','ColorBarAxis','Color',[1 1 1], 'XColor',[1 1 1], 'YColor',[1 1 1]);
    zoom on

    % also plot in separate window
    if get(findobj(gcf,'Tag','plotnew'),'Value')
        figure
        imagesc([xmin xmax],[ymax ymin],OUT);
        hold on
        grid on
        load coast;
        plot (long,lat,'k-');

        title(titel);
        set(gca,'YDir','normal','Color',[1 1 1], 'XColor',[0 0 0], 'YColor',[0 0 0], ...
            'Xgrid','on','Ygrid','on', 'Xtick',xmin:inter:xmax,'Ytick',ymin:inter:ymax, ...
            'Xlim',[xmin xmax],'Ylim',[ymin ymax],...
            'Box','on','Visible','on');
        colorbar
    end
else
    if out == 9, NO_SV = OUTt; end
    if get(findobj(gcf,'Tag','saveout'),'Value')
        filename = get(findobj(gcf,'Tag','filename'),'String');
        if out == 6
            save(filename,'MDE','SR','BSR','NO_SV','xmin','ymin', ...
                'starttime','endtime')
            return
        else
            save(filename,'OUTt','NO_SV','xmin','ymin','starttime','endtime')
        end
    end
    
    xmin = 1; xmax = max(time);
    ymin = min(OUTt); ymax= max(OUTt);
    figure
    x=xmin:xmax;
    st=stairs(x,NO_SV,'g');
    set(st,'Linewidth',2)
    set(gca,'Ytick',min(NO_SV):1:max(NO_SV),'Ylim',[min(NO_SV)-1 max(NO_SV)+1]);
    ds = max(NO_SV) - min(NO_SV) + 2;
    ylabel('number of satellites (green)');
    if out ~= 9
        set(gca,'YaxisLocation','right','xlim',[xmin xmax],'xtick',[],'xticklabel',[]);
        hold on
        h=axes('position',get(gca,'position'));

        hold on
        set(h,'color','none')
        set(h,'xgrid','off','ygrid','off');
        set(h,'YaxisLocation','left')
        plot(x,OUTt,'r','Linewidth',2);
        yl=get(gca,'Ylim');
        dy=yl(2)-yl(1);
        yy = 0;
        if ismember(out,[4 5])
            yl(2)=1;
            dy   = yl(2)-yl(1);
            yy   = 0.05*dy/ds;
        end
        if dy<1e-2, yl(1)=yl(1)-1; yy = 0.01; dy = dy+1; end
        set(gca,'xlim',[xmin xmax],'Ytick',yl(1):dy/ds:yl(2),'Ylim',[yl(1) yl(2)+yy]);
        ylabel('output (red)');


    end
    grid on

    xlabel('starting epoch ');

    title(titel);
end
watchoff(figNumber)
