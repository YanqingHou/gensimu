% 
% This user interface is meant to visualize geodetic design parameters
% for GNSS.
% The user can specify system, observations, locations, etc.
% Initially, the GPS configuration at the current time is shown, together
% with the satellite orbits for the current day.
% The user can choose whether to compute:
% - a time series for one location:  start and end time together with time
%                                    interval specify the times for which  
%                                    output is computed; number of epochs  
%                                    must be specified by the user 
%                                    (same for all times)!
% - the output for an area or world: start and end time together with the  
%                                    number of epochs specify observation 
%                                    interval,
%                                    the satellite configuration at start  
%                                    time is used! 
% _________________________________________________________________________
%
% TIPS AND WARNINGS
% _________________________________________________________________________
% 
% A manual with a description of underlying models / equations is
% available as PDF-document (VisualManual.pdf) and is titled:
% " Visualization of GNSS-related design parameters "
%
% The results will be saved in the file TMP.MAT
% Computations may take a while, especially for worldwide computations and
% high resolution (grid size < 5 degrees)!
%  
% _________________________________________________________________________
%
% POSSIBLE SETTINGS and <DEFAULTS>
% _________________________________________________________________________
% 
% system               : <GPS> / GLONASS / GALILEO / GPS + GALILEO
%
% almanac file         : <yumaGPS.txt>
%                        file from which satellite ephemerides is read,
%                        should be a recent yuma file, file name must start 
%                        with yuma.
%                        Yuma files can be obtained for GPS from 
%                        http://www.navcen.uscg.gov/ftp/GPS/almanacs/Yuma/
%                        For GLONASS from ftp://polaris.feld.cvut.cz/pub/almanac/
%                        For GALILEO a yuma almanac is available generated based
%                        on the nominal design constellation: yumaGAL.txt
%                        Also a RINEX navigation file may be used.
%                        In either case the file must be in a directory
%                        in the MATLAB path!
% almanac file 2       : <yumaGAL.txt>
%                        if system = GPS + GALILEO the Galileo almanac file
%                        must be specified here.
%
% start date & time    : <current day 0:00>, 
%                        syntax example: 01-Jan-2002 12:00:00
%
% end date & time      : <current day 24:00>, defines the observation interval                        
%
% number of epochs     : <1>
%
% time interval        : <300> seconds, for computations for one point
%                        results are computed for each time stamp
%                        defined by this interval, assuming that
%                        the measurements start at that times 
%                        and with the number of epochs as given
%
% cutoff elevation     : <15> degrees 
%
% scenario             : single point / <single baseline> / geometry free
%
% receiver             : <stationary> / roving
%
% ____________________
%
% ionosphere           : <fixed> / float / weighted                       
%
% ionospheric weight   : <standard deviation> : in meters
%                        baseline length      : in kilometers, the standard
%                                               deviation is computed with: 
%                                               sigma = 0.68 * baseline/1000
%                        only when the weighted ionosphere model is chosen!
%
% troposphere          : <fixed> / float
%                        With troposphere 'fixed' the double difference zenith 
%                        wet delays are assumed known / absent. With
%                        troposphere 'float' these delays are estimated
%                        (every epoch).
% 
% mapping function     : <1/cos(z)> / Ifadis
%
% ____________________
%
% observations         : user can specify which observations (phase/code)  
%                        on which frequencies will be considered.  
%                        In case of dual-frequency GPS it is possible to 
%                        choose cross-correlated code observations.
%                        Also the undifferenced standard deviations (meters)
%                        must be specified. Single difference observations  
%                        are assumed to be uncorrelated! 
%                        If elevation-dependent model is to be used,
%                        non-zero values for 'a' must be chosen. The
%                        standard deviation 's(el)' for a satellite with 
%                        elevation angle 'el' is then given as:
%                     
%                           s(el) = st.dev * ( 1 + a * exp{-el/e0} )
%
%                        Depending on system choice the observations
%                        on the following frequencies [MHz] can be used:
%
%                        GPS    : L1 (1575.42), L2 (1227.60), L5 (1176.45)
%                        GLONASS: L1 (1602.00), L2 (1246.00)
%                        GALILEO: E1 (1575.42), E5a(1176.45), 
%                                 E5b(1207.14), E6 (1278.75)
%
% ____________________
%
% location             : one point: the entire time interval as defined by  
%                                   start and end time will be considered
%                        area     : computation for start time for user 
%                                   defined area
%                        <world>  : computation for start time for entire 
%                                   earth
%
%                        <point in map>: user can point desired point or 
%                                        area in map after all settings are  
%                                        made and the APPLY button is pushed
%                        type          : user can define desired point or  
%                                        area by specifying the minimum and  
%                                        maximum longitude and latitude 
%
% resolution           : <5> degrees, not required for 'one point'
%
% station height       : <0> meters (equal for all locations)
%
% ____________________
%
% output               : MDB             : minimal detectable biases [meter] 
%                        MDE             : minimal detectable errors [meter]
%                        BNR             : bias-to-noise ratios
%                        success rates   : probability of correct ambiguity 
%                                          fixing
%                        bias-affected success rates.
%                        MDE+success rates (no output plotted)
%                        PDOP            : position dilution of precision
%                        GDOP            : geometry dilution of precision
%                        <number of 
%                           satellites>  : number of visible satellites
%                        satellite tracks: tracks from start to end time and  
%                                          position at start time 
%
%                        if MDB, MDE, BNR or biased success rates is
%                        chosen, the user must specify the following
%                        parameters: 
%                        
%                        type of error     : <carrier slip> / code outlier   
%                        on which frequency: <L1>/L2/L5 (<E1>/E5a/E5b/E6)
%                        in epoch          : <1>, arbitrary as long as it 
%                                            is smaller than the number of 
%                                            epochs
%                                            for a cycle slip, it specifies
%                                            the first epoch of occurance
%                        
%                        the output is then given for the 'worst' satellite 
%
%                        DOP values are plotted with a maximum of 10.
%
% save output to file : optionally the output is stored in a .MAT file named
%                       <tmp.mat> (in working directory)
%
% plot output in separate window: optional if 'world' or 'area' is chosen
%                       for location
%
% _____________________________________________________________________________
%
% BUTTONS
% ______________________________________________________________________________
%
% APPLY               : when all settings are OK, the computations can be 
%                       carried out by pushing this button.
% 
% DEFAULTS            : user interface will be restarted (may take a while)
% 
% HELP                : shows help window
%
% CLOSE               : closes the interface
%
% ________________________________________________________________________________
%
% File                 : visual.m
% Creator              : Sandra Verhagen (A.A.Verhagen@TUDelft.nl)
% Date                 : February 2006
% Version              : 1.0
% ________________________________________________________________________________
%
% Hint for people who want to read and/or modify the code:
% the eq.nrs. in the comment refer to equations in the documentation that comes
% with this user-interface: " Visualization of GNSS-related design parameter ".
% ________________________________________________________________________________
%
% Matlab-routines needed for this user-interface:
%                      - the SKYLAB routines
%                      - plh2xyzwgs.m
%                      - cpelev1.m
%                      - plworldini.m
%                      - cplos.m
%                      - decorrel.m (from the LAMBDA-routines)
%                      - ldldecom.m ( "        "         "   )
%
% Available Yuma almanacs: 
%                      - yumaGPS.txt (GPS)
%                      - yumaGLO.txt (GLONASS)
%                      - yumaGAL.txt (GALILEO)
% 
% See also: SKYLAB

