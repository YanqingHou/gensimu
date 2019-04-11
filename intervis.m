function intervis (file,tsat,station1,station2,eph,cutoff);
%INTERVIS: Find satellites in common view between two stations
%
% Syntax:
%    intervis (file,tsat,station1,station2,eph,cutoff);
%
% The function create a table showing the number of satellites visible at
% each of the stations. Also the satellites (number and PRN's) of the
% satellites in common view are shown in this table.
%
% Input arguments:
%    file     - Filename of file to be created
%               empty string ==> write to screen
%    tsat     - Times on which positions are computed
%    station1 - Coordinates (WGS84/XYZ) of first station
%    station2 - Coordinates (WGS84/XYZ) of second station
%    eph      - Satellite ephemeris
%    cutoff   - Selected cutoff elevation
%
% Output arguments:
%    none

% ----------------------------------------------------------------------
% File.....: intervis.m
% Date.....: 01-JUN-1999
% Version..: 1.0
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

% -------------------------------------------------------
% --- Compute satellite positions (XYZ WGS84 & AZ/EL) ---
% -------------------------------------------------------

[xsat,ysat,zsat,azim1,elev1,rotmat] = cpaziele (tsat,eph,station1);
[xsat,ysat,zsat,azim2,elev2,rotmat] = cpaziele (tsat,eph,station2);

% -----------------
% --- Open file ---
% -----------------

if length(file) ~= 0;
  fid = fopen (file,'wt'); 
else;
  fid = 1;
end;

% ---------------------------------
% --- Write station information ---
% ---------------------------------

plh1 = xyz2plh(station1);
plh2 = xyz2plh(station2);

fprintf (fid,['|---------------------|---------------------|----------------------|\n', ...
	      '| Station 1           | Station 2           | Difference           |\n', ...
	      '|---------------------|---------------------|----------------------|\n', ...
	      '| X | %15.3f | X | %15.3f | DX | %15.3f |\n', ...
	      '| Y | %15.3f | Y | %15.3f | DY | %15.3f |\n', ...
	      '| Z | %15.3f | Z | %15.3f | DZ | %15.3f |\n', ...
	      '|---------------------|---------------------|----------------------|\n', ...
	      '| P | %15.6f | P | %15.6f | DP | %15.6f |\n', ...
	      '| L | %15.6f | L | %15.6f | DL | %15.6f |\n', ...
	      '| H | %15.3f | H | %15.3f | DH | %15.3f |\n', ...
	      '|---------------------|---------------------|----------------------|\n\n\n'], ...
	      station1(1),station2(1),station2(1)-station1(1), ...
	      station1(2),station2(2),station2(2)-station1(2), ...
	      station1(3),station2(3),station2(3)-station1(3), ...
	      rad2deg(plh1(1)),rad2deg(plh2(1)),rad2deg(plh1(1)-plh2(1)), ...
	      rad2deg(plh1(2)),rad2deg(plh2(2)),rad2deg(plh1(2)-plh2(2)),  ...
	      plh1(3),plh2(3),plh1(3)-plh2(3));
	 
% ------------------------------------
% --- Write visibility information ---
% ------------------------------------

fprintf (fid,'|--------------------------|---------|---------|----------|-----------------------------------------------|\n');
fprintf (fid,'|          Date/time [UTC] | # stat1 | # stat2 | # common | Common satellites                             |\n');
fprintf (fid,'|--------------------------|---------|---------|----------|-----------------------------------------------|\n');

for i = 1:length(tsat);

    str1 = gpst2str (tsat(i));
    num1 = length(find (elev1(:,i) > cutoff));
    num2 = length(find (elev2(:,i) > cutoff));
    i    = find (elev1(:,i)>cutoff & elev2(:,i)>cutoff);
    str2 = sprintf (' %2d',eph(i,1));
    str2 = [str2 char(32*ones(1,45-3*length(i)))];

    fprintf (fid,'| %24s | %7d | %7d | %8d | %45s |\n',str1,num1,num2,length(i),str2);
    
end;

fprintf (fid,'|--------------------------|---------|---------|----------|-----------------------------------------------|\n');

% ------------------
% --- Close file ---
% ------------------

if fid ~= 1;
  fclose (fid);
end;

% --------------------------------
% --- End of function INTERVIS ---
% --------------------------------