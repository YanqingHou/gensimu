function xefc = qic2efc(tod,xinc)
%QIC2EFC Pseudo "Inertial" coordinates to Earth Fixed coordinates.
%        Converts inertial alike coordinates to Earth fixed 
%        coordinates by applying the rotation of the system
%        due to the rotation of the Earth, tod is the time of
%        day in seconds. 

%        H. van der Marel, LGR, 07-05-95
%        (c) Geodetic Computing Centre, TU Delft

omgedo = 7292115.1467e-11;

c = cos(omgedo*tod);
s = sin(omgedo*tod);

xefc = [ [c,s;-s,c]*xinc(1:2) ; xinc(3) ];

