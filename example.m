% sample use of the wgs84egsa87 and egsa87wgs84 functions

% https://en.wikipedia.org/wiki/Morosini_Fountain
% latitude (phi) and longitude (lambda) in WGS84 (degrees) 
phi=35.3391370;
lambda=25.1331727;

% to radians
lambda2=lambda/180.0*pi;
phi2=phi/180.0*pi;

% convert to EGSA87
[x y] = wgs84egsa87(phi2, lambda2);
fprintf('Point in EGSA87: %f %f\n', x, y);

% back to WGS84
[phi2 lambda2] = egsa87wgs84(x, y);

% to degrees
phi2=phi2/pi*180.0;
lambda2=lambda2/pi*180.0;

fprintf('Point in WGS84: %f %f (error %g)\n', phi2, lambda2, norm([phi, lambda]-[phi2, lambda2]));
