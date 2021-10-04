
function [InView, range, el] = GS_View(GS, R, time, Re, f)
% Rsez = C*Reci
% GS Input [Lat, Long, Altitude]

% Day fraction
D = time/3600/24;

% Long at prime meridian
thetag = 1.0027379093*2*pi*D;

L = GS(1);
% Long of GS
theta = GS(2) + thetag;

% X & Y of GS in ECEF
x = (Re/(sqrt(1-f^2*sin(L)^2)) + GS(3))*cos(L);
z = (Re*(1-f^2)/(sqrt(1-f^2*sin(L)^2)) + GS(3))*sin(L);

% Station radius vector
Rstat_eci = [x*cos(theta); x*sin(theta); z];

% Transformation Matrix from ECEF to ECI
C = [sin(L)*cos(theta) sin(L)*sin(theta) -cos(L);
    -sin(theta)        cos(theta)         0;
     cos(L)*cos(theta) cos(L)*sin(theta)  sin(L)];

rho = R - Rstat_eci; % Range
rhoSEZ = C*rho;
el = asin(rhoSEZ(3)/norm(rhoSEZ)); % Elevation

range = norm(rho);

% If above horizon
if el >= deg2rad(5) 
    InView = 1;
else
    InView = 0;
end

end