
function [InView, range, el] = GS_View(GS, R, time)
% Rsez = C*Reci
% GS Input [Lat, Long, Altitude]

D = time/3600/24;

thetag = 1.0027379093*2*pi*D;

L = GS(1);
theta = GS(2) + thetag;

Re = 6378;
mu = 3.986e5;
f = 0.08182;

x = (Re/(sqrt(1-f^2*sin(L)^2)) + GS(3))*cos(L);
z = (Re*(1-f^2)/(sqrt(1-f^2*sin(L)^2)) + GS(3))*sin(L);

Rstat_eci = [x*cos(theta); x*sin(theta); z];

C = [sin(L)*cos(theta) sin(L)*sin(theta) -cos(L);
    -sin(theta)        cos(theta)         0;
     cos(L)*cos(theta) cos(L)*sin(theta)  sin(L)];

rho = R - Rstat_eci;
rhoSEZ = C*rho;
el = asin(rhoSEZ(3)/norm(rhoSEZ));

range = norm(rho);

if el >= deg2rad(5)
    InView = 1;
else
    InView = 0;
end

end

