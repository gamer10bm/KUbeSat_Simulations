function [lat, lon] = ECEF2latlon(R_IJK)
% [lat, lon] = ECEF2latlon(R_IJK)
%
% Inputs: R_IJK - Position Vector in Earth Centered Earth Fixed (km)
%
% Outputs: lat - Latitude (deg)
%                lon - Longitude (deg)
%
% Created by Bailey Miller 2/12/2020
%

%Parse R vector
ri = R_IJK(1); rj = R_IJK(2); rk = R_IJK(3);
%Calculate rdsat
rdsat = sqrt(ri^2+rj^2);
%Find alpha
alpha = acos(ri/rdsat); %radians
if rj < 0
    alpha = -alpha;
end
%Set longitude to alpha
lon = rad2deg(alpha);
%Find d
d = atan2(rk,rdsat);
%Iterate to find latitude
epsilon = 1e-12;
phi_gd = d;
phi_gd_old = phi_gd+epsilon*8;
Re = 6378.1363; %km from Vallado
e_earth=0.0033528131; %~ from Vallado
i = 0;
while abs(phi_gd-phi_gd_old)>epsilon
    i = i+1;
    phi_gd_old = phi_gd;
    Cbody = Re/(sqrt(1-e_earth^2*sin(phi_gd)^2));
    phi_gd = atan2((rk+Cbody*e_earth^2*sin(phi_gd)),rdsat);
end
lat = rad2deg(phi_gd);
end
    

