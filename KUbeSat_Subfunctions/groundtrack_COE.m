function [lat, lon] = groundtrack_COE(a,e,nu,raan,w,i,lon0)
% [lat, lon] = groundtrack(a,e,nu,raan,w,i)
% 
% This function is used to calculate the groundtrack of an orbiting
% satellite.
%
% Inputs: a - Semimajor axis (km)
%         e - Eccentricity (~)
%         nu - True Anomaly (deg)
%         raan - Right Ascension of the Ascending Node (deg)
%         w - Argument of Perigee (deg)
%         i - Inclination (deg)
%         lon0 - Starting Longitude (deg)
%
% Outputs: lat - Latitude (deg)
%          lon - Longitude (deg)
%
% Algorithm based on
% http://fgg-web.fgg.uni-lj.si/~/mkuhar/Pouk/SG/Seminar/Vrste_tirnic_um_Zemljinih_sat/Orbit_and_Ground_Track_of_a_Satellite-Capderou2005.pdf 
% Created by Bailey Miller 2/11/2020

%Calculate orbital radius (r)
r = a.*(1-e^2)./(1+e.*cosd(nu));

%Calculate transformation matrix terms
term1 = cosd(raan)*cosd(w+nu)-sind(raan)*sind(w+nu)*cosd(i);
term2 = sind(raan)*cosd(w+nu)+cosd(raan)*sind(w+nu)*cosd(i);
term3 = sind(w+nu)*sind(i);

%Determine X Y Z coordinates
coordmat = r.*[term1;term2;term3];
x = coordmat(1); y = coordmat(2); z = coordmat(3);

%Use coordinates to solve for lat and lon
lat = asind(z/r);
if y>0
    lon = acosd(x/(r*cosd(lat)));
else
    lon = lon0+acosd(x/(r*cosd(lat)));
end
end