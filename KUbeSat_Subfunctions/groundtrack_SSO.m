function [lat, lon] = groundtrack_SSO(i,t,p_SC,p_body,lat0,lon0)
% [lat, lon] = groundtrack(a,e,nu,raan,w,i)
% 
% This function is used to calculate the groundtrack of an orbiting
% satellite in a Sun Synchronous orbit
%
% Inputs: i - Inclination (deg)
%         t - Time since simulation start (s)
%         p_SC - Period of the spacecraft orbit (s)
%         p_body - Period of rotation of orbiting body (s)
%         lon0 - Starting Longitude (deg)
%
% Outputs: lat - Latitude (deg)
%          lon - Longitude (deg)
%
% Algorithm based on
% http://fgg-web.fgg.uni-lj.si/~/mkuhar/Pouk/SG/Seminar/Vrste_tirnic_um_Zemljinih_sat/Orbit_and_Ground_Track_of_a_Satellite-Capderou2005.pdf 
% Created by Bailey Miller 2/11/2020

%Calculate new latitude
lat = lat0+asind(sind(i)*sin(2*pi/p_SC*t));

%Calculate new longitude
raan_dot = 2*pi*(1/p_body-1/p_SC);
lon = lon0+atand(cosd(i)*tan(2*pi/p_SC*t))-rad2deg(2*pi/p_body*t)+rad2deg(raan_dot*t);
end