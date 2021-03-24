function [COEStruct,hmag,i,Omega,emag,lilomega,theta,a,T] = RV2coe(r,v,mu)
% function [COEStruct,hmag,i,Omega,emag,lilomega,theta,a,T] = RV2coe(r,v,mu)
%
% Inputs: r - position vector in km
%             v - veloctiy vector in km/s
%             mu - gravitational parameter (km^3/s^2)
%
% Outputs: COEStruct - Structure filled with COE outputs
%                hmag - specific angular momentum (km^2/s)
%                i - inclination (deg)
%                Omega - Right Ascension of the Ascending Node (deg)
%                emag - Eccentricity (~)
%                lilomega - Argument of Perigee (deg)
%                theta - True Anomaly (deg)
%                a - Semimajor Axis (km)
%                T - Orbital Period (sec)
%
% Created by Bailey Miller 2/12/2020
% see also: coe2RV.m 

if nargin<3
    mu = 3.986004e5; %km^3/s^2
    if nargin < 2
        v = r(4:6);
        r = r(1:3);
    end
end

%% Calculate distance and velocity magnitude 
rmag = norm(r);
vmag = norm(v);

%% Calculate radial velocity
rad_v = dot(r,v)/rmag;

%% Calculate specific angular momentum
hcross = cross(r,v); %km^2/s

%% Calculate h magnitude
hmag = sqrt(dot(hcross,hcross)); %km^2/s verified

%% Calculate inclination
i = acosd(hcross(3)/hmag); %deg verified

%% Calculate vector defining the node line
Ncross = cross([0 0 1],hcross);

%% Calculate N magnitude
Nmag = norm(Ncross);

%% Calculate the right ascension of the ascending node
if Ncross(2)>=0
    Omega = acosd(Ncross(1)/Nmag); %deg Verified
elseif Ncross(2)<0
    Omega = 360-acosd(Ncross(1)/Nmag); %deg Verified
end

%% calculate eccentricity vector
evec = mu^(-1)*((vmag^2-mu/rmag).*r-rmag*rad_v.*v);

%% Calculate e mag
emag = norm(evec);%verified

%% Calculate argument of perigee
if evec(3)>=0
    lilomega = acosd(dot(Ncross,evec)/(Nmag*emag)); %deg verified
elseif evec(3)<0
    lilomega = 360 - acosd(dot(Ncross,evec)/(Nmag*emag)); %deg verified
end
%% Calculate true anomaly
if rad_v>=0
    theta = acosd(dot(evec,r)/(emag*rmag)); %deg verified
elseif rad_v<0
    theta = 360-acosd(dot(evec,r)/(emag*rmag)); %deg verified
end

%% Calculate orbital period
r_p = hmag^2/mu*(1/(1+emag));
r_a = hmag^2/mu*(1/(1-emag));

%semimajor axis
a = .5*(r_p+r_a); %km verified
%period
T = 2*pi/sqrt(mu)*a^(1.5); %sec verified

%% Load output structure
COEStruct.h_AngularMomentum = hmag; %km^2/s
COEStruct.e_Eccentricity = emag; 
COEStruct.i_Inclination = deg2rad(i);%radians
COEStruct.RA_RightAscension = deg2rad(Omega); %radians
COEStruct.w_ArgumentofPerigee = deg2rad(lilomega); %radians
COEStruct.nu_TrueAnomaly = deg2rad(theta);%radians
COEStruct.a_SemiMajorAxis = a;%km
COEStruct.T_Period = T; %sec
COEStruct.Rp_RadiusofPerigee = r_p; %km
COEStruct.Ra_RadiusofApogee = r_a; %km
COEStruct.mu = mu; %km^3/s^2
[COEStruct.M_MeanAnomaly, COEStruct.n_MeanMotion, COEStruct.E_EccAnom, ...
    COEStruct.t_TimeSincePer] = KeplerAnomaly(COEStruct);
end

