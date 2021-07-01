close all
clear
clc
%%% INPUTS %%%
% Constants %
R_e = 6378.14 ; % Radius of Earth (km)
solarFlux = 1376 ; % Max. Solar Flux (W/m^2)
earthIR_max = 258 ; % Max. Earth Infrared (W/m^2)
earthIR_min = 216 ; % Min. Earth Infrared (W/m^2)
albedo = 0.35 ; % Earth Albedo (~)
stefan = 5.67e-8 ; % Stefan-Boltzman Constant (W/m^2-K^4)
T_surr = 2.725 ; % Ambient Temperature of Space (K)
k = 237 ; % Thermal Conductivity of Aluminum @ 200-300 K (W/m-K)
c = 903 ; % Specific Heat of Aluminum @ 300 K (J/kg-K)
rho = 2702 ; % Density of Aluminum @ 300 K (kg/m^3)
% Design Choices %
h = 500 ; % Orbital Altitude (km)
sc_long = 0.3405 ; % Spacecraft length (m)
sc_wide = 0.10 ; % Spacecraft width (m)
Q_min = 0.02 ; % Minimum Power Draw (W)
Q_max = 5.8 ; % Maximum Power Draw (W)
epsIR = 0.8 ; % Infrared Emissivity (~)
alpha_metal = 0.316 ; % Solar Absorptivity of Degraded Aluminumized Teflon(~)
alpha_panel = 0.805 ; % Solar Absorptivity of Photovoltaic Panel (~)
mass = 3.3 ; % Spacecraft mass (kg)
min_power = 3.4 ; % Minimum power expenditure (W)
max_power = 4.3 ; % Maximum power expenditure (W)
%%% CALCULATIONS %%%
% Spacecraft Geometry %
V_sc = sc_wide*sc_wide*sc_long ; % Spacecraft Volume (m^3)
A_surf = 2*(sc_wide^2) + 4*(sc_long*sc_wide) ; % Spacecraft surface area (m^2)
D_sph = sqrt(A_surf/pi) ; % Equivalent sphere diameter (m)
V_sph = 4/3*pi*(D_sph/2)^3 ; % Equivalent sphere volume (m^3)
A_panel = 2*(sc_wide^2) + 2*(1.5*sc_wide^2) + 1*(sc_long*sc_wide) ; % Solar Panel Area (m^2)
A_other = A_surf - A_panel ; % Non-Solar Panel Surface Area (m^2)
alphaS = (A_panel*alpha_panel + A_other*alpha_metal)/(A_panel + A_other) ; % Area-Weighted Average Solar % Thermal Equilibrium %
Ka = 0.657 + 0.54*(R_e/(R_e+h)) - 0.196*(R_e/(R_e+h))^2 ; % Collimated Solar Radiation Factor (~)
F_se = 0.5*(1 - sqrt(h^2 + 2*h*R_e)/(h+R_e)) ; % Spacecraft-Earth View Factor (~)
Tmaxs = zeros(1,20) ; % Initialize worst-case hot temperature matrix
Tmins = zeros(1,20) ; % Initialize worst-case cold temperature matrix
Qs = linspace(min_power,max_power,20) ; % Power dissipation ranges
for i = 1:length(Qs)
Q = Qs(i) ;
T_max = ((0.25*solarFlux*alphaS + earthIR_max*epsIR*F_se + solarFlux*albedo*F_se*alphaS*Ka + (Q/A_surf))/(stefan*epsIR))^0.25;
T_min = ((earthIR_min*epsIR*F_se + (Q/A_surf))/(stefan*epsIR))^0.25 ;
Tmaxs(i) = T_max ;
Tmins(i) = T_min ;
end
% Lumped Capacitance Model %
L = D_sph/2 ; % Characteristic Length (m)
h_max = epsIR*stefan*(T_max+T_surr)*(T_max^2 + T_surr^2) ;
h_min = epsIR*stefan*(T_min+T_surr)*(T_min^2 + T_surr^2) ;
Bi_max = L*h_max/k ; % Biot Number at Hottest (~)
Bi_min = L*h_min/k ; % Biot Number at Coldest (~)
t_hc = (mass*c)/(h_max*A_surf) ; % Time to reach equilibrium passing into shadow (s)
t_ch = (mass*c)/(h_min*A_surf) ; % Time to reach equilibrium passing into sunlight (s)
%%% OUTPUTS %%%
fprintf('Max. Biot Number: %g \n',Bi_max)
fprintf('Min. Biot Number: %g \n',Bi_min)
fprintf('LCM Time-to-Equilibrium (entering eclipse): %g hr\n',t_hc/3600)
fprintf('LCM Time-to-Equilibrium (leaving shadow): %g hr\n',t_ch/3600)
subplot(1,2,1)
hold on
plot(Qs,Tmins,'linewidth',2)
plot(Qs,Tmaxs,'linewidth',2)
xlabel('Internal Power Dissipation (W)')
ylabel('Steady-State Equilibrium Temperature (K)')
legend('Worst-Case Cold','Worst-Case Hot','location','best')
subplot(1,2,2)
hold on
plot(Qs,Tmins-273,'linewidth',2)
plot(Qs,Tmaxs-273,'linewidth',2)
xlabel('Internal Power Dissipation (W)')
ylabel('Steady-State Equilibrium Temperature (^oC)')
legend('Worst-Case Cold','Worst-Case Hot','location','best')

function [F_se] = sphereViewFactor(R,h)
%SPHEREVIEWFACTOR Computes the sphere-to-sphere view factor for a satellite
%in orbit of a planet.
% R - radius of the planet (km)
% h - altitude above the planet (km)
F_se = 0.5*(1 - sqrt(h^2 + 2*h*R)/(h+R)) ; % Spacecraft-Earth View Factor (~)
end

function [Ka] = collimatingFactorh(R_e,h)
%COLLIMATINGFACTOR Computed the Collimating Factor for a satellite in orbit
%of Earth
% R_e = radius of Earth (km)
% h = altitude above the planet (km)
Ka = 0.657 + 0.54*(R_e/(R_e+h)) - 0.196*(R_e/(R_e+h))^2 ; % Collimated Solar Radiation Factor (~)
end