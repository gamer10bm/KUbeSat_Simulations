clc
clearvars
close all

addpath('KUbeSat_Subfunctions')

%% KUbeSat1 Matlab Simulations
%Define settings for simulation
perigee_altitude = 550; %km
N =3; %Number of orbits
timesteps = 100; %Time steps per orbit
yawsteps = 100; %Yaw steps per time step
raansteps = 360; %Right Ascension steps
ideal_raan = -89.0352; %Right ascension used in simulator
if ideal_raan <0
    ideal_raan = 360+ideal_raan;
end

%% Initialization Steps
% initialize vehicle model
SC = initiate_SC_model();

% initialize simulation epoch
epoch = initial_epoch(2022, 1, 1, 0, 0, 0);

% determine SSO orbit for given altitude and spacecraft
OE0 = SSO_Earth(perigee_altitude,SC);

% compute orbit period and initial mean anomaly
p = OE0(1)^2/SC.mu; T = 2*pi/sqrt(SC.mu)*p^(3/2);
E0 = theta_to_E(OE0(2),OE0(6)); M0 = E0 - OE0(2)*sin(E0);

% determine time array parameters
n = timesteps*N; % Basing array on number of orbits and points per orbit
time = linspace(0, N*T, n); % define time array
dt = time(2) - time(1); % save time step

% identify the yaw angle array (eliminate redundant boundary)
n2 = yawsteps;
Yaw = linspace(0, 2*pi, n2+1); Yaw = Yaw(1:n2);

% Initialize RAAN vector
RAANmat = linspace(0, 360,raansteps); %degrees

% Initialize generated power matrix (also zdat)
inst_genpow = zeros(n,raansteps);

%% Iterate over RAANs
refstate = []; refraan = [];
for rid = 1:raansteps
    %Set the new RAAN
    OE0(4) = (RAANmat(rid))*pi/180; % Sets RAAN
    % initialize vehicle state structure
    if exist('state','var')
        clearvars state
    end
    state(n) = initiate_SC_state(SC);

    % preallocate states structure memory
    % This is necessary because although we told Matlab it would have a certain
    % number of elements, the actual state elements are stored as pointers. So,
    % we need to assign values to these elements to hold the data memory.
    % Matlab might still throw a warning when accessing structure elements in 
    % loops but the memory allocation should be fine to reduce perfomance hit.
    for i = n-1:-1:1
        state(i) = state(n);
    end

    %% Do Keplerian Motion Simulation
    for i = 1:n
        % store time from array
        state(i).t = time(i);

        % find and store vehicle translation states
        M = mod(M0 + 2*pi*state(i).t/T,2*pi);
        E = kepler_ellipse_L(OE0(2), M);
        TA = E_to_theta(OE0(2),E);
        OE = OE0; OE(6) = TA;
        state(i).OE = OE;
        [R, V] = OE2SV(SC.mu, OE);
        state(i).R = R; state(i).V = V;

        % current routines do not use a separate estimation state but save
        % it anyway in case we do in the future
        state(i).R_est = R; state(i).V_est = V;
    end

    %% Do Runga-Kutta Simulation with J2 and Drag
    %Initial conditions
    rstart =state(1).R; %km
    vstart = state(1).V; %km/s

    Xstart = [rstart;vstart];

    dervfunc = @(t,X)OrbitDerivFunc(X,SC.Re,SC.mu,SC.J2,SC.CD,SC.A,SC.m,SC.rho0,SC.r0,SC.H,SC.thetadot);
    %Do the Runge-Kutta
    tstart = 0; tend = N*T;
    [tRK,XRK]=RungeKutta(dervfunc,Xstart,dt,tstart,tend+dt);

    %Put RK results into state matrix
    for i = 1:n
        Rnow = XRK(1:3,i); Vnow = XRK(4:6,i);
        state(i).R = Rnow;
        state(i).V = Vnow;
    end

    for i = 1:n
        %Determine latitude and longitude (no earth motion)
        R_IJK = state(i).R;
        [state(i).lat,lon_nomo] = ECEF2latlon(R_IJK);

        %Change longitude based on earth's rotation
        P_s = 24*3600; %seconds length of mean solar day
        P_e = 86164; %seconds rotation period of earth
        lon_mo = lon_nomo+360*(1/P_e-1/P_s)*state(i).t;
        if lon_mo>=180
            lon_mo = lon_mo-360*floor(i/n*N);
        end
        state(i).lon = lon_mo;
    end
    
    %Grab ideal raan state
    if isempty(refraan) || abs(refraan-ideal_raan)>abs(RAANmat(rid)-ideal_raan)
        %Store raan and state
        refstate = state;
        refraan = RAANmat(rid);
    end

    %% Do solar panel generation 
    %Iterate over generated states to determine power generation from solar
    %panels
    
    % We need to set the reference value for the quaternion conversion
    % routine. In the initial step, set this to identity [0,0,0,1].
    % Otherwise, set this to the value that produced the optimum last time.
    Q_near = [0; 0; 0; 1];
    for i = 1:n
        % find sun direction at this time
        RS = sun_vector(epoch, state(i).t);

        % initialize maximum solar power
        SP_max = 0;

        % cycle through possible yaw angles
        for j = 1:n2
            % determine applicable quaternion attitude based on yaw angle from
            % array and set pitch and roll to zero. This is based on the local
            % RTN coordinate system, so we need position and velocity vectors
            %
            % This routine takes in a reference quaternion to pick the branch
            % since our possible answer are non-unique. Using the most recent
            % value for quaternion attitude keeps this continuous.
            Q = quaternionC(0, 0, Yaw(j), state(i).R, state(i).V, Q_near);

            % predict solar collection at this yaw angle
            SP = solar_panel(SC, state(i).R, Q, RS);

            % check if maximum and save
            if (SP > SP_max)
                SP_max = SP ;
                Qmax = Q;
            end
        end

        % save maximum power quaternion to state array attitude command
        if (SP_max > eps) 
            state(i).Q_com = Qmax; 
        else
            state(i).Q_com = Q_near;
        end

        % set maximum value to reference quaternion if no eclipse
        if (SP_max > eps) 
            Q_near = Qmax;
        end

        % save maximum power generated
        inst_genpow(i,rid) = SP_max;
    end
end

%% Plot RAAN as heatmap
zdat = inst_genpow'; %2D matrix of Solar power vs time vs. RAAN
ydat = RAANmat; %Array of right ascension values
xdat = time./60; %Array of time values
figure(1)
contourf(xdat,ydat,zdat);
xlabel('Time (min)')
ylabel('Right Ascension \Omega (deg)')
c = colorbar;
c.Label.String = 'Solar Power Generated (W)';
title(sprintf('Generated Power over %.0f Orbit Sweeping %s',N,'\Omega'))
FontWidthandPos

%% Plot best and worst case of solar power generation
meanpow = mean(inst_genpow);
meanmeanpow = mean(meanpow);
[bestcase_pow, best_id] = max(meanpow);
[worstcase_pow, worst_id] = min(meanpow);

figure(2)
bplot = plot(xdat,inst_genpow(:,best_id));
hold on
plot([xdat(1) xdat(end)], [bestcase_pow bestcase_pow],':','Color',bplot.Color)
wplot = plot(xdat,inst_genpow(:,worst_id));
plot([xdat(1) xdat(end)], [worstcase_pow worstcase_pow],':','Color',wplot.Color)
plot([xdat(1) xdat(end)], [meanmeanpow meanmeanpow],':')
avg_consum = 5.5176; %W from KUbeSat1_Simulator
plot([xdat(1) xdat(end)], [avg_consum avg_consum])
hold off
grid on
legend({'Best Profile','Best Average','Worst Profile','Worst Average','Average Average','Consumed Power Average'})
ylabel('Generated Power (W)')
xlabel('Time (min)')
title(sprintf('Best (%s=%.0f^o) and Worst (%s=%.0f^o) Solar Profiles','\Omega',RAANmat(best_id),'\Omega',RAANmat(worst_id)))
FontWidthandPos

%% Plot best and worst case battery storage based on ideal_raan state
%Get inst_torque for ref state
% Momentum buildup
inst_torque = zeros(1,n);
for i = 1:n
    R = refstate(i).R; V = refstate(i).V;
    eclip_chk = false; %Default in sunlight
    if inst_genpow(i) <= 0
        eclip_chk = true; %Change to in eclipse
    end
    inst_torque(i) = momentumbuildup(SC,R,V,eclip_chk); %N-m
end
%Get the best and worst case battery stored
[inst_usepow_best, batt_pow_best,~, ~, ~, SubPowStruct_best] = ...
KUbeSat1_PowerandData(SC,refstate,inst_torque,inst_genpow(:,best_id));
[inst_usepow_worst, batt_pow_worst,~, ~, ~, SubPowStruct_worst] = ...
KUbeSat1_PowerandData(SC,refstate,inst_torque,inst_genpow(:,worst_id));
%Get average instant use power
avg_consum_best = mean(inst_usepow_best);
avg_consum_worst = mean(inst_usepow_best);
%Make battery state plot
figure(20)
plot(time./60,batt_pow_best(1:end-1)./3600) %Whr
hold on
plot(time./60,batt_pow_worst(1:end-1)./3600) %Whr
grid on
title('Battery Power')
xlabel('Time (min)')
% xlim(shorttimelim)
ylabel('Power (Whr)')
legend({'Best Solar Profile','Worst Solar Profile'})
FontWidthandPos
%Make battery state plot with two y axes
if 1 %Probably to busy to present on
    figure(22)
    yyaxis('left')
    plot(time./60,batt_pow_best(1:end-1)./3600) %Whr
    hold on
    plot(time./60,batt_pow_worst(1:end-1)./3600) %Whr
    hold off
    grid on
    ylabel('Energy Stored (Whr)')
    yyaxis('right')
    plot(time./60,inst_usepow_best,'r--')
    hold on
    plot(time./60,inst_usepow_worst,'k--')
    hold off
    grid on
    ylabel('Instant Power Draw (W)')
    title('Battery Power')
    xlabel('Time (min)')
    legend({'Best Solar Profile','Worst Solar Profile','Best Use Power','Worst Use Power'})
    FontWidthandPos
end

%Plot groundtrack as chk
if 0
    %Make geocircles
    lat_gs = 38.971669; %(deg) [Lawrence, KS]
    lon_gs = -95.23525; %(deg) [Lawrence, KS]
    gs_ant_BW = 131; %Groundstation Antenna Beamwidth (deg)
    capture_radius = perigee_altitude*tand(gs_ant_BW/2);
    %Get geocircle for transmit
    [telemcirc_lats, telemcirc_lons] = geocircle(lat_gs,lon_gs,capture_radius);

    lat_mcmurdo = -82.406; lon_mcmurdo = 0; %(deg) McMurdo Station
    mcmurdo_ant_BW = 90;
    capture_radius = perigee_altitude*tand(mcmurdo_ant_BW/2);
    %Get geocircle for HiCalX pulses
    [HiCalXcirc_lats, HiCalXcirc_lons] = geocircle(lat_mcmurdo,lon_mcmurdo,capture_radius);
    %Pull lat and lon 
    statecell = struct2cell(refstate); statefields = fieldnames(refstate);
    latname = 'lat';
    latid = find(strcmp(statefields,latname));
    lonid = latid+1;
    lats = reshape(cell2mat(statecell(latid,:,:)),[1 n]);
    lons = reshape(cell2mat(statecell(lonid,:,:)),[1 n]);
    %Plot ground track of refstate
    chunk = 1:length(refstate);
    figure(30)
    geoaxes();
    geoscatter(lats(chunk),lons(chunk),'.b')
    hold on
    geoscatter(telemcirc_lats,telemcirc_lons,'.r')
    geoscatter(HiCalXcirc_lats,HiCalXcirc_lons,'.m')
    geolimits([-90 90], [-180 180])
    title('Ground Track')
    % title(sprintf('Ground Track for %.0f Orbits',N))
    FontWidthandPos
end