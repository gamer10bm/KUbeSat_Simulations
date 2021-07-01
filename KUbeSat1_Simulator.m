clc
close all

addpath('KUbeSat_Subfunctions')
addpath('3D_Shape')

%% KUbeSat1 Matlab Simulations
%Define settings for simulation
perigee_altitude = 550; %km
N =2; %Number of orbits
timesteps = 50; %Time steps per orbit
yawsteps = 100; %Yaw steps per time step
RAANinit = 100.0330; %deg From STK
% RAANinit = 360+lon_gs-97; %deg First pass is over lawrence


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

%% Initialization Steps
% initialize vehicle model
SC = initiate_SC_model();

% initialize simulation epoch
yr_init = 2022; mnth_init = 1; day_init = 1; hr_init = 0;
min_init = 0; sec_init = 0;
init_utcvec = [yr_init, mnth_init, day_init, hr_init, min_init, sec_init];
epoch = initial_epoch(yr_init, mnth_init, day_init, hr_init, min_init, sec_init);

% determine SSO orbit for given altitude and spacecraft
OE0 = SSO_Earth(perigee_altitude,SC);
OE0(4) = (RAANinit)*pi/180; % Sets RAAN and centers groundtrack on Lawrence
[~,Rstart,Vstart] = coe2RV(OE0(1),OE0(2),OE0(3),OE0(6),OE0(4),OE0(5),SC.mu);
COEstruct0 = RV2coe(Rstart,Vstart,SC.mu);

% compute orbit period and initial mean anomaly
T = COEstruct0.T_Period;

% determine time array parameters
n = timesteps*N; % Basing array on number of orbits and points per orbit
time = linspace(0, N*T, n); % define time array
dt = time(2) - time(1); % save time step

% identify the yaw angle array (eliminate redundant boundary)
n2 = yawsteps;
Yaw = linspace(0, 2*pi, n2+1); Yaw = Yaw(1:n2);

% initialize vehicle state structure
 if exist('state','var')
    clearvars state
end
state(n) = initiate_SC_state(SC);

% preallocate states structure memory
for i = n-1:-1:1
    state(i) = state(n);
end

%% Do Keplerian Motion Simulation
COEstructnow = COEstruct0;
for i = 1:n
    % store time from array
    state(i).t = time(i);
    
    % find and store vehicle translation states    
    [~,COEstructnow] = FutureAnomaly(state(i).t,COEstructnow);
    state(i).OE = COEstructnow;
    [~,R, V] = coe2RV(COEstructnow);
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
    [state(i).lat, state(i).lon] =ECI2latlon(R_IJK,state(i).t,init_utcvec);
end

%% Do solar panel generation 
%Iterate over generated states to determine power generation from solar
%panels
inst_genpow = zeros(n,1);
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
    inst_genpow(i) = SP_max; %W
end

%% Momentum buildup
inst_torque = zeros(1,n);
for i = 1:n
    R = state(i).R; V = state(i).V;
    eclip_chk = false; %Default in sunlight
    if inst_genpow(i) <= 0
        eclip_chk = true; %Change to in eclipse
    end
%     inst_torque(i) = momentumbuildup(SC,R,V,eclip_chk); %N-m
    inst_torque(i) = momentum_Daniel(SC,R,V,epoch, state(i).t); %N-m
end

%% Power usage and data stored
[inst_usepow, batt_pow, data_stored, data_rate, torque_build, SubPowStruct] = ...
    KUbeSat1_PowerandData(SC,state,inst_torque,inst_genpow);

%Get and print duty cycles from SubPowstruct
fnames = fieldnames(SubPowStruct);
fduty = zeros(1,length(fnames)-1);
for fid = 2:length(fnames)
    %Get the maximum value for the field
    fmax = max(SubPowStruct.(fnames{fid}));
    %Get number of values at max
    fmax_occurs = length(find(SubPowStruct.(fnames{fid})==fmax));
    %Get duty cycle based on total SubPow length
    fduty(fid) = fmax_occurs./length(SubPowStruct.(fnames{fid}));
    %Print it
    fprintf('\n%s duty cycle = %.2f%s\n',fnames{fid},fduty(fid)*100,'%')
end
%% Plot Section
shorttimelim = [0 2*T./60];
longtimelim = [0 state(end).t/60];
%%%%%%%%%%% Instant Power vs. Time %%%%%%%%%%%%%
figure(3)
gp = plot(time./60,inst_genpow);
hold on
plot([time(1) time(end)]./60,[mean(inst_genpow) mean(inst_genpow)],':','Color',gp.Color)
ip = plot(time./60,inst_usepow);
plot([time(1) time(end)]./60,[mean(inst_usepow) mean(inst_usepow)],':','Color',ip.Color)
hold off
grid on
legend({'Generated','Avg. Generated','Used','Avg. Used'})
title('Instantaneous Power')
xlabel('Time (min)')
xlim(shorttimelim)
ylabel('Power (W)')
FontWidthandPos

%%%%%%%% Subsystem Power vs. Time %%%%%%%%%%%%%
oneperlim = [0 T./60];
subfields = fieldnames(SubPowStruct);
% subfields = subfields(2:end);
figure(10)
for ids = 1:length(subfields)
    plot(time./60,SubPowStruct.(subfields{ids}),'-.')
    hold on
end
hold off
grid on
legend(subfields)
title('Instant Power by Subsystem')
xlabel('Time (min)')
ylabel('Power (W)')
xlim(oneperlim)
FontWidthandPos

%%%%%%%%%%% Battery Storage vs. Time %%%%%%%%%%
figure(4)
plot(time./60,batt_pow(1:end-1)./3600) %Whr
grid on
title('Battery Power')
xlabel('Time (min)')
xlim(longtimelim)
ylabel('Power (Whr)')
FontWidthandPos

%%%%%%%%%% Data Storage vs. Time %%%%%%%%%%%
figure(2)
plot(time./60,data_stored(1:end-1)./8)
grid on
title('Instantaneous Data Storage')
ylabel('Data Stored (kB)')
xlabel('Time (min)')
xlim(longtimelim)
FontWidthandPos

%%%%%%%%%% Momentum Buildup vs. Time %%%%%%%%%
figure(5)
plot(time./60,torque_build(1:end-1)./1e3)
grid on
title('Reaction Wheel Momentum Buildup')
xlabel('Time (min)')
xlim(shorttimelim)
ylabel('Momentum (N-m-sec)')
FontWidthandPos

%%%%%%%%%%%%% Ground Track %%%%%%%%%%%%%%%%
statecell = struct2cell(state); statefields = fieldnames(state);
latname = 'lat';
latid = find(strcmp(statefields,latname));
lonid = latid+1;
lats = reshape(cell2mat(statecell(latid,:,:)),[1 n]);
lons = reshape(cell2mat(statecell(lonid,:,:)),[1 n]);
%Plot the ground track
N_off = 2; %Orbit offsets in the groundtrack plot
N_now = 1;
minichunk = 1:timesteps;
chunk = [];
while N_now <=N
    if N_now == 1 || rem(N_now,N_off)==0
        chunk = [chunk minichunk+(N_now-1)*timesteps];
    end
    N_now = N_now+1;
end
figure(1)
itvec = 1:10:length(chunk);
if itvec(end) ~= length(chunk)
    itvec(end+1) = length(chunk);
end
for i = length(chunk)% 1:10:length(chunk)
    clf
    subchunk = chunk(1:i);
    geoaxes();
    geoscatter(lats(subchunk),lons(subchunk),'.b')
    hold on
    geoscatter(telemcirc_lats,telemcirc_lons,'.r')
    geoscatter(HiCalXcirc_lats,HiCalXcirc_lons,'.m')
    geolimits([-90 90], [-180 180])
    title('Ground Track')
    % title(sprintf('Ground Track for %.0f Orbits',N))
    FontWidthandPos
    drawnow
end