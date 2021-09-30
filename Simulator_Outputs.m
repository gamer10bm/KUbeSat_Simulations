close all
clearvars
%This script is used to generate output figures based on the outputs of the
%simulink model

addpath('KUbeSat_Subfunctions')
addpath('3D_Shape')
addpath('SimV2')

%% Run Simulink Model with below settings
Num_orb = 1;
Period = 95.65; %min (approximately)
time_step = 5; %sec 
perigee_altitude = 550; %km %Change to get from R values

time_sim = Num_orb*Period*60; %sec
model_name = fullfile('SimV2','SatStates');

%Run simulation
simout = sim(model_name,'StartTime','0','StopTime',num2str(time_sim),'FixedStep',num2str(time_step));
%% Translate outputs

%Grab vector of time points
time = simout.tout; %seconds

%Determine number of timesteps per orbit
timesteps = length(time);

%Grab generated power
inst_genpow = simout.Pwr_SP.Data;

%Grab battery level
batt_pow = simout.Batt_Charge.Data;

%Grab data storage
data_stored = simout.data.Data;

%Grab torque build up
torque_build = 0;

%% Generate lat and lon
R = simout.R.Data;
V = simout.V.Data;
% initialize simulation epoch
yr_init = 2022; mnth_init = 1; day_init = 1; hr_init = 0;
min_init = 0; sec_init = 0;
init_utcvec = [yr_init, mnth_init, day_init, hr_init, min_init, sec_init];
lats = zeros(1,size(R,1)); lons = zeros(1,size(R,1));
for i = 1:length(R)
    
    [lats(i), lons(i)] =ECI2latlon(R(i,:)',time(i),init_utcvec);
end

%Get list of output fields
outfields = simout.who;
%Find valid fields for power
pow_pat = 'Pwr_';
clearvars SubPowStruct
for id = 1:length(outfields)
    curfield = outfields{id};
    if contains(curfield,pow_pat) && ~contains(curfield,'SP')
        split_str = strsplit(curfield,pow_pat);
        sub_str = split_str{2}; %2 is the subsystem name;
        %Load to SubPowStruct
        SubPowStruct.(sub_str) = simout.(curfield).Data;
    end
end

%% Print Duty Cycles and Generate Instant Use Power
%Get and print duty cycles from SubPowstruct
fnames = fieldnames(SubPowStruct);
fduty = zeros(1,length(fnames)-1);
inst_usepow = [];
for fid = 1:length(fnames)
    %Get the maximum value for the field
    fmax = max(SubPowStruct.(fnames{fid}));
    %Get number of values at max
    fmax_occurs = length(find(SubPowStruct.(fnames{fid})==fmax));
    %Get duty cycle based on total SubPow length
    fduty(fid) = fmax_occurs./length(SubPowStruct.(fnames{fid}));
    %Print it
    fprintf('\n%s duty cycle = %.2f%s\n',fnames{fid},fduty(fid)*100,'%')
    %Save to inst_usepow
    if fid == 1 || isempty(inst_usepow)
        inst_usepow = SubPowStruct.(fnames{fid});
    else
        inst_usepow = inst_usepow+SubPowStruct.(fnames{fid});
    end
end

%% Plot output figures
%Initialize figure saving structures
fig_save = struct('fignum','','savename','');
fignum = 0;

%Set time axis limits
shorttimelim = [0 2*Period./60]; %minutes
longtimelim = [0 time(end)/60]; %minutes

%%%%%%%%%%% Instant Power vs. Time %%%%%%%%%%%%%
if 1
    fignum = fignum +1;
    %Send figure to save structure
    fig_save(end+1) = struct('fignum',fignum,'savename','KUbeSat1_Instant_Power');
    %Plot it
    figure(fignum)
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
end

%%%%%%%% Subsystem Power vs. Time %%%%%%%%%%%%%
if 1
    fignum = fignum +1;
    %Send figure to save structure
    fig_save(end+1) = struct('fignum',fignum,'savename','KUbeSat1_Subsystem_Power');
    %Plot it
    figure(fignum)
    oneperlim = [0 Period./60];
    subfields = fieldnames(SubPowStruct);
    % subfields = subfields(2:end);
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
end

%%%%%%%%%%% Battery Storage vs. Time %%%%%%%%%%
if 1
    fignum = fignum +1;
    %Send figure to save structure
    fig_save(end+1) = struct('fignum',fignum,'savename','KUbeSat1_Battery_Storage');
    %Plot it
    figure(fignum)
    plot(time./60,batt_pow./3600) %Whr
    grid on
    title('Battery Power')
    xlabel('Time (min)')
    xlim(longtimelim)
    ylabel('Power (Whr)')
    FontWidthandPos
end

%%%%%%%%%% Data Storage vs. Time %%%%%%%%%%%
if 1
    fignum = fignum +1;
    %Send figure to save structure
    fig_save(end+1) = struct('fignum',fignum,'savename','KUbeSat1_Data_Stored');
    %Plot it
    figure(fignum)
    plot(time./60,data_stored./8)
    grid on
    title('Instantaneous Data Storage')
    ylabel('Data Stored (kB)')
    xlabel('Time (min)')
    xlim(longtimelim)
    FontWidthandPos
end

%%%%%%%%%% Momentum Buildup vs. Time %%%%%%%%%
if 0
    fignum = fignum +1;
    %Send figure to save structure
    fig_save(end+1) = struct('fignum',fignum,'savename','KUbeSat1_Instant_Power');
    %Plot it
    figure(fignum)
    plot(time./60,torque_build(1:end-1)./1e3)
    grid on
    title('Reaction Wheel Momentum Buildup')
    xlabel('Time (min)')
    xlim(shorttimelim)
    ylabel('Momentum (N-m-sec)')
    FontWidthandPos
end

%%%%%%%%%%%%% Ground Track %%%%%%%%%%%%%%%%
if 1
    %Generate geocircles for groundstation and mcmurdo
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
    %Plot it
    fignum = fignum +1;
    %Send figure to save structure
    fig_save(end+1) = struct('fignum',fignum,'savename','KUbeSat1_Instant_Power');
    %Plot it
    figure(fignum)
    %Plot the ground track
    N_off = 2; %Orbit offsets in the groundtrack plot
    N_now = 1;
    minichunk = 1:timesteps;
    chunk = [];
    while N_now <=Num_orb
        if N_now == 1 || rem(N_now,N_off)==0
            chunk = [chunk minichunk+(N_now-1)*timesteps];
        end
        N_now = N_now+1;
    end
    itvec = 1:10:length(chunk);
    if itvec(end) ~= length(chunk)
        itvec(end+1) = length(chunk);
    end
    %Can be made to make a movie
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
end

%% Save figures
if 0
    save_dir = 'outputs';
    save_type = '.fig';
    for save_id = 2:length(fig_save)
        fn = fullfile(save_dir,[fig_save(save_id).savename save_type]);
        fprint('Saving Figure %d: %s\n',fig_save(save_id).fignum,fn)
        saveas(figure(fig_save(save_id).fignum),fn);
    end
end