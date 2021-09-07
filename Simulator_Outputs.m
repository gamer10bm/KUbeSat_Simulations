%This script is used to generate output figures based on the outputs of the
%simulink model

%% Translate outputs
%Grab number of orbits simulated
Num_orb = 0;

%Grab vector of time points
time = 0;

%Determine number of timesteps per orbit
timesteps = 0;

%Grab generated power
inst_genpow = 0;

%Grab used power
inst_usepow = 0;

%Grab battery level
batt_pow = 0;

%Grab data storage
data_stored = 0;

%Grab torque build up
torque_build = 0;

%Grab lat and lon
lats = 0;
lons = 0;

%Send to SubPowStruct
SubPowStruct = 0;



%% Print Duty Cycles
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

%% Plot output figures
%Initialize figure saving structures
fig_save = struct('fignum','','savename','');
fignum = 0;

%Set time axis limits
period = 0;
shorttimelim = [0 2*period./60]; %minutes
longtimelim = [0 time(end)/60]; %minutes

%%%%%%%%%%% Instant Power vs. Time %%%%%%%%%%%%%
if 1
    fignum = fignum +1;
    %Send figure to save structure
    fig_save(end+1) = structure('fignum',fignum,'savename','KUbeSat1_Instant_Power');
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
    fig_save(end+1) = structure('fignum',fignum,'savename','KUbeSat1_Subsystem_Power');
    %Plot it
    figure(fignum)
    oneperlim = [0 period./60];
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
    fig_save(end+1) = structure('fignum',fignum,'savename','KUbeSat1_Battery_Storage');
    %Plot it
    figure(fignum)
    plot(time./60,batt_pow(1:end-1)./3600) %Whr
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
    fig_save(end+1) = structure('fignum',fignum,'savename','KUbeSat1_Data_Stored');
    %Plot it
    figure(fignum)
    plot(time./60,data_stored(1:end-1)./8)
    grid on
    title('Instantaneous Data Storage')
    ylabel('Data Stored (kB)')
    xlabel('Time (min)')
    xlim(longtimelim)
    FontWidthandPos
end

%%%%%%%%%% Momentum Buildup vs. Time %%%%%%%%%
if 1
    fignum = fignum +1;
    %Send figure to save structure
    fig_save(end+1) = structure('fignum',fignum,'savename','KUbeSat1_Instant_Power');
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
    fignum = fignum +1;
    %Send figure to save structure
    fig_save(end+1) = structure('fignum',fignum,'savename','KUbeSat1_Instant_Power');
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
save_dir = 'outputs';
save_type = '.fig';
for save_id = 2:length(fig_save)
    fn = fullfile(save_dir,[fig_save(save_id).savename save_type]);
    fprint('Saving Figure %d: %s\n',fig_save(save_id).fignum,fn)
    saveas(figure(fig_save(save_id).fignum),fn);
end