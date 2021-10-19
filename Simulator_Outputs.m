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

base_dir = 'SimV2';
model_name = 'SatStates';
model_path = fullfile(base_dir,model_name);
if 1
    %Load custom simulator settings (open file for examples)
    Simulator_Settings
end

%Run simulation
simout = sim(model_path,'StartTime','0','StopTime',num2str(time_sim),'FixedStep',num2str(time_step));
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

%Grab lat and lon in model
lats2 = simout.lat.Data;
lons2 = simout.lon.Data;

%Grab torque build up
torque_build = 0;

%Grab quaternions
quats = cell(1,length(time));
for iq = 1:length(quats)
    quats{iq} = [simout.quaternion.Data(iq,1); simout.quaternion.Data(iq,2);...
        simout.quaternion.Data(iq,3); simout.quaternion.Data(iq,4)];
end

%% Generate lat and lon
R = simout.R.Data;
V = simout.V.Data;
% initialize simulation epoch
yr_init = 2022; mnth_init = 1; day_init = 1; hr_init = 0;
min_init = 0; sec_init = 0;
init_utcvec = [yr_init, mnth_init, day_init, hr_init, min_init, sec_init];
lats = zeros(1,size(R,1)); lons = zeros(1,size(R,1));
for i = 1:length(R)
    
    [lats(i), lons(i)] =ECEF2latlon(R(i,:)',time(i));
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
shorttimelim = [0 time(end)/60]; %minutes
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
    oneperlim = [0 Period];
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
    plot(time./60,batt_pow) %Whr
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

%     lat_mcmurdo = -82.406; lon_mcmurdo = 0; %(deg) McMurdo Station
%     mcmurdo_ant_BW = 90;
%     capture_radius = perigee_altitude*tand(mcmurdo_ant_BW/2);
%     %Get geocircle for HiCalX pulses
%     [HiCalXcirc_lats, HiCalXcirc_lons] = geocircle(lat_mcmurdo,lon_mcmurdo,capture_radius);
    if exist('const_struc','var')
        const_cells = struct2cell(const_struc);
        HiCal_ind = find(strcmp(const_cells(1,1,:),'Lat_HiCal_deg'),1,'first');
    else
        HiCal_ind = [];
    end
    if isempty(HiCal_ind)
        HiCal_lat_deg = -60;
    else
        HiCal_lat_deg = str2num(const_struc(HiCal_ind).val);
    end
    HiCalXcirc_lats = ones(1,100)*HiCal_lat_deg;
    HiCalXcirc_lons = linspace(-180, 180,length(HiCalXcirc_lats));
    %Plot it
    fignum = fignum +1;
    %Send figure to save structure
    fig_save(end+1) = struct('fignum',fignum,'savename','KUbeSat1_Instant_Power');
    %Plot it
    figure(fignum)
    %Plot the ground track
    N_off = 2; %Orbit offsets in the groundtrack plot
    N_now = 1;
    chunk = 1:timesteps;
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
        if exist('lats2','var') && ~isempty(lats2)
            geoscatter(lats2(subchunk),lons2(subchunk),'.g')
        end
        geoscatter(telemcirc_lats,telemcirc_lons,'.r')
        geoscatter(HiCalXcirc_lats,HiCalXcirc_lons,'.m')
        geolimits([-90 90], [-180 180])
        title('Ground Track')
        % title(sprintf('Ground Track for %.0f Orbits',N))
        FontWidthandPos
        drawnow
    end
end

%%%%%%%%%%%%% Orientation vs. Time %%%%%%%%%%%%%%%%
if 1
    %% Make 3d matrix representing the KUbeSat
    globalmax_m = 1;%m
    globalN = 300;

    SC_sizevec_m = [0.3 0.1 0.1]; %m
    SC_strtvec_m = [0.1 0.1 0.1]; %m

    [globalmat_SC, spacevec_m] = make3drect(SC_sizevec_m, SC_strtvec_m, globalmax_m,globalN);

    %Cut hole through middle with some thickness on either side
    wall_thick = 0.01; %m
    yz_sub_sizevec_m = SC_sizevec_m-wall_thick*2;
    yz_sub_sizevec_m(1) = globalmax_m; %Cuts all the way through x-axis
    yz_sub_strtvec_m = SC_strtvec_m+wall_thick;
    yz_sub_strtvec_m(1) = 0; 
    [globalmat_yz_sub] = make3drect(yz_sub_sizevec_m,yz_sub_strtvec_m,globalmax_m,globalN);

    %Make 2 xy cuts 
    x_off = wall_thick/2; y_off = wall_thick;
    xy_inner_sub_sizevec_m = [SC_sizevec_m(1)./2-x_off*2 SC_sizevec_m(2)-y_off*2 0];
    xy_inner_sub_sizevec_m(3) = globalmax_m;
    xy1_inner_sub_strtvec_m = SC_strtvec_m +[x_off y_off 0];
    xy1_inner_sub_strtvec_m(3) = 0;
    xy2_inner_sub_strtvec_m = xy1_inner_sub_strtvec_m;
    xy2_inner_sub_strtvec_m(1) = xy1_inner_sub_strtvec_m(1) + SC_sizevec_m(1)/2;
    xy2_inner_sub_strtvec_m(3) =0;
    [globalmat_xy1_inner_sub] = make3drect(xy_inner_sub_sizevec_m,xy1_inner_sub_strtvec_m,globalmax_m,globalN);
    [globalmat_xy2_inner_sub] = make3drect(xy_inner_sub_sizevec_m,xy2_inner_sub_strtvec_m,globalmax_m,globalN);


    %Cut hole in globalmat with subtraction mats
    globalmat = globalmat_SC;
    globalmat = and(globalmat,~globalmat_xy1_inner_sub);
    globalmat = and(globalmat,~globalmat_xy2_inner_sub);
    globalmat = and(globalmat,~globalmat_yz_sub);

    %% Extract external geometry and center about origin
    [XYfaces, XZfaces, YZfaces, cg_xyz] = extract3dshape(globalmat,spacevec_m,1);

    %Center each face
    XYfaces_cent(1,:) = XYfaces(1,:)-cg_xyz(1);
    XYfaces_cent(2,:) = XYfaces(2,:)-cg_xyz(2);
    XYfaces_cent(3,:) = XYfaces(3,:)-cg_xyz(3);
    XZfaces_cent(1,:) = XZfaces(1,:)-cg_xyz(1);
    XZfaces_cent(2,:) = XZfaces(2,:)-cg_xyz(2);
    XZfaces_cent(3,:) = XZfaces(3,:)-cg_xyz(3);
    YZfaces_cent(1,:) = YZfaces(1,:)-cg_xyz(1);
    YZfaces_cent(2,:) = YZfaces(2,:)-cg_xyz(2);
    YZfaces_cent(3,:) = YZfaces(3,:)-cg_xyz(3);


    %% Plot the result 
    %Initialize for movie
    moviestep = 60; %sec
    itvec = floor(1:moviestep/time_step:length(time));
    timestoit = time(itvec);

    %Get rotations from statestoit
    roll = zeros(1,length(timestoit));
    pitch = zeros(1,length(timestoit));
    yaw = zeros(1,length(timestoit));
    Rotmat = cell(1,length(timestoit));
    for i = 1:length(timestoit)
        %Get the rotation matrix and 
        [roll(i), pitch(i), yaw(i), Rotmat{i}] = quaternion2euler(quats{i});
    end
    keyboard
    %Make a moooovie
    fignum = fignum + 1;
    h = figure(fignum); clf
    axsat = axes('NextPlot','replacechildren','OuterPosition',[0 0 .8 1]);
    ax2 = axes('Box','on','OuterPosition',[0.6 .4 .4 .6],'NextPlot','replacechildren');
    M(length(timestoit)) = struct('cdata',[],'colormap',[]);
    for i = 1:length(timestoit)
        %Rotate all the faces
        XZrot = Rotmat{i}*XZfaces_cent;
        XYrot = Rotmat{i}*XYfaces_cent;
        YZrot = Rotmat{i}*YZfaces_cent;
        set(h,'currentaxes',axsat);
        cla;
        scatter3(axsat,XZrot(1,:),XZrot(2,:),XZrot(3,:),'filled')
        hold(axsat,'on');
        scatter3(axsat,XYrot(1,:),XYrot(2,:),XYrot(3,:),'filled')
        scatter3(axsat,YZrot(1,:),YZrot(2,:),YZrot(3,:),'filled')
        %Plot axis for reference
        n = 50;
        ax1vec = linspace(0,globalmax_m/2,n);
        ax23vec = zeros(1,n);
        scatter3(axsat,ax1vec,ax23vec,ax23vec,12,'ko','filled')
        scatter3(axsat,ax23vec,ax1vec,ax23vec,12,'ko','filled')
        scatter3(axsat,ax23vec,ax23vec,ax1vec,12,'ko','filled')
        n = n/2;
        ax1vec = linspace(0,-globalmax_m/2,n);
        ax23vec = zeros(1,n);
        scatter3(axsat,ax1vec,ax23vec,ax23vec,12,'ks','filled')
        scatter3(axsat,ax23vec,ax1vec,ax23vec,12,'ks','filled')
        scatter3(axsat,ax23vec,ax23vec,ax1vec,12,'ks','filled')
        hold(axsat,'off');
        grid on
        % xlim([-globalmax_m globalmax_m])
        % ylim([-globalmax_m globalmax_m])
        % zlim([-globalmax_m globalmax_m])
        title(axsat,sprintf('Roll = %.2f^o, Pitch =%.2f^o, Yaw = %.2f^o',roll(i),pitch(i),yaw(i)))
        xlabel(axsat,'X Axis (m)')
        ylabel(axsat,'Y Axis (m)')
        zlabel(axsat,'Z Axis (m)')
        axis(axsat,'equal')
        FontWidthandPos
        %Make over lay plot
        set(h,'currentaxes',ax2);
        plot(ax2,[timestoit]./60,roll)
        hold on
        plot(ax2,[timestoit]./60,pitch)
        plot(ax2,[timestoit]./60,yaw)
        plot(ax2,[timestoit(i) timestoit(i)]./60,[-300 300],'k')
        hold off
        grid on
        ylim(ax2,[-200 200])
        ylabel(ax2,'Degrees')
        xlabel(ax2,'Minutes')
        legend(ax2,{'Roll','Pitch','Yaw'},'Location','northoutside')
        FontWidthandPos
        drawnow
        M(i) = getframe(gcf);
        if i == 1
    %         h.Visible = 'off';
        end
        pause(.001)
    end
    if 0
        %Send movie to gif
        filename = 'Rotation_GIF.gif';
        % Capture the plot as an image 
              frame = getframe(h); 
        for n = 1:length(M)
            frame = M(n);
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            % Write to the GIF File 
            if n == 1 
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
            else 
                imwrite(imind,cm,filename,'gif','WriteMode','append'); 
            end
        end
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