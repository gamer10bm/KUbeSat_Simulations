%% Mess with orientation sequence

%Fill a vector with rotvec angles and sequence of seqvec axes
angsweep_rad = {0};
seqvec = {'x'};
timevec = [0];
dt = 15; %seconds per step
sweepn = 180; 
finaltime = 360; %sec
rpy_vec_deg= [0 0 0];

%Add angle commands
yawrate_deg = 2; %degrees per second
yawgoal1_deg = 720; %degrees in yaw for next maneuver
yawgoal2_deg=360;
pitchrate_deg = 1; %degrees per second
pitchgoal1_deg = 90;
pitchgoal2_deg = -90;
n = 1;
while n<sweepn && yawrate_deg*timevec(n) < abs(yawgoal1_deg)
    n = n+1;
    timevec(n) = timevec(n-1)+dt;
    angstep = [];
    seqstep = '';
    %Do yaw command
    if yawrate_deg*timevec(n) <= abs(yawgoal1_deg)
        angstep = [angstep sign(yawgoal1_deg)*deg2rad(yawrate_deg)*dt];
        seqstep = [seqstep 'z'];
    end
    
    %Do pitch command
%     if pitchrate_deg*timevec(n) < pitchgoal_deg
%         angstep = [angstep deg2rad(pitchrate_deg)*dt];
%         seqstep = [seqstep 'y'];
%     end
    
    %Load command vectors
    angsweep_rad{n} = angstep;
    seqvec{n} = seqstep;
end

% 
timestrt = timevec(n);
while pitchrate_deg*(timevec(n)-timestrt) < abs(pitchgoal1_deg) && n<sweepn
    n = n+1;
    timevec(n) = timevec(n-1)+dt;
    angsweep_rad{n} = sign(pitchgoal1_deg)*deg2rad(pitchrate_deg)*dt;
    seqvec{n} = ['y'];
end

timestrt = timevec(n);
while n<sweepn && yawrate_deg*(timevec(n)-timestrt) < abs(yawgoal2_deg)
    n = n+1;
    timevec(n) = timevec(n-1)+dt;
    angstep = [];
    seqstep = '';
    
    %Do yaw command
    angstep = [angstep sign(yawgoal2_deg)*deg2rad(yawrate_deg)*dt];
    seqstep = [seqstep 'z'];
    
    %Load command vectors
    angsweep_rad{n} = angstep;
    seqvec{n} = seqstep;
end

timestrt = timevec(n);
while pitchrate_deg*(timevec(n)-timestrt) <  abs(pitchgoal2_deg) && n<sweepn
    n = n+1;
    timevec(n) = timevec(n-1)+dt;
    angsweep_rad{n} = sign(pitchgoal2_deg)*deg2rad(pitchrate_deg)*dt;
    seqvec{n} = ['y'];
end
sweepn = n;

%Initialize for quaternion calculations
q = cell(1,sweepn);
roll = zeros(1,sweepn);
pitch = zeros(1,sweepn);
yaw = zeros(1,sweepn);
Rotmat = cell(1,sweepn);
rpy_vec_rad = deg2rad(rpy_vec_deg); %Initial orientation
for i = 1:sweepn
    q{i}= rotseq2quaternion(angsweep_rad{i},seqvec{i},rpy_vec_rad); %angsweep_rad(i)
    [roll(i), pitch(i), yaw(i),Rotmat{i}] = quaternion2euler(q{i});
    rpy_vec_rad= deg2rad([roll(i) pitch(i) yaw(i)]);
end
% keyboard

%Make a moooovie
h = figure(4); clf
axsat = axes('NextPlot','replacechildren','OuterPosition',[0 0 .8 1]);
ax2 = axes('Box','on','OuterPosition',[0.6 .4 .4 .6],'NextPlot','replacechildren');
clearvars M
M(sweepn) = struct('cdata',[],'colormap',[]);
for i = 1:sweepn
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
    plot(ax2,timevec./60,roll)
    hold on
    plot(ax2,timevec./60,pitch)
    plot(ax2,timevec./60,yaw)
    plot(ax2,[timevec(i) timevec(i)]./60,[-300 300],'k')
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

%Send movie to gif
filename = 'HorizontalSpin_Reorient_VerticalSpin_GIF.gif';
if exist(filename,'file')
    delete(filename)
end
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