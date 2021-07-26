function [lat, lon, utcvec] = ECI2latlon(R_IJK,t_sec,init_utcvec)
% [lat, lon, utcvec] = ECI2latlon(R_IJK,t_sec,init_utcvec)
%
% This function is used to convert from the Earth Centered Inertial frame
% to latitude and longitude
%
% Inputs: R_IJK - Position Vector in Earth Centered Inertial (km)
%            t_sec - Seconds since start of init_utcvec (Sec)
%            init_utcvec - Vector of [year, month, day, hour, minute, sec]
%
% Outputs: lat - Latitude (deg)
%                lon - Longitude (deg)
%
% Created by Bailey Miller 3/31/2021
% See also: ECEF2latlon.m, dcmeci2ecef.m

%Check for inputs
if ~exist('t_sec','var') || isempty(t_sec)
    warning('t_sec set to 0')
    t_sec = 0;
end
default_utcvec = [2020 1 1 0 0 0];
if ~exist('init_utcvec','var')
    init_utcvec = default_utcvec;
elseif isempty(init_utcvec) || length(init_utcvec)<6
    warning('init_utcvec is the wrong length. Setting to default 2020/1/1')
    init_utcvec = default_utcvec;
end

%Make UTC vector based on t_sec since start
    sec_now = init_utcvec(6)+ t_sec; 
    min_now = init_utcvec(5); hr_now = init_utcvec(4); day_now = init_utcvec(3);
    mnth_now = init_utcvec(2); yr_now = init_utcvec(1);
    %Modulate based on seconds
    while sec_now >= 60
        %Subtract 60 seconds
        sec_now = sec_now - 60;
        %Add one minute
        min_now = min_now +1;
        %Check minutes
        if min_now >= 60
            %Subtract 60 minutes
            min_now = min_now - 60;
            %Add 1 hr
            hr_now = hr_now+1;
        end
        %Check hours
        if hr_now >=24
            %Subtract 24 hours
            hr_now = hr_now - 24;
            %Add 1 day
            day_now = day_now +1;
        end
        %Check days
        if any(mnth_now==[1 3 5 7 8 10 12]) && day_now > 31
            day_now = 1;
            %Add one month
            mnth_now = mnth_now+1;
        elseif mnth_now == 2
            if rem(yr_now,4) == 0 && day_now >29 %leap year
                day_now = 1;
                %Add one month
                mnth_now = mnth_now+1;
            elseif day_now >= 28
                day_now = 1;
                mnth_now = mnth_now+1;
            end
        elseif any(mnth_now== [4 6 9 11]) && day_now > 30
            day_now = 1;
            mnth_now = mnth_now +1;
        end
        %Check months
        if mnth_now >12
            mnth_now = 1;
            %Add one year
            yr_now = yr_now +1;
        end
    end
    
    %Make utc vec and determine dcm
    utcvec = [yr_now, mnth_now, day_now, hr_now, min_now, sec_now];
    dcm = dcmeci2ecef('IAU-2000/2006',utcvec);
    %Do coordinate transformation from ECI to ECEF then solve
    R_xyz = R_IJK.'*dcm;
    [lat,lon] = ECEF2latlon(R_xyz);
    
    %     %Change longitude based on earth's rotation
%     P_s = 24*3600; %seconds length of mean solar day
%     P_e = 86164; %seconds rotation period of earth
%     lon_mo = lon-360*(1/P_e-1/P_s)*state(i).t;
%     if lon_mo>=180
%         lon_mo = lon_mo-360*floor(i/n*N);
%     end
%     lon = lon_mo;
end