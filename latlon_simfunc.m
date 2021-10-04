function [lat, lon] = fcn(R,time,year,month,day)

yr_init = year; mnth_init = month; day_init = day;
hr_init = 0; min_init = 0; sec_init = 0;

init_utcvec = [yr_init, mnth_init, day_init, hr_init, min_init, sec_init];

[lat, lon] = ECI2latlon(R',time,init_utcvec);
end