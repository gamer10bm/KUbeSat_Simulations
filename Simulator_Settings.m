%This script can programmatically set the constant values throughout the
%simulink model
if ~exist('model_name','var')
    model_name = 'SatStates';
end
%Set constant values as "" strings
    %Earth Constants
    const_struc = struct('name','mu_km3s2','val',"3.986e5");
    const_struc(end+1) = struct('name','J2','val',"1.082626638e-3");
    const_struc(end+1) = struct('name','eps_earth_deg','val',"23.5");
    const_struc(end+1) = struct('name','Re_km','val',"6378");
    const_struc(end+1) = struct('name','flat_rat','val',"0.08182");
    const_struc(end+1) = struct('name','p0','val',"4.58e-6");
    
    %Orbit Constants
    const_struc(end+1) = struct('name','a_km','val',"555");
    const_struc(end+1) = struct('name','e','val',"0.001");
    const_struc(end+1) = struct('name','i_deg','val',"97.6123");
    const_struc(end+1) = struct('name','RAAN_deg','val',"90");
    const_struc(end+1) = struct('name','argPeri_deg','val',"0");
    const_struc(end+1) = struct('name','TA_deg','val',"0");
    const_struc(end+1) = struct('name','year','val',"2022");
    const_struc(end+1) = struct('name','month','val',"1");
    const_struc(end+1) = struct('name','day','val',"1");
    
    %Location Constants
    const_struc(end+1) = struct('name','Lat_GS_deg','val',"39.96");
    const_struc(end+1) = struct('name','Long_GS_deg','val',"-95.25");
    const_struc(end+1) = struct('name','alt_GS_km','val',"3.6");
    const_struc(end+1) = struct('name','Lat_HiCal_deg','val',"-60");
    
    %Satellite Constants
    const_struc(end+1) = struct('name','Sol_Axyz_m2','val',"[0.03 0.01 0]");
    const_struc(end+1) = struct('name','S_W_m2','val',"1370");
    const_struc(end+1) = struct('name','nSolar','val',"0.28"); %0.185
    const_struc(end+1) = struct('name','nPath','val',"0.9");  
    const_struc(end+1) = struct('name','batt_init_Whr','val',"40");
    const_struc(end+1) = struct('name','pow_PCRD_W','val',"3");
    const_struc(end+1) = struct('name','data_PCRD_kbps','val',"20");
    const_struc(end+1) = struct('name','pow_HiCal_W','val',"0.5");
    const_struc(end+1) = struct('name','data_HiCal_kbps','val',"0.1");
    const_struc(end+1) = struct('name','pow_CDH_W','val',"0.5");
    const_struc(end+1) = struct('name','pow_Comms_W','val',"0.3");
    const_struc(end+1) = struct('name','data_Comms_kbps','val',"250");
    const_struc(end+1) = struct('name','pow_EPS_W','val',"0.2");
    const_struc(end+1) = struct('name','pow_BattHeat_W','val',"0.8");
    const_struc(end+1) = struct('name','Cd','val',"1.28");
    const_struc(end+1) = struct('name','COM','val',"[1.14122; 19.92782; 0.54375]*0.001");
    const_struc(end+1) = struct('name','Area','val',"[0.0341; 0.01; 0.03141]");
    const_struc(end+1) = struct('name','J_Inert','val',"diag([1082793.489*20; 4354.008742*223; 1085536.436*20])/1000^3");
    const_struc(end+1) = struct('name','H0','val',"[.3e-3;0.01e-3;.3e-3]");
    


%Load model
open_system(model_name)
%Edit model parameters
init_name = sprintf('%s/%s',model_name,'Initialize');
for c_id = 1:length(const_struc)
    %Loop through const_struc
    set_param(sprintf('%s/%s',init_name,const_struc(c_id).name),'Value',const_struc(c_id).val);
end
%Save model
save_system(model_name);