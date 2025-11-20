clear 
addpath /discover/nobackup/qliu/gdelanno_RTM/MATLAB_LDASSA

data_path = '/discover/nobackup/projects/land_da/Evaluation/IVs/output/';

mod_version = 'OLv8_M36_cd';

mod_resolution = 'M36';

% number of days in lag
Nlag = 2;

Nmin = 100;

start_time.year = 2018;
start_time.month = 8; 
start_time.day = 1;
start_time.hour = 0;
start_time.min = 0;
start_time.sec = 0;

end_time.year = 2024; 
end_time.month = 6;  
end_time.day = 30;
end_time.hour = 0;
end_time.min = 0;
end_time.sec = 0;

if end_time.month == 1
    time_tag = [num2str(start_time.year),num2str(start_time.month,'%2.2d'),'_',...
        num2str(end_time.year-1),'12'];
else
    time_tag = [num2str(start_time.year),num2str(start_time.month,'%2.2d'),'_',...
        num2str(end_time.year),num2str(end_time.month-1,'%2.2d')];
end

% get EASEv2 coord
if strcmp(mod_resolution,'M36')
    [lat,lon] = smapeasev2_ind2latlon([0:405],[0:963],'M36');
    lon_L4 = repmat(lon',[1,length(lat)]);
    lat_L4 = repmat(lat,[length(lon),1]);
    clear lat lon
else
    f_EASE = '/home/qliu/smap/SMAP_L4/L4_SM/Vv5014/lmc/SMAP_L4_SM_lmc_00000000T000000_Vv5014_001.h5';
    lon_L4 = double(h5read(f_EASE,'/cell_lon'));
    lat_L4 = double(h5read(f_EASE,'/cell_lat'));
end

Nlon = size(lon_L4,1);
Nlat = size(lon_L4,2);

mod_sm_sum =  zeros(Nlon*Nlat,1);
mod_sm2_sum =  zeros(Nlon*Nlat,1);

obs_sm_sum =  zeros(Nlon*Nlat,1);
obs_sm2_sum =  zeros(Nlon*Nlat,1);

modxobs_sm_sum =  zeros(Nlon*Nlat,1);
modxlag1_sm_sum =  zeros(Nlon*Nlat,1);
obsxlag1_sm_sum =  zeros(Nlon*Nlat,1);
obslag1xmod_sm_sum = zeros(Nlon*Nlat,1);
modlag1xobs_sm_sum = zeros(Nlon*Nlat,1);

N_sm =  zeros(Nlon*Nlat,1);

% load climatology file

f_clim= [data_path,'ASCL4_',mod_version,'_clim_pentad_',time_tag,'_w31.mat'];
tmp = load(f_clim,'mod_sm_clim','obs_sm_clim');
mod_clim = tmp.mod_sm_clim;
obs_clim = tmp.obs_sm_clim;
clear tmp

start_time = augment_date_time(Nlag*86400,start_time); % start from day 2 for computing lag 1 autocov
date_time = start_time;
date_time = get_dofyr_pentad(date_time);
date_time_pre = augment_date_time(-Nlag*86400,date_time);
date_string_pre = [num2str(date_time_pre.year,'%4.4d'), ...
        num2str(date_time_pre.month,'%2.2d'),num2str(date_time_pre.day,'%2.2d')];
    
while 1
    
    if (date_time.year ==end_time.year   && ...
            date_time.month==end_time.month  && ...
            date_time.day  ==end_time.day    )
        break
        
    end
    
    date_string = [num2str(date_time.year,'%4.4d'), ...
        num2str(date_time.month,'%2.2d'),num2str(date_time.day,'%2.2d')];
    
    date_time_pre = augment_date_time(-Nlag*86400,date_time);
    date_string_pre = [num2str(date_time_pre.year,'%4.4d'), ...
        num2str(date_time_pre.month,'%2.2d'),num2str(date_time_pre.day,'%2.2d')];
    
    fname = [data_path,'ASCL4_H119_SMSF_L4_',mod_version,'_QC_1_',date_string,'.mat'];
    fname_pre = [data_path,'ASCL4_H119_SMSF_L4_',mod_version,'_QC_1_',date_string_pre,'.mat'];
    
    if exist(fname,'file') && exist(fname_pre,'file')
        disp(fname)
        tmp = load(fname,'sm_mod','sm_obs','idx_EASEv2_lonxlat');
        tmp_pre = load(fname_pre,'sm_mod','sm_obs','idx_EASEv2_lonxlat'); tmp_pre.sm_obs = tmp_pre.sm_obs;
        
        [idx,Inow,Ipre] = intersect(tmp.idx_EASEv2_lonxlat,tmp_pre.idx_EASEv2_lonxlat);
        
        if ~isempty(idx)
            % compute anomaly
            tmp_clim = mod_clim(:, date_time.pentad);
            sm_mod = tmp.sm_mod(Inow) - tmp_clim( idx); clear tmp_clim
            
            tmp_clim = obs_clim(:, date_time.pentad);
            sm_obs = tmp.sm_obs(Inow) - tmp_clim(idx); clear tmp_clim
            
            tmp_clim = mod_clim(:, date_time_pre.pentad);
            sm_mod_pre = tmp_pre.sm_mod(Ipre) - tmp_clim( idx); clear tmp_clim
            
            tmp_clim = obs_clim(:, date_time_pre.pentad);
            sm_obs_pre = tmp_pre.sm_obs(Ipre) - tmp_clim(idx); clear tmp_clim
            
            
            iv = find(~isnan(sm_mod));
            % convert compressed data back on global
            mod_sm_sum(idx(iv)) = mod_sm_sum(idx(iv)) + sm_mod(iv);
            mod_sm2_sum(idx(iv)) = mod_sm2_sum(idx(iv)) + sm_mod(iv).^2;
            
            obs_sm_sum(idx(iv)) = obs_sm_sum(idx(iv)) + sm_obs(iv);
            obs_sm2_sum(idx(iv)) = obs_sm2_sum(idx(iv)) + sm_obs(iv).^2;
            
            modxobs_sm_sum(idx(iv)) = modxobs_sm_sum(idx(iv)) + sm_mod(iv) .* sm_obs(iv);
            modxlag1_sm_sum(idx(iv)) = modxlag1_sm_sum(idx(iv)) + sm_mod(iv) .* sm_mod_pre(iv);
            obsxlag1_sm_sum(idx(iv)) = obsxlag1_sm_sum(idx(iv)) + sm_obs(iv) .* sm_obs_pre(iv);
            obslag1xmod_sm_sum(idx(iv)) = obslag1xmod_sm_sum(idx(iv)) + sm_mod(iv) .* sm_obs_pre(iv);
            modlag1xobs_sm_sum(idx(iv)) = modlag1xobs_sm_sum(idx(iv)) + sm_mod_pre(iv) .* sm_obs(iv);
            
            N_sm(idx(iv)) = N_sm(idx(iv)) + 1;
            
            clear sm_mod sm_obs tmp iv sm_mod_pre sm_obs_pre
            
        end
        
        clear idx fname fname_pre
    end
    
    date_time = augment_date_time(86400, date_time);
    
end

NN = N_sm;
NN(NN < Nmin) = NaN;

mod_sm_mean = mod_sm_sum ./NN;
obs_sm_mean = obs_sm_sum ./NN;

mod_sm_mean(isinf(mod_sm_mean)) = NaN;
obs_sm_mean(isinf(obs_sm_mean)) = NaN;

C_mod_mod = mod_sm2_sum ./NN - mod_sm_mean .^2;
C_mod_mod(C_mod_mod < 0.) = NaN;
C_obs_obs = obs_sm2_sum ./NN - obs_sm_mean .^2;
C_obs_obs(C_obs_obs < 0.) = NaN;

C_mod_obs = modxobs_sm_sum ./NN  - mod_sm_mean .* obs_sm_mean; 
%C_mod_obs(C_mod_obs < 0.) = NaN;

R_mod_obs = C_mod_obs ./ sqrt(C_mod_mod .*C_obs_obs);
R_mod_obs(isnan(R_mod_obs)) = -9999.;
 
C_mod_obs(C_mod_obs < 0.) = NaN;

C_mod_modlag1 = modxlag1_sm_sum ./NN - mod_sm_mean .^2;  C_mod_modlag1(C_mod_modlag1 <0) = NaN;

C_obs_obslag1 = obsxlag1_sm_sum ./NN - obs_sm_mean .^2; C_obs_obslag1(C_obs_obslag1 <0) = NaN;

C_mod_obslag1 = obslag1xmod_sm_sum./NN - mod_sm_mean .* obs_sm_mean;

C_modlag1_obs = modlag1xobs_sm_sum./NN - mod_sm_mean .* obs_sm_mean;

S_ivd = sqrt(C_mod_modlag1./C_obs_obslag1);

S_ivs_obs = C_mod_obslag1 ./C_obs_obslag1 ;
S_ivs_mod = C_mod_modlag1 ./C_modlag1_obs ;

R2_ivd_mod = C_mod_obs .* S_ivd ./ C_mod_mod; R2_ivd_mod(R2_ivd_mod <0 ) = NaN;
R2_ivd_mod(R2_ivd_mod > 1.) = 1.;
R2_ivd_obs = C_mod_obs ./C_obs_obs ./ S_ivd; R2_ivd_obs(R2_ivd_obs <0 ) = NaN;
R2_ivd_obs(R2_ivd_obs > 1.) = 1.;

R2_ivs_mod = C_mod_obs .* S_ivs_obs ./ C_mod_mod; R2_ivs_mod(R2_ivs_mod <0 ) = NaN;
R2_ivs_mod(R2_ivs_mod > 1.) = 1.;
R2_ivs_obs = C_mod_obs ./C_obs_obs ./ S_ivs_obs; R2_ivs_obs(R2_ivs_obs <0 ) = NaN;
R2_ivs_obs(R2_ivs_obs > 1.) = 1.;

fname_out = [data_path,'ASCL4_',mod_version,'_IVD_IVS_stats_lag',num2str(Nlag),'day_',time_tag,'.mat'];
save(fname_out,'N_sm','Nmin','Nlag','R2_ivd_mod','R2_ivd_obs','R2_ivs_mod','R2_ivs_obs','R_mod_obs')
