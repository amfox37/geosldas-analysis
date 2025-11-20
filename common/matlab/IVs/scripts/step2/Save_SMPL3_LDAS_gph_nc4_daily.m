clear
addpath /discover/nobackup/qliu/gdelanno_RTM/MATLAB_LDASSA/
addpath /home/qliu/projects/matlab_code/SMAP/

obs_path = '/home/qliu/smap/SMAP_Nature/SMAP/SPL3SMP.008/';

%mod_path = '/discover/nobackup/amfox/Experiments/DAv7_M36_ASCAT_type_2/';
%mod_run = 'DAv7_M36_ASCAT_type_2';
%mod_version = 'DAv7ASCt2_M36';

mod_path = '/discover/nobackup/amfox/Experiments/DAv7_M36_ASCAT_type_13_comb_fp_scaled/';
mod_run = 'DAv7_M36_ASCAT_type_13_comb_fp_scaled';
mod_version = 'DAv7_M36_ASCAT_type_13_comb_fp_scaled';

out_collection = '.SMAP_L4_SM_gph.';
domain = 'SMAP_EASEv2_M36_GLOBAL';

qc_yes = 1;

out_path = '/discover/nobackup/amfox/IVs_output/';

orbit = {'AM','PM'};

start_time.year = 2015;
start_time.month = 4; 
start_time.day = 1;
start_time.hour = 0;
start_time.min = 0;
start_time.sec = 0;

end_time.year = 2021; 
end_time.month = 4; 
end_time.day = 1;
end_time.hour = 0;
end_time.min = 0;
end_time.sec = 0;

[lat,lon] = smapeasev2_ind2latlon([0:405],[0:963],'M36');
lon_L4 = repmat(lon',[1,length(lat)]);
lat_L4 = repmat(lat,[length(lon),1]);
clear lat lon

f_tilecoord = [mod_path,mod_run,'/output/',domain,'/rc_out/',mod_run,'.ldas_tilecoord.bin'];
tile_coord = read_tilecoord(f_tilecoord);
tc = tile_coord;

date_time = start_time;

while 1
    
    if (date_time.year ==end_time.year   && ...
            date_time.month==end_time.month  && ...
            date_time.day  ==end_time.day    )
        break
        
    end
    
    date_string = [num2str(date_time.year,'%4.4d'), ...
        num2str(date_time.month,'%2.2d'),num2str(date_time.day,'%2.2d')];
    
    flist = dir([obs_path, '/Y',num2str(date_time.year,'%4.4d'), ...
        '/SMAP_L3_SM_P_',date_string,'_R18290_*.h5']);
    
    if length(flist) >= 1
        f_obs = [obs_path, '/Y',num2str(date_time.year,'%4.4d'), ...
            '/',flist(end).name];
        clear flist
        
        clear sm_tmp
        for iorb = 1:length(orbit)
            
            grp_tag = orbit{iorb};
            
            if contains(orbit{iorb},'PM')
                var_tag =  '_dca_pm';
                var_tag1 = '_pm';
            else
                var_tag ='';
                var_tag1 ='';
            end
            
            L3_sm = h5read(f_obs,['/Soil_Moisture_Retrieval_Data_',grp_tag,'/soil_moisture',var_tag]);
            L3_sm(L3_sm < 0 ) = NaN;
            % quality flag
            L3_qf = h5read(f_obs,['/Soil_Moisture_Retrieval_Data_',grp_tag,'/retrieval_qual_flag',var_tag]);
            % surface status land = 0, nonland = 1
            L3_ss = h5read(f_obs,['/Soil_Moisture_Retrieval_Data_',grp_tag,'/grid_surface_status',var_tag1]);
            % surface flag
            L3_sf = h5read(f_obs,['/Soil_Moisture_Retrieval_Data_',grp_tag,'/surface_flag',var_tag1]);
            L3_rfi_h = h5read(f_obs,['/Soil_Moisture_Retrieval_Data_',grp_tag,'/tb_qual_flag_h',var_tag1]);
            L3_rfi_v = h5read(f_obs,['/Soil_Moisture_Retrieval_Data_',grp_tag,'/tb_qual_flag_v',var_tag1]);
            
            if qc_yes
                % QC based on retrieval quality flag
                %L3_rt = bitget(L3_qf, 1); % only use retrieval_recommended
                L3_rt = bitget(L3_qf, 3); % use retrieval_succeeded
                
                % QC b ased on surface flag
                L3_frozen_model = bitget(L3_sf, 9); % model frozen ground
                L3_snow = bitget(L3_sf,6); % snow and ice
                L3_pice = bitget(L3_sf,7); % permanent snow and ice
                L3_rfi_h_qf = bitget(L3_rfi_h,1); % quality flag RFI H
                L3_rfi_v_qf = bitget(L3_rfi_v,1); % quality flag RFI V
                idx_bad = find(L3_rt ~= 0 | L3_frozen_model ~=0 | ...
                    L3_snow ~= 0 | L3_pice ~=0 | L3_ss ~= 0 | ...
                    L3_rfi_h_qf ~= 0 | L3_rfi_v_qf ~= 0);
                
                % only keep data/coord that pass QC
                L3_sm(idx_bad) = NaN;
                clear idx_bad
                
            end
            sm_tmp(:,:,iorb) = L3_sm; clear L3_sm
        end
        sm_obs = nanmean(sm_tmp,3);
        
        clear sm_tmp
        
    else
        sm_obs = NaN *ones(size(lon_L4));
    end
   
    % check aggregated daily files first 
    fname = [mod_path, mod_run,'/output/',domain,'/cat/ens_avg/Y', ...
        num2str(date_time.year, '%4.4d'), ...
        '/M', num2str(date_time.month, '%2.2d'), ...
        '/',mod_run, out_collection, ...
        num2str(date_time.year, '%4.4d'), ...
        num2str(date_time.month, '%2.2d'), ...
        num2str(date_time.day,'%2.2d'), '.nc4'];
    
    if exist(fname,'file')
        disp(fname)
        smsf_tmp = ncread(fname,'sm_surface');
        smsf_tmp(smsf_tmp < 0 ) = NaN;
        sm_mod_3hr = smsf_tmp; clear smsf_tmp
        
    else
        
        hr_list = [1:3:22];
        
        for this_hr = 1:8
            date_time_string =[num2str(date_time.year, '%4.4d'), ...
                num2str(date_time.month, '%2.2d'), ...
                num2str(date_time.day, '%2.2d'), ...
                '_',num2str(hr_list(this_hr), '%2.2d'), ...
                num2str(30, '%2.2d')];
            
            fname = [mod_path, mod_run,'/output/',domain,'/cat/ens_avg/Y', ...
                num2str(date_time.year, '%4.4d'), ...
                '/M', num2str(date_time.month, '%2.2d'), ...
                '/',mod_run, out_collection, ...
                date_time_string,'z.nc4'];
            
            if this_hr == 1, disp(fname), end
            
            if exist(fname,'file')
                
                smsf_tmp = ncread(fname,'sm_surface');
                smsf_tmp(smsf_tmp < 0 ) = NaN;
                sm_mod_3hr(:,this_hr) = smsf_tmp; clear smsf_tmp
                
            else
                
                disp(['File doesn''t exist, stop ',fname])
                return
            end
        end
    end % end hourly file loop
   
    % compute daily mean 
    sm_mod = nanmean(sm_mod_3hr,2); clear sm_mod_3hr
    if(length(sm_mod) ~= tc.N_tile)
        error('something is wrong reading LDAS data')
    end
    
    sm_mod_glob_2d = NaN * sm_obs;
    
    for k = 1:length(sm_mod)
        sm_mod_glob_2d(tc.i_indg(k)+1,tc.j_indg(k)+1) = sm_mod(k);
    end
    
    clear sm_mod
    sm_mod = sm_mod_glob_2d; clear sm_mod_glob_2d
    
    sm_obs(isnan(sm_mod)) = NaN;
    sm_mod(isnan(sm_obs)) = NaN;
    
    idx_EASEv2_lonxlat = find(~isnan(sm_obs(:)));
    
    if ~isempty(idx_EASEv2_lonxlat)
        
        sm_obs = sm_obs(idx_EASEv2_lonxlat);
        sm_mod = sm_mod(idx_EASEv2_lonxlat);
        
    else
        sm_obs = [];
        sm_mod = [];
    end
    
    f_out = [out_path,'SMPL3_R18_SMSF_MOD_',mod_version,'_QC_1_',date_string,'.mat'];
    
    save(f_out,'sm_mod','sm_obs','idx_EASEv2_lonxlat');
    
    clear sm_obs sm_mod idx_EASEv2_lonxlat
    
    date_time = augment_date_time(86400, date_time);
end


