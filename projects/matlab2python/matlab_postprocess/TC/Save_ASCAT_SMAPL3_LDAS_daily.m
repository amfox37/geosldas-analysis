clear
addpath /home/qliu/projects/matlab_code/SMAP/

ASCAT_path = '/discover/nobackup/qliu/merra_land/DATA/ASCAT_HSAF/';
ASCAT_version = 'H119_H120';

SMAPL3_path = '/home/qliu/smap/SMAP_Nature/SMAP/SPL3SMP.008/';
SMAPL3_version = 'R18290';
orbit = {'AM','PM'};
qc_yes = 1;

%MOD_path = '/home/qliu/smap/SMAP_Nature/SMAP_Nature_v11/';
%MOD_run = 'SMAP_Nature_v11.4_M36';
%MOD_version = 'NRv11.4_M36';

MOD_path = '/discover/nobackup/borescan/par/';
MOD_run = 'SMAP_Nature_v11.4_M36_exp';
MOD_version = 'NRv11.4exp_M36';
%MOD_path = '/discover/nobackup/amfox/Experiments/OLv7_M36_ascat/';
%MOD_run = 'OLv7_M36_ascat';
%MOD_version = 'OLv7_M36';

out_collection = '.SMAP_L4_SM_gph.';
domain = 'SMAP_EASEv2_M36_GLOBAL';

resolution = 'M36';

out_path = '/discover/nobackup/qliu/merra_land/matlab/SMAP/L3L4ASCAT/';

start_time.year = 2015;
start_time.month = 4;
start_time.day = 1;
start_time.hour = 0;
start_time.min = 0;
start_time.sec = 0;

end_time.year = 2023; %start_time.year+1;
end_time.month = 4;
end_time.day = 1;
end_time.hour = 0;
end_time.min = 0;
end_time.sec = 0;

% get ASCAT GPI coord
f_gpi = [ASCAT_path,'Auxiliary/TUW_WARP5_grid_info_2_2.nc'];

land_flag = ncread(f_gpi, 'land_flag');
gpi = ncread(f_gpi,'gpi');
gpi_land = gpi(land_flag == 1);
lon_gpi = ncread(f_gpi,'lon');
lon_gpi = double(lon_gpi(land_flag ==1));
lat_gpi = ncread(f_gpi,'lat');
lat_gpi = double(lat_gpi(land_flag ==1));
clear gpi

% get EASEv2 coord
if strcmp(resolution,'M09')
    f_EASE = '/css/smapl4/public/L4_Products/L4_SM/Vv7030/lmc/SMAP_L4_SM_lmc_00000000T000000_Vv4030_001.h5';
    lon_L4 = double(h5read(f_EASE,'/cell_lon'));
    lat_L4 = double(h5read(f_EASE,'/cell_lat'));
elseif strcmp(resolution,'M36')
    [lat,lon] = smapeasev2_ind2latlon([0:405],[0:963],'M36');
    lon_L4 = repmat(lon',[1,length(lat)]);
    lat_L4 = repmat(lat,[length(lon),1]);
    clear lat lon
else
    error('incorrect resolution')
end

f_tilecoord = [MOD_path,MOD_run,'/output/',domain,'/rc_out/',MOD_run,'.ldas_tilecoord.bin'];
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
    
    % read daily ASCAT data
    f_ASCAT = [ASCAT_path, ASCAT_version,'_processed/Y',num2str(date_time.year,'%4.4d'), ...
        '/M',num2str(date_time.month,'%2.2d'),'/' ...
        'ASCAT_HSAF_H119_SM_',date_string,'_AD.mat']
    
    if exist(f_ASCAT,'file')
        ASC = load(f_ASCAT);
        
        sm_asc = double(ASC.sm_tile);
        sm_asc(sm_asc > 100) = NaN;
 
        % H115 more aggressive QC: only keep data with conf_flag == 0
        % i.e. remove data with larege noise, wetland frac, top complexity
        % and sensitivity level (see PUM for flag details)
        % the new QC results in ~10% fewer data compared to previous QC
        conf_flag = ASC.conf_flag_tile;
        sm_asc(conf_flag >= 1) = NaN;
        
        clear ASC
        
        sm_ASC  = griddata(lon_gpi,lat_gpi,sm_asc,lon_L4,lat_L4,'natural'); clear sm_asc
    else
        sm_ASC = NaN *ones(size(lon_L4));
    end
    
    % read SMAP L3 data
    flist = dir([SMAPL3_path, '/Y',num2str(date_time.year,'%4.4d'), ...
        '/SMAP_L3_SM_P_',date_string,'_',SMAPL3_version,'_*.h5']);
    
    if length(flist) >= 1
        f_obs = [SMAPL3_path, '/Y',num2str(date_time.year,'%4.4d'), ...
            '/',flist(end).name];
        disp(f_obs)
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
        sm_L3 = nanmean(sm_tmp,3);
        
        clear sm_tmp
        
    else
        sm_L3 = NaN *ones(size(lon_L4));
    end
    
    %read model data
    
    fname = [MOD_path, MOD_run,'/output/',domain,'/cat/ens0000/Y', ...
        num2str(date_time.year, '%4.4d'), ...
        '/M', num2str(date_time.month, '%2.2d'), ...
        '/',MOD_run, out_collection, ...
        num2str(date_time.year, '%4.4d'), ...
        num2str(date_time.month, '%2.2d'), ...
        num2str(date_time.day,'%2.2d'), '_1200z.nc4'];
    
    if ~exist(fname,'file')
        fname = strrep(fname,'/ens0000/','/ens_avg/');
    end
    
    if ~exist(fname,'file')
        fname = [MOD_path, MOD_run,'/output/',domain,'/cat/ens0000/Y', ...
            num2str(date_time.year, '%4.4d'), ...
            '/M', num2str(date_time.month, '%2.2d'), ...
            '/',MOD_run, out_collection, ...
            num2str(date_time.year, '%4.4d'), ...
            num2str(date_time.month, '%2.2d'), ...
            num2str(date_time.day,'%2.2d'), '.nc4'];
    end
    
    if ~exist(fname,'file')
        fname = strrep(fname,'/ens0000/','/ens_avg/');
    end
    
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
            
            fname = [MOD_path, MOD_run,'/output/',domain,'/cat/ens0000/Y', ...
                num2str(date_time.year, '%4.4d'), ...
                '/M', num2str(date_time.month, '%2.2d'), ...
                '/',MOD_run, out_collection, ...
                date_time_string,'z.nc4'];
            
            if ~exist(fname,'file')
                fname = strrep(fname,'/ens0000/','/ens_avg/');
            end
            
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
    
    sm_mod = nanmean(sm_mod_3hr,2); clear sm_mod_3hr
    
    if(length(sm_mod) ~= tc.N_tile)
        error('something is wrong reading LDAS data')
    end
    
    sm_mod_glob_2d = NaN * sm_L3;
    
    for k = 1:length(sm_mod)
        sm_mod_glob_2d(tc.i_indg(k)+1,tc.j_indg(k)+1) = sm_mod(k);
    end
    clear sm_mod
    
    sm_mod = sm_mod_glob_2d; clear sm_mod_glob_2d
    
    sm_ASC(isnan(sm_L3)) = NaN;
    sm_L3(isnan(sm_ASC)) = NaN;
    sm_mod(isnan(sm_ASC)) = NaN;
    
    idx_lonxlat = find(~isnan(sm_ASC(:)));
    
    if ~isempty(idx_lonxlat)
        
        sm_ASC = sm_ASC(idx_lonxlat);
        sm_L3 = sm_L3(idx_lonxlat);
        sm_mod = sm_mod(idx_lonxlat);
        
    else
        sm_ASC = [];
        sm_L3 = [];
        sm_mod = [];
    end
    
    f_out = [out_path,'ASCAT_',ASCAT_version(1:4), ...
        '_L3_',SMAPL3_version,'_',MOD_version,'_SMSF_',date_string,'.mat'];
    
    save(f_out,'sm_mod','sm_L3','sm_ASC','idx_lonxlat');

    clear sm_ASC sm_L3 sm_mod idx_lonxlat
    
    date_time = augment_date_time(86400, date_time);
    
end




