clear 
addpath /discover/nobackup/qliu/gdelanno_RTM/MATLAB_LDASSA/
addpath /home/qliu/projects/matlab_code/SMAP/

ASCAT_path = '/discover/nobackup/qliu/merra_land/DATA/ASCAT_HSAF/';

L4_path = '/discover/nobackup/projects/land_da/CYGNSS_Experiments/OLv8_M36_cd/';
L4_run = 'OLv8_M36_cd'; 

out_collection = '.tavg24_1d_lnd_Nt.';
domain = 'SMAP_EASEv2_M36_GLOBAL';

out_path = '/discover/nobackup/projects/land_da/Evaluation/IVs/output/';

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
if ~isempty(strfind(domain,'M09'))
    f_EASE = '/css/smapl4/public/L4_Products/L4_SM/Vv7032/lmc/SMAP_L4_SM_lmc_00000000T000000_Vv7032_001.h5';
    lon_L4 = double(h5read(f_EASE,'/cell_lon'));
    lat_L4 = double(h5read(f_EASE,'/cell_lat'));
else
    [lat,lon] = smapeasev2_ind2latlon([0:405],[0:963],'M36');
    lon_L4 = repmat(lon',[1,length(lat)]);
    lat_L4 = repmat(lat,[length(lon),1]);
    clear lat lon
end

f_tilecoord = [L4_path,L4_run,'/output/',domain,'/rc_out/',L4_run,'.ldas_tilecoord.bin'];
tile_coord = read_tilecoord_GEOS(f_tilecoord);
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
    
        
        f_ASCAT = [ASCAT_path, 'H119_H120_processed/Y',num2str(date_time.year,'%4.4d'), ...
            '/M',num2str(date_time.month,'%2.2d'),'/' ...
            'ASCAT_HSAF_H119_SM_',date_string,'_AD.mat']
        
        if exist(f_ASCAT,'file')
            ASC = load(f_ASCAT);

            sm_asc = double(ASC.sm_tile);
            sm_asc(sm_asc > 100) = NaN;

            % aggressive QC: only keep data with conf_flag == 0
            % i.e. remove data with larege noise, wetland frac, top complexity
            % and sensitivity level (see PUM for flag details)
            % the new QC results in ~10% fewer data compared to previous QC
            conf_flag = ASC.conf_flag_tile;
            sm_asc(conf_flag >= 1) = NaN;
            
            clear ASC
            
            sm_obs  = griddata(lon_gpi,lat_gpi,sm_asc,lon_L4,lat_L4,'natural'); clear sm_asc
        else
            sm_obs = NaN *ones(size(lon_L4));
        end

        fname = [L4_path, L4_run,'/output/',domain,'/cat/ens0000/Y', ...
                     num2str(date_time.year, '%4.4d'), ...
                    '/M', num2str(date_time.month, '%2.2d'), ...
                    '/',L4_run, out_collection, ...
                    num2str(date_time.year, '%4.4d'), ...
                    num2str(date_time.month, '%2.2d'), ...
                    num2str(date_time.day,'%2.2d'), '_1200z.nc4'];

        if ~exist(fname,'file')
            fname = strrep(fname,'/ens0000/','/ens_avg/');
        end

        if ~exist(fname,'file')
            fname = [L4_path, L4_run,'/output/',domain,'/cat/ens0000/Y', ...
                num2str(date_time.year, '%4.4d'), ...
                '/M', num2str(date_time.month, '%2.2d'), ...
                '/',L4_run, out_collection, ...
                num2str(date_time.year, '%4.4d'), ...
                num2str(date_time.month, '%2.2d'), ...
                num2str(date_time.day,'%2.2d'), '.nc4'];
        end

        if ~exist(fname,'file')
            fname = strrep(fname,'/ens0000/','/ens_avg/');
        end

        
        if contains(L4_run,'MERRA2_')
            fname = strrep(fname,'.SMAP_L4_SM_gph.','.tavg1_2d_lnd_Nx.');
        end

        if exist(fname,'file')
                disp(fname)
                smsf_tmp = ncread(fname,'SFMC');
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
                
                fname = [L4_path, L4_run,'/output/',domain,'/cat/ens0000/Y', ...
                    num2str(date_time.year, '%4.4d'), ...
                    '/M', num2str(date_time.month, '%2.2d'), ...
                    '/',L4_run, out_collection, ...
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
        
        f_out = [out_path,'ASCL4_H119_SMSF_L4_',L4_run,'_QC_1_',date_string,'.mat'];
        
        save(f_out,'sm_mod','sm_obs','idx_EASEv2_lonxlat');
        
        
        clear sm_obs sm_mod idx_EASEv2_lonxlat
        
   
    date_time = augment_date_time(86400, date_time);
    
end



