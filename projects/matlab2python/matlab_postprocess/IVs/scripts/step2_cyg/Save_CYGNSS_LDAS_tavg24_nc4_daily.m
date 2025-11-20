clear 
addpath /discover/nobackup/qliu/gdelanno_RTM/MATLAB_LDASSA/
addpath /home/qliu/projects/matlab_code/SMAP/

CYG_path = '/discover/nobackup/projects/gmao/smap/SMAP_Nature/CYGNSS';

L4_path = '/discover/nobackup/projects/land_da/CYGNSS_Experiments/DAv8_M36_cd/';
L4_run = 'DAv8_M36_cd'; 

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

f_tilecoord = [L4_path,L4_run,'/output/',domain,'/rc_out/',L4_run,'.ldas_tilecoord.bin'];
tile_coord = read_tilecoord_GEOS(f_tilecoord);
tc = tile_coord;

date_time = start_time;

[latM36, lonM36] = smapeasev2_ind2latlon(0:405, 0:963, 'M36');   % requires your SMAP toolbox path
[lon_L4, lat_L4] = ndgrid(double(lonM36), double(latM36));  

while true
    
    % Build date strings early (used many times)
    yyyy = num2str(date_time.year,'%04d');
    mm   = num2str(date_time.month,'%02d');
    dd   = num2str(date_time.day,'%02d');
    date_string = [yyyy, mm, dd];

        
% --- Build CYGNSS filename for this day ---
% assumes variables: yyyy, mm, dd as strings (e.g., '2018','08','03')
    cyg_fn = fullfile(CYG_path, ['Y' yyyy], ['M' mm], ...
    sprintf('cyg.ddmi.s%s%s%s-030000-e%s%s%s-210000.l3.grid-soil-moisture-36km.a32.d33.nc', ...
            yyyy, mm, dd, yyyy, mm, dd));

    if ~isfile(cyg_fn)
        % If missing, fill obs with NaNs on EASE grid
        sm_obs_2d = nan(964, 406, 'single');
    else
        % --- Read CYGNSS subdaily (timeslices x lat x lon) ---
        % --- Read native lon/lat (curvilinear) and cast to double ---
        lon2d = double(ncread(cyg_fn, 'longitude'));   % [lat,lon]
        lat2d = double(ncread(cyg_fn, 'latitude'));    % [lat,lon]
        if min(lon2d(:)) >= 0 && max(lon2d(:)) > 180, lon2d = lon2d - 360; end

        % --- Read SM_subdaily (unknown dimension order); cast to double ---
        sm_raw = double(ncread(cyg_fn, 'SM_subdaily'));   % dims could be [T,LAT,LON] or [LAT,LON,T] etc.
        sz = size(sm_raw);
        tgt = size(lat2d);  % [nlat nlon], should equal size(lon2d)

        % Figure out which two dims match (nlat,nlon), and which is time
        % Try all permutations to find [LAT LON T] layout
        perm_found = [];
        perms3 = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];
        for ii = 1:6
            p = perms3(ii,:);
            tmp = permute(sm_raw, p);   % size becomes sz(p)
            if isequal(size(tmp,1), tgt(1)) && isequal(size(tmp,2), tgt(2))
                perm_found = p;  % tmp is [LAT LON T]
                break
            end
        end
        if isempty(perm_found)
            error('SM_subdaily dims %s do not align with lat/lon dims %dx%d.', mat2str(sz), tgt(1), tgt(2));
        end

        % Reorder to [LAT LON T] and average over T
        sm_sub = permute(sm_raw, perm_found);          % -> [LAT LON T]
        clear sm_raw
        % QC (range + sigma) on the reordered array
        vmin = 0.0; vmax = 1.0;
        try, vmin = ncreadatt(cyg_fn,'SM_subdaily','valid_min'); end
        try, vmax = ncreadatt(cyg_fn,'SM_subdaily','valid_max'); end

        % Optional sigma handling: read and permute to match [LAT LON T] if present
        has_sigma = false;
        try
            sig_raw = double(ncread(cyg_fn,'SIGMA_subdaily'));   % unknown order too
            sig_sub = permute(sig_raw, perm_found);              % align to [LAT LON T]
            has_sigma = true; clear sig_raw
        catch
        end

        bad_range = (sm_sub < vmin) | (sm_sub > vmax) | ~isfinite(sm_sub);
        if has_sigma
            bad_sigma = ~isfinite(sig_sub);  % or sig_sub > 0.08 (tune)
            sm_sub(bad_sigma) = NaN;
        end
        sm_sub(bad_range) = NaN;

        % Daily mean & contributing-count on native grid (LAT×LON)
        sm_obs_2d = mean(sm_sub, 3, 'omitnan');        % size(sm_obs_2d) == size(lat2d)
        Nsub      = uint8(sum(isfinite(sm_sub), 3));   % 0..T slices
        clear sm_sub

        % --- Now do your griddata to EASE-M36 (lon_L4/lat_L4 are 964×406) ---
        sm_obs_2d = double(sm_obs_2d);
        sm_obs_L4 = griddata(lon2d, lat2d, sm_obs_2d, double(lon_L4), double(lat_L4), 'natural');
        sm_obs_2d = single(sm_obs_L4);  % keep memory modest
        clear sm_obs_L4
        
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
    
    sm_mod_glob_2d = NaN * sm_obs_2d;

    for k = 1:length(sm_mod)
        sm_mod_glob_2d(tc.i_indg(k)+1,tc.j_indg(k)+1) = sm_mod(k);
    end
    
    clear sm_mod
    sm_mod = sm_mod_glob_2d; clear sm_mod_glob_2d

    clear sm_obs
    sm_obs = sm_obs_2d; clear sm_obs_2d

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
    
    f_out = [out_path,'CYGL3_SMSF_',L4_run,'_QC_1_',date_string,'.mat'];
    
    save(f_out,'sm_mod','sm_obs','idx_EASEv2_lonxlat', '-v7.3');
    
    
    clear sm_obs sm_mod idx_EASEv2_lonxlat
    
     % --- exit if this was the end date; else advance one day ---
     if (date_time.year == end_time.year && date_time.month == end_time.month && date_time.day == end_time.day)
        break
    end   

    date_time = augment_date_time(86400, date_time);
    
end



