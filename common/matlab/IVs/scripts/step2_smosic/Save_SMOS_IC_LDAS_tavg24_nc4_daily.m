clear
addpath /discover/nobackup/qliu/gdelanno_RTM/MATLAB_LDASSA/
addpath /home/qliu/projects/matlab_code/SMAP/

SMOSIC_nc_path = '/discover/nobackup/projects/land_da/SMOS_IC/preprocessed_m36_daily/';

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

[lat,lon] = smapeasev2_ind2latlon(0:405, 0:963, 'M36');
lon_L4 = repmat(lon',[1,length(lat)]);
lat_L4 = repmat(lat,[length(lon),1]); %#ok<NASGU>
Nlon = size(lon_L4,1);
Nlat = size(lon_L4,2);
clear lat lon lon_L4

f_tilecoord = [L4_path,L4_run,'/output/',domain,'/rc_out/',L4_run,'.ldas_tilecoord.bin'];
tile_coord = read_tilecoord_GEOS(f_tilecoord);
tc = tile_coord;

date_time = start_time;

while true

    date_string = [num2str(date_time.year,'%4.4d'), ...
        num2str(date_time.month,'%2.2d'), num2str(date_time.day,'%2.2d')];

    % --- Read SMOS-IC daily sparse cache on M36 ---
    f_obs_nc = [SMOSIC_nc_path, 'smos_ic_sm_m36_', date_string, '.nc'];
    sm_obs = NaN(Nlon, Nlat);

    if exist(f_obs_nc, 'file')
        idx0 = int64(ncread(f_obs_nc, 'idx_EASEv2_lonxlat'));
        obs_vals = double(ncread(f_obs_nc, 'sm_obs'));
        lin = double(idx0(:)) + 1; % 0-based -> MATLAB 1-based linear indexing
        obs_vals = obs_vals(:);
        good = lin >= 1 & lin <= (Nlon * Nlat) & isfinite(obs_vals);
        sm_obs(lin(good)) = obs_vals(good);
    end

    % --- Read daily model SFMC (same logic as ASCAT/CYGNSS scripts) ---
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
        % Guard SFMC orientation: require [tile x time] before averaging.
        if isvector(smsf_tmp)
            smsf_tmp = smsf_tmp(:);
        else
            [n1, n2] = size(smsf_tmp);
            if n1 == tc.N_tile
                % already [tile x time]
            elseif n2 == tc.N_tile
                smsf_tmp = smsf_tmp.';
            else
                error('Unexpected SFMC shape [%d x %d], cannot map to N_tile=%d.', n1, n2, tc.N_tile)
            end
        end
        sm_mod_3hr = smsf_tmp; clear smsf_tmp

    else

        hr_list = 1:3:22;

        for this_hr = 1:8
            date_time_string = [num2str(date_time.year, '%4.4d'), ...
                num2str(date_time.month, '%2.2d'), ...
                num2str(date_time.day, '%2.2d'), ...
                '_', num2str(hr_list(this_hr), '%2.2d'), ...
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
    end

    sm_mod = nanmean(sm_mod_3hr,2); clear sm_mod_3hr
    if length(sm_mod) ~= tc.N_tile
        error('something is wrong reading LDAS data')
    end

    sm_mod_glob_2d = NaN(Nlon, Nlat);
    for k = 1:length(sm_mod)
        sm_mod_glob_2d(tc.i_indg(k)+1, tc.j_indg(k)+1) = sm_mod(k);
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

    f_out = [out_path,'SMOSIC_SMSF_MOD_',L4_run,'_QC_1_',date_string,'.mat'];
    save(f_out,'sm_mod','sm_obs','idx_EASEv2_lonxlat','-v7.3');

    clear sm_obs sm_mod idx_EASEv2_lonxlat

    if (date_time.year == end_time.year && ...
            date_time.month == end_time.month && ...
            date_time.day == end_time.day)
        break
    end

    date_time = augment_date_time(86400, date_time);
end
