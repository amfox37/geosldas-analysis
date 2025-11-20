% This script computes the stats against sparse network INSITU data

clear

% addpath('/discover/nobackup/amfox/metop_sm_GEOSldas/GEOSldas/src/Applications/LDAS_App/util/shared/matlab/')
% addpath('/gpfsm/dnb34/amfox/GEOSldas_diagnostics/Jupyter/Matlab_functions')
addpath /discover/nobackup/qliu/gdelanno_RTM/MATLAB_LDASSA/
addpath /home/qliu/projects/matlab_code/SMAP/
addpath /home/amfox/for_Andy

INSITU_tag ='SCAN'; % 'USCRN', 'Oznet', 'Msnet' , 'SMOSM', 'CalVal_M33', 'SCAN'

add_anomR = 1;  % no need to compute anom R

fout_path = '/discover/nobackup/projects/land_da/Evaluation/InSitu/output/';    

if contains(INSITU_tag, 'CalVal')
    rzmc_tag = '';
elseif strcmp(INSITU_tag, 'Oznet') || strcmp(INSITU_tag, 'SMOSM')
    rzmc_tag = 'c3smv';
else
    if  strcmp(INSITU_tag, 'Msnet')
        rzmc_tag = 'c123smv';
    else
        rzmc_tag = 'c1234smv';
    end
end

% exp_path = ['/discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_OLv8_M36/'];
exp_path = ['/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v2/'];

exp_run_name= {...
        'LS_OLv8_M36'; ...
        'LS_DAv8_M36'
     };

kkk = 2;

% tmp mat file to store the time series
fout_name = [fout_path,'/',exp_run_name{kkk},'_',INSITU_tag, '_SM_1d_',rzmc_tag,'_9yr_']

if exist([fout_name,'raw_timeseries.mat'],'file')
    
    load([fout_name,'raw_timeseries.mat'])
    
else
    
    % read INSITU coord file
    [INSITU_lat, INSITU_lon, INSITU_id] = get_INSITU_coord(INSITU_tag);
    
    if strcmp(INSITU_tag, 'SCAN')
        INSITU_path = '/discover/nobackup/qliu/merra_land/DATA/SCAN/data/PROCESSED_202501_QL/';
    elseif strcmp(INSITU_tag, 'USCRN')
        INSITU_path = '/discover/nobackup/qliu/merra_land/DATA/USCRN/data/PROCESSED_202501/';
    elseif strcmp(INSITU_tag, 'Oznet')
        INSITU_path = '/discover/nobackup/qliu/merra_land/DATA/Murrumbidgee/QLiu_202211/data/PROCESSED/';
    elseif strcmp(INSITU_tag, 'Msnet')
        INSITU_path = '/discover/nobackup/qliu/merra_land/DATA/Oklahoma_Mesonet/data/PROCESSED_201805/';
    elseif strcmp(INSITU_tag, 'SMOSM')
        INSITU_path = '/discover/nobackup/qliu/merra_land/DATA/ISMN/data/PROCESSED_20230101/SMOSMANIA/';
    elseif contains(INSITU_tag, 'CalVal')
        INSITU_path = '/discover/nobackup/qliu/merra_land/DATA/CalVal_insitu/SMAP_refpix/202104/';
    else
        error(['do not recognize ', INSITU_tag ]);
    end
    


    exp_run= exp_run_name{kkk};

    domain = 'SMAP_EASEv2_M36_GLOBAL';

    out_collection = '.tavg24_1d_lnd_Nt.';

    start_time.year  = 2015;
    start_time.month = 4;
    start_time.day   = 1;
    start_time.hour  = 12;
    start_time.min   = 0;
    start_time.sec   = 0;
    
    end_time.year  =2024;
    end_time.month = 6; 
    end_time.day   = 1;
    end_time.hour  = start_time.hour;
    end_time.min   = start_time.min;
    end_time.sec   = start_time.sec;
    
    start_time = get_dofyr_pentad(start_time);
    end_time   = get_dofyr_pentad(end_time);
    
    file_tag = 'ldas_tile_daily_out';
    
    if contains(file_tag, 'daily')
        dtstep = 86400;
    else
        dtstep = 10800;
    end
    
    % file INSITU site corresponded tile index
    
    fname_tilecoord = [ exp_path, '/', exp_run, '/output/', domain, '/rc_out/', ...
        exp_run, '.ldas_tilecoord.bin'];
    
    tile_coord = read_tilecoord(fname_tilecoord);
    
    tc = tile_coord;
    
    max_distance = .1;
    
    clear ind_tile
    
    for i=1:length(INSITU_id)
        
        this_lat = INSITU_lat(i);
        this_lon = INSITU_lon(i);
        
        distance = (( this_lat - tc.com_lat) .^2 + ...
            ( this_lon - tc.com_lon) .^2     );
        
        ind_tile(i) = find(min(distance) == distance);
        
        distance_min(i) = distance(ind_tile(i));
        
        if distance(ind_tile(i)) > max_distance
            
            disp(['could not find tile for In Situ station ', INSITU_id(i)])
            disp(['distance=', num2str(distance(ind_tile(i))), ...
                ' exceeds max_distance=', num2str(max_distance) ])
            
            ind_tile(i) = NaN;
            
        end
    end
    
    % eliminate In Situ stations for which tile could not be found
    
    INSITU_id(    isnan(ind_tile) ) = [];
    INSITU_lat(   isnan(ind_tile) ) = [];
    INSITU_lon(   isnan(ind_tile) ) = [];
    distance_min( isnan(ind_tile) ) = [];
    
    ind_tile(   isnan(ind_tile) ) = [];
    
    % -----------------------------------------------------------------------
    % read time series of tavg fields
    clear tile_coord tc LDAS_sm_org LDAS_prcp_org
    
        
        % re-read tilecoord file for each run in case the tile orders are different
        
        fname_tilecoord = [ exp_path,'/', exp_run, '/rc_out/',exp_run, '.ldas_tilecoord.bin'];
        
        if ~exist(fname_tilecoord,'file')
            fname_tilecoord = [ exp_path, '/', exp_run, '/output/', domain, '/rc_out/', ...
                exp_run, '.ldas_tilecoord.bin'];
        end
        
        clear tile_coord tc
        
        tile_coord = read_tilecoord(fname_tilecoord);
        
        tc = tile_coord;
        
        clear ind_tile
        
        for i=1:length(INSITU_id)
            
            this_lat = INSITU_lat(i);
            this_lon = INSITU_lon(i);
            
            distance = (( this_lat - tc.com_lat) .^2 + ...
                ( this_lon - tc.com_lon) .^2     );
            
            ind_tile(i) = find(min(distance) == distance);
            
        end
        
        date_time = start_time;
        
        t_ind = 0;
        
        while 1
            
            if (date_time.year ==end_time.year   && ...
                    date_time.month==end_time.month  && ...
                    date_time.day  ==end_time.day    )
                break
            end
            
            t_ind = t_ind + 1;
            
            date_time_vec(t_ind) = date_time;
            date_time_string = get_date_time_string(date_time,file_tag);
            
            % load LDAS ens_avg
            
            file_ext = '.nc4';
            
            fname = [exp_path, exp_run, '/output/', domain,          ...
                    '/cat/ens_avg',                                       ...
                    '/Y',   num2str(date_time.year, '%4.4d'),         ...
                    '/M',   num2str(date_time.month,'%2.2d'),         ...
                    '/', exp_run, out_collection,                          ...
                    date_time_string, '_1200z',file_ext];
                
            if ~exist(fname,'file')
                fname = strrep(fname,'/ens_avg/','/ens0000/');
            end

            if ~exist(fname,'file')
                fname = [exp_path, '/', exp_run, '/output/', domain,         ...
                        '/cat/ens_avg',                                       ...
                        '/Y',   num2str(date_time.year, '%4.4d'),         ...
                        '/M',   num2str(date_time.month,'%2.2d'),         ...
                        '/', exp_run, out_collection,                         ...
                        date_time_string(1:8), file_ext];
            end

            if ~exist(fname,'file')
                fname = strrep(fname,'/ens_avg/','/ens0000/');
            end
            % read data
            disp(fname)
            clear tile_data_tmp
            tile_data_tmp = ncread(fname, 'SFMC');
            if size(tile_data_tmp,2) > 1 && dtstep == 10800
                 ind_hour = ceil((date_time.hour+1)/3.);
            else
                 ind_hour = 1;                    
            end
            LDAS_sm_org(t_ind,1,:) = tile_data_tmp(ind_tile,ind_hour);  %sfmc
            clear tile_data_tmp
            
            tile_data_tmp = ncread(fname, 'RZMC');
            LDAS_sm_org(t_ind,2,:) = tile_data_tmp(ind_tile,ind_hour); 
            clear tile_data_tmp

            tile_data_tmp = ncread(fname, 'TSURFLAND');
            LDAS_st_org(t_ind,1,:) = tile_data_tmp(ind_tile,ind_hour);
            clear tile_data_tmp

            tile_data_tmp = ncread(fname, 'TSOIL1');
            LDAS_st_org(t_ind,2,:) = tile_data_tmp(ind_tile,ind_hour);
            clear tile_data_tmp

            tile_data_tmp = ncread(fname,'PRECTOTCORRLAND');
            LDAS_prcp_org(t_ind,1,:) = tile_data_tmp(ind_tile,ind_hour); 
            clear tile_data_tmp
                
            date_time = augment_date_time(dtstep, date_time);
            clear fname
        end
        
    
    % read INSITU data
    
    clear INSITU_st INSITU_sm date_time_vec1 INSITU_prcp

    % For CalVal daily: keep a reference [Y M D] axis from the first station
    ref_YMD = [];
    ref_T   = [];
    
    for i=1:length(INSITU_id)
        
        INSITU_id_string = INSITU_id{i};
        
        if strcmp(INSITU_tag, 'USCRN')
            tag1 = [INSITU_tag,'_'];
            tmp_tag = '_2009_2024';
        elseif strcmp(INSITU_tag, 'Oznet')
            tag1 = '';
            tmp_tag = '_2015_2022';
        elseif strcmp(INSITU_tag, 'Msnet')
            tag1 = ['Mesonet','_'];
            tmp_tag = '_201504_201804';
        elseif strcmp(INSITU_tag, 'SMOSM')
            tag1 = '';
            tmp_tag = '_2010_2022';
        elseif contains(INSITU_tag, 'CalVal')
            tag1 = '';
            tmp_tag = '_201504_202103';
        else
            tag1 = [INSITU_tag,'_'];
            tmp_tag = '';
        end
        
        if contains(INSITU_tag, 'CalVal')
            fname_tmp = [INSITU_path, '/', tag1,INSITU_id_string,tmp_tag];
        else
            fname_tmp = [INSITU_path, '/',INSITU_id{i},'/', tag1,INSITU_id_string,tmp_tag];
        end
        clear tmp_tag tag1

        if contains(INSITU_tag, 'CalVal')
            data_ext = '1h';
        else
            if contains(file_tag, 'daily')
                data_ext = '1d';
            else
                data_ext = '3h';
            end
        end        
        
        fname = [fname_tmp,'_',data_ext, '.mat'];
        disp(['loading ', fname])
        load(fname);
        
        % --- load/prepare Sdata and map into INSITU_* arrays ---
        if contains(file_tag, 'daily')
            
            if strcmp(INSITU_tag, 'USCRN')
                Sdata = Udata_1d; clear Udata_1d

            elseif contains(INSITU_tag, 'CalVal')
                % CalVal hourly -> aggregate to DAILY
                % regular_data columns expected:
                % [year, month, day, hour, min, sec, doy, pentad, prcp, sm1, sm2, tempC, ...]
                S = regular_data;

                % Group by calendar day
                Y = S(:,1); M = S(:,2); D = S(:,3);
                G = findgroups(Y, M, D);

                % Daily SM means (two SM cols), daily Temp mean (°C), daily Precip sum
                sm1_daily = splitapply(@(x) mean(x, 'omitnan'), S(:,10), G);
                sm2_daily = splitapply(@(x) mean(x, 'omitnan'), S(:,11), G);

                tC_daily  = splitapply(@(x) mean(x, 'omitnan'), S(:,12), G);

                pr = S(:,9); pr(pr < -1e-5) = NaN;                         % QC: negatives -> NaN
                pr_daily  = splitapply(@(x) sum(x, 'omitnan'), pr, G);

                % Day keys (one per group)
                Yd = splitapply(@(x) x(1), Y, G);
                Md = splitapply(@(x) x(1), M, G);
                Dd = splitapply(@(x) x(1), D, G);

                % Build daily time metadata (set 12:00 to match model's daily timestamp convention)
                dt  = datetime(Yd, Md, Dd, 12, 0, 0);
                doy = day(dt, 'dayofyear');
                pentad = floor((doy - 1) / 5) + 1;

                % Recreate Sdata in the same column schema the rest of the code expects
                % cols: [year month day hour min sec doy pentad prcp sm1 sm2 tempC]
                Sdata = [Yd, Md, Dd, repmat(12, size(Yd)), zeros(size(Yd)), zeros(size(Yd)), ...
                         doy, pentad, pr_daily, sm1_daily, sm2_daily, tC_daily];

                clear S Y M D G sm1_daily sm2_daily tC_daily pr pr_daily Yd Md Dd dt doy pentad

            else
                Sdata = Sdata_1d; clear Sdata_1d
            end

        else
            % ----- Sub-daily (e.g., 3-hourly) path: unchanged -----
            if strcmp(INSITU_tag, 'USCRN')
                Sdata = Udata_3h; clear Udata_3h
            elseif contains(INSITU_tag, 'CalVal')
                Sdata = regular_data;   % hourly CalVal used as-is for sub-daily runs
            else
                Sdata = Sdata_3h; clear Sdata_3h
            end
        end

        % Map Sdata into unified INSITU_* arrays (by network conventions)
        if strcmp(INSITU_tag, 'Msnet')
            INSITU_sm(:,:,i) = Sdata(:,10:12);
            INSITU_st(:,:,i) = Sdata(:,13:15);

        elseif strcmp(INSITU_tag, 'SMOSM')
            INSITU_sm(:,:,i) = Sdata(:, 9:12);
            INSITU_st(:,:,i) = Sdata(:,13:16);

        elseif contains(INSITU_tag, 'CalVal') && contains(file_tag, 'daily')
            % CalVal DAILY: reindex each station to the FIRST station's [Y M D] axis
        
            if isempty(ref_YMD)
                % First station sets the reference calendar days
                ref_YMD = Sdata(:,1:3);             % [year month day]
                ref_T   = size(ref_YMD,1);
        
                INSITU_sm(:,:,i)   = Sdata(:,10:11);
                INSITU_st(:,:,i)   = Sdata(:,12) + 273.15;
                INSITU_prcp(:,:,i) = Sdata(:,9);
        
            else
                % Map this station's days to the reference axis
                this_YMD = Sdata(:,1:3);
                [tf, loc] = ismember(this_YMD, ref_YMD, 'rows');
        
                sm_pad = NaN(ref_T, 2);
                tK_pad = NaN(ref_T, 1);
                pr_pad = NaN(ref_T, 1);
        
                sm_pad(loc(tf), :) = Sdata(tf, 10:11);
                tK_pad(loc(tf))    = Sdata(tf, 12) + 273.15;
                pr_pad(loc(tf))    = Sdata(tf, 9);
        
                INSITU_sm(:,:,i)   = sm_pad;
                INSITU_st(:,:,i)   = tK_pad;
                INSITU_prcp(:,:,i) = pr_pad;
        
                clear this_YMD tf loc sm_pad tK_pad pr_pad
            end
        
        elseif contains(INSITU_tag, 'CalVal')
            % CalVal SUB-DAILY (e.g., hourly used for 3-hourly runs): original mapping
            INSITU_sm(:,:,i)   = Sdata(:,10:11);
            INSITU_st(:,:,i)   = Sdata(:,12) + 273.15;
            INSITU_prcp(:,:,i) = Sdata(:,9);
        

        else
            INSITU_sm(:,:,i) = Sdata(:,10:14);
            INSITU_st(:,:,i) = Sdata(:,15:19);
        end

        % Build/refresh the observation time vector on first station only
        if i == 1
            date_time_vec1.year   = Sdata(:, 1);
            date_time_vec1.month  = Sdata(:, 2);
            date_time_vec1.day    = Sdata(:, 3);
            date_time_vec1.hour   = Sdata(:, 4);
            date_time_vec1.min    = Sdata(:, 5);
            date_time_vec1.sec    = Sdata(:, 6);
            date_time_vec1.doy    = Sdata(:, 7);
            date_time_vec1.pentad = Sdata(:, 8);
        end

        clear Sdata
    end    
    
% --- REPLACED BLOCK: align model/obs by true overlap (no forced start/end) ---

    % Build datetime arrays from the existing time vectors
    model_dt = datetime( ...
        [ [date_time_vec.year]'  [date_time_vec.month]'  [date_time_vec.day]' ...
          [date_time_vec.hour]'  [date_time_vec.min]'    [date_time_vec.sec]' ], ...
        'Format','yyyy-MM-dd HH:mm:ss');

    obs_dt = datetime( ...
        date_time_vec1.year, date_time_vec1.month, date_time_vec1.day, ...
        date_time_vec1.hour, date_time_vec1.min,   date_time_vec1.sec, ...
        'Format','yyyy-MM-dd HH:mm:ss');

    if contains(file_tag, 'daily')
        % Match by calendar day (robust even if hours differ, e.g., 00 vs 12)
        model_key = dateshift(model_dt, 'start', 'day');
        obs_key   = dateshift(obs_dt,   'start', 'day');
    else
        % Sub-daily: match exact timestamps.
        % For CalVal (hourly), bin to 3-hour boundaries to match LDAS cadence.
        model_key = model_dt;
        if contains(INSITU_tag, 'CalVal')
            obs_key = dateshift(obs_dt, 'start', 'hour');
            obs_key.Minute(:) = 0;  obs_key.Second(:) = 0;
            % Floor to 3-hour bins anchored at midnight: 0,3,6,...,21
            obs_key.Hour = obs_key.Hour - mod(obs_key.Hour, 3);
        else
            obs_key = obs_dt;
        end
    end

    % Keep only timestamps present in BOTH series (preserve model order)
    [common_key, im, io] = intersect(model_key, obs_key, 'stable');

    if isempty(common_key)
        error('No overlapping timestamps between model and %s observations after alignment.', INSITU_tag);
    end

    % Subset the MODEL arrays to overlap
    date_time_vec = date_time_vec(im);             % keep model time metadata in sync
    LDAS_sm_org   = LDAS_sm_org(im,:,:);
    if exist('LDAS_st_org','var'),   LDAS_st_org   = LDAS_st_org(im,:,:);   end
    if exist('LDAS_prcp_org','var'), LDAS_prcp_org = LDAS_prcp_org(im,:,:); end

    % Subset the OBS arrays to overlap
    INSITU_sm_tmp = INSITU_sm(io,:,:);
    INSITU_st_tmp = INSITU_st(io,:,:);
    if exist('INSITU_prcp','var'), INSITU_prcp_tmp = INSITU_prcp(io,:,:); end

    % CalVal + 3-hourly: apply a 3-point temporal smoothing (centered)
    if contains(INSITU_tag,'CalVal') && (dtstep == 10800)
        INSITU_sm_tmp = movmean(INSITU_sm_tmp, 3, 1, 'omitnan');
        INSITU_st_tmp = movmean(INSITU_st_tmp, 3, 1, 'omitnan');
        if exist('INSITU_prcp_tmp','var')
            INSITU_prcp_tmp = movmean(INSITU_prcp_tmp, 3, 1, 'omitnan');
            INSITU_prcp_tmp(INSITU_prcp_tmp < -1e-5) = NaN;
            INSITU_prcp = INSITU_prcp_tmp;
            clear INSITU_prcp_tmp
        end
    end

    INSITU_sm_tmp(INSITU_st_tmp < 277.16) = NaN;
    INSITU_sm_tmp(INSITU_sm_tmp < 1e-4) = NaN;

    % Finalize obs arrays
    INSITU_sm = INSITU_sm_tmp;
    INSITU_st = INSITU_st_tmp;
    clear INSITU_sm_tmp INSITU_st_tmp common_key im io model_dt obs_dt model_key obs_key
    
    % -----------------------------------------------------------
    % add fields to INSITU_sm
    
    % c12smv  (weighted average of c1smv and c2smv)
    % c123smv (weighted average of c1smv and c2smv and c3smv)
    % ...
    
    if strcmp(INSITU_tag, 'Msnet')
        tavg_tag_INSITU = {'c1smv';'c2smv';'c3smv'};
        
        ind_tmp = [1  2  3  ];    % c1smv, c2smv, ..., c5smv
        w_tmp   = [2  10 24  ];    % corresponding layer depths in inches
        
        tmpstr  = '123';
    elseif strcmp(INSITU_tag, 'SMOSM')
        tavg_tag_INSITU = {'c1smv';'c2smv';'c3smv';'c4smv'};
        
        ind_tmp = [1  2  3  4];    % c1smv, c2smv, ..., c5smv
        w_tmp   = [3  3  8 16 ];    % corresponding layer depths in inches
        
        tmpstr  = '1234';
        
    else
        
        tavg_tag_INSITU = {'c1smv';'c2smv';'c3smv';'c4smv';'c5smv'};
        
        ind_tmp = [1 2 3   4 5 ];    % c1smv, c2smv, ..., c5smv
        w_tmp   = [  3  3  8 16 10 ];    % corresponding layer depths in inches
        
        tmpstr  = '12345';
        
    end
    
    if ~contains(INSITU_tag, 'CalVal')
        w_mat_tmp = ones( size(INSITU_sm,1), 1) * w_tmp;
        
        for k=1:size(INSITU_sm,3)
            w_mat_tmp3d(:,:,k) = w_mat_tmp;
        end
        
        for j=2:length(ind_tmp)
            % compute weighted average over cXsmv
            INSITU_sm(:,end+1,:) = ...
                sum( INSITU_sm(:,ind_tmp(1:j),:).*w_mat_tmp3d(:,1:j,:), 2) / ...
                sum(w_tmp(1:j));
            
            tavg_tag_INSITU{end+1} = [ 'c', tmpstr(1:j), 'smv'];
        end
        
        % only keep surface and selected rootzone SM for INSITU_sm
        % the following 2 lines should be changed together
        
        if strcmp(INSITU_tag, 'Oznet') || strcmp(INSITU_tag, 'SMOSM')
            rzmc_tag = 'c3smv';
        else
            if contains(file_tag, 'daily')
                rzmc_tag = 'c1234smv';
            else
                if strcmp(INSITU_tag, 'Msnet')
                    rzmc_tag = 'c123smv';
                else
                    rzmc_tag = 'c1234smv';
                end
                %rzmc_tag = 'c3smv';
            end
        end
        
        ind_tmp = find(strcmp(tavg_tag_INSITU, rzmc_tag));
        
        INSITU_sm(:,setxor([1 ind_tmp], [1:size(INSITU_sm,2)]),:) = [];
        
    end
    
    %INSITU_sm = permute(INSITU_sm, [1, 3,2]);
    
    % Minum points for statistics (R, RMSE, etc)
    
    if contains(file_tag, 'daily')
        Nmin = 200;
    else
        Nmin = 480;
    end
    
    save([fout_name,'raw_timeseries.mat'])
    
end

for it = 1:length(date_time_vec)
    doy_vec(it) = date_time_vec(it).dofyr;
    year_vec(it) = date_time_vec(it).year;
    month_vec(it) = date_time_vec(it).month;
    day_vec(it) = date_time_vec(it).day;
end

% Compute SM skills

nv = 2;
nn = size(LDAS_sm_org,3);
nc = size(LDAS_sm_org,4);

R   = NaN*ones(nn,nv,nc);
RLO = NaN*ones(nn,nv,nc) ;
RUP = NaN*ones(nn,nv,nc);

if add_anomR
    anomR   = NaN*ones(nn,nv,nc);
    anomRLO = NaN*ones(nn,nv,nc);
    anomRUP = NaN*ones(nn,nv,nc);
end

Bias   = NaN*ones(nn,nv,nc);
BiasLO = NaN*ones(nn,nv,nc);
BiasUP = NaN*ones(nn,nv,nc);

absBias = NaN*ones(nn,nv,nc);
absBiasLO = NaN*ones(nn,nv,nc);
absBiasUP = NaN*ones(nn,nv,nc);

RMSE =   NaN*ones(nn,nv,nc);
ubRMSE =  NaN*ones(nn,nv,nc);
ubRMSELO = NaN*ones(nn,nv,nc);
ubRMSEUP = NaN*ones(nn,nv,nc);

for kk = kkk % 1:length(exp_run)
    
    clear LDAS_sm
    
    LDAS_sm(:,1,:) = LDAS_sm_org(:,1,:);
    LDAS_sm(:,2,:) = LDAS_sm_org(:,2,:);  %rzmc
    LDAS_sm(LDAS_sm == 0) = NaN;
    
    %Cross mask data
    LDAS_sm(isnan(INSITU_sm)) = NaN;
    INSITU_sm(isnan(LDAS_sm)) = NaN;
    
    for i = 1:size(LDAS_sm,3)
        for j = 1:size(LDAS_sm,2)
            
            tmpdata = [INSITU_sm(:,j,i), LDAS_sm(:,j,i)];
            [ stats, stats_tags ] = get_validation_stats( tmpdata, 1 , 'complete', ...
                1, [1:2], Nmin, 1 );
            
            N_data = sum(all(~isnan(tmpdata),2));
            nodata_val = -9999;
            nodata_tol =  1E-4;
            
            stats(abs(stats-nodata_val)<nodata_tol ) = NaN;
            stats(stats == 0 ) = NaN;
            
            R(i,j)   = stats(strcmp(stats_tags,'R'));
            RLO(i,j) = stats(strcmp(stats_tags,'RLO')) - R(i,j) ;
            RUP(i,j) = stats(strcmp(stats_tags,'RUP')) - R(i,j);
            
            Bias(i,j)   = -stats(strcmp(stats_tags,'bias'));
            BiasLO(i,j) = -stats(strcmp(stats_tags,'CI_bias')) ;
            BiasUP(i,j) =  stats(strcmp(stats_tags,'CI_bias')) ;
            
            absBias(i,j) =  abs(stats(strcmp(stats_tags,'bias')));
            absBiasLO(i,j) = -stats(strcmp(stats_tags,'CI_bias')) ;
            absBiasUP(i,j) =  stats(strcmp(stats_tags,'CI_bias')) ;
            
            RMSE(i,j) =     sqrt( stats(strcmp(stats_tags,'MSE')) );
            RMSELO(i,j) = -sqrt(stats(strcmp(stats_tags,'CI_MSE'))) ;
            RMSEUP(i,j) =  sqrt(stats(strcmp(stats_tags,'CI_MSE'))) ;
            ubRMSE(i,j) =   sqrt(stats(strcmp(stats_tags,'MSE')) - stats(strcmp(stats_tags,'bias')).^2);
            ubRMSELO(i,j) = -sqrt(stats(strcmp(stats_tags,'CI_MSE'))+stats(strcmp(stats_tags,'MSE')) - ...
                stats(strcmp(stats_tags,'bias')).^2)+ubRMSE(i,j);
            ubRMSEUP(i,j) =  sqrt(stats(strcmp(stats_tags,'CI_MSE'))+stats(strcmp(stats_tags,'MSE')) - ...
                stats(strcmp(stats_tags,'bias')).^2)-ubRMSE(i,j);
            
            if add_anomR
                
                insitu_data = INSITU_sm(:,j,i);
                model_data  = LDAS_sm(:,j,i);
                
                % cross mask data
                insitu_data(isnan(model_data)) = NaN;
                model_data(isnan(insitu_data)) = NaN;
                
                % window of temporal smoothing
                Nday_window = 31; % days
                Nday_shift = (Nday_window-1)/2;
                
                % minumum data requirement for computing the daily climatology
                Nmin_day = 150 % 240 ; % 8 perday *30 days for one clim window
                
                disp('computing anomalies...')
                
                clear dofyr_list
                for dofyr = 1:365
                    if dofyr <= 15
                        dofyr_list(:,dofyr) = [1:dofyr+Nday_shift 365-(Nday_shift-dofyr):365];
                    elseif dofyr >= 351
                        dofyr_list(:,dofyr) = [dofyr-Nday_shift:365 1:Nday_shift-(365-dofyr)];
                    else
                        dofyr_list(:,dofyr) = [dofyr-Nday_shift : dofyr+Nday_shift];
                    end
                end
                
                for doy = 1:365
                    
                    ind = [];
                    for id = 1:length(dofyr_list(:,doy))
                        ind = [ind, find(doy_vec(:) == dofyr_list(id,doy))'];
                    end
                    
                    tmp = insitu_data(ind);
                    tmp = tmp(~isnan(tmp));
                    if length(tmp) >= Nmin_day
                        tmp = mean(tmp);
                    else
                        tmp = NaN;
                    end
                    insitu_clim(doy) = tmp; clear tmp
                    
                    tmp = model_data(ind);
                    tmp = tmp(~isnan(tmp));
                    if length(tmp) >= Nmin_day
                        tmp = mean(tmp);
                    else
                        tmp = NaN;
                    end
                    model_clim(doy) = tmp; clear tmp
                    
                end
                
                model_anom = NaN * ones(size(model_data));
                insitu_anom = NaN * ones(size(model_data));
                
                for doy = 1:366
                    ind = find(doy_vec(:) == doy);
                    if ~isempty(ind)
                        insitu_anom(ind) = insitu_data(ind) - insitu_clim(min(365,doy));
                        model_anom(ind) = model_data(ind) - model_clim(min(365,doy));
                    end
                end
                
                clear stats stats_tags
                [stats, stats_tags] = get_validation_stats(...
                    [insitu_anom model_anom], 1, 'complete', 1, [1:2], Nmin, 1 );
                
                nodata_val = -9999;
                nodata_tol =  1E-4;
                
                stats(abs(stats-nodata_val)<nodata_tol ) = NaN;
                stats(stats == 0 ) = NaN;
                
                
                anomR(i,j)   = stats(strcmp(stats_tags,'R'));
                anomRLO(i,j) = stats(strcmp(stats_tags,'RLO')) - anomR(i,j) ;
                anomRUP(i,j) = stats(strcmp(stats_tags,'RUP')) - anomR(i,j);
                
            end
        end
    end
end

clear LDAS_sm LDAS_sm_org LDAS_prcp_org INSITU_sm INSITU_prcp

disp('writing output ...')
save([fout_name,'stats.mat']);
