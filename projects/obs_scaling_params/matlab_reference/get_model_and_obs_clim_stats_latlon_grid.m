function [] = get_model_and_obs_clim_stats_latlon_grid( species_names, ...
        run_months, exp_path, exp_run, domain, start_year, end_year,   ...
        dt_assim, t0_assim, species, combine_species_stats,            ...
        resol, w_days, Ndata_min, prefix, print_each_DOY,              ...
        print_each_pentad, print_all_pentads, out_dir, enable_dedup )

% Default for optional switch
if nargin < 20 || isempty(enable_dedup)
    enable_dedup = true;   % ON by default
end

% ---------------- user inputs / constants ----------------
nodata     = -9999;
nodata_tol = 1e-4;
overwrite  = 1;
Nf         = 7;     % [o_mean, o_std, m_mean, m_std, N, m_min, m_max]
N_pentads  = 73;

% Rounding tolerances for duplicate key
LONLAT_Q   = 1e3;   % 0.001° for lon/lat
OBS_Q      = 1e4;   % 1e-4 for obs value

disp('ASSUMING obs are not on Earth-fixed regular grid (e.g., ASCAT)');
disp(['Calculating scaling parameters on grid with resolution = ', num2str(resol) , ' degrees']);

if combine_species_stats
  N_species = 1;
else
  N_species = length(species);
end

inpath  = [ exp_path, '/', exp_run, '/output/', domain ];
outpath = [ inpath, '/stats/', out_dir ];
if ~exist(outpath, 'dir'); mkdir(outpath); end

% Output name bases
ind  = start_year == min(start_year); mi_m = min(run_months(ind));
ind  = end_year   == max(end_year);   ma_m = max(run_months(ind));
D(1) = 1; P(1) = 1;
if mi_m > 1
  D(1) = sum(days_in_month(2014, 1:mi_m-1)) + 1;
  P(1) = ceil(D(1) / 5);
end
D(2) = sum(days_in_month(2014, 1:ma_m));
P(2) = floor(D(2) / 5);

fname_out_base_d = [outpath,  '/', prefix,                                            ...
                    num2str(min(start_year)), '_doy', num2str(D(1)), '_',             ...
                    num2str(max(end_year)), '_doy', num2str(D(2)),                    ...
                    '_W_', num2str(w_days), 'd_Nmin_', num2str(Ndata_min)];
fname_out_base_p = [outpath, '/', prefix,                                             ...
                    num2str(min(start_year)), '_p', num2str(P(1)), '_',               ...
                    num2str(max(end_year)), '_p', num2str(P(2)),                      ...
                    '_W_', num2str(round(w_days/5)), 'p_Nmin_', num2str(Ndata_min)];

% ---------------- lat/lon grid ----------------
ll_lon = -180; ll_lat = -90;
d_lon = resol; d_lat = resol;
n_lon = round(360 / d_lon);
n_lat = round(180 / d_lat);
ll_lons = linspace(ll_lon, ll_lon + (n_lon-1)*d_lon, n_lon);
ll_lats = linspace(ll_lat, ll_lat + (n_lat-1)*d_lat, n_lat);

obsnum         = (1:n_lon*n_lat)';
[i_out, j_out] = ind2sub([n_lon, n_lat], obsnum);
lon_out        = ll_lons(i_out)'; %#ok<NASGU>
lat_out        = ll_lats(j_out)'; %#ok<NASGU>
N_gridcells    = length(obsnum);

% ---------------- accumulators ----------------
o_data_sum  = NaN(N_species, N_gridcells, w_days);
m_data_sum  = NaN(N_species, N_gridcells, w_days);
o_data_sum2 = NaN(N_species, N_gridcells, w_days);
m_data_sum2 = NaN(N_species, N_gridcells, w_days);
m_data_min  = NaN(N_species, N_gridcells, w_days);
m_data_max  = NaN(N_species, N_gridcells, w_days);
N_data      = NaN(N_species, N_gridcells, w_days);

data_out    = NaN(N_species, Nf, N_gridcells, N_pentads);
data2D      = NaN(Nf, N_gridcells);

% For pentad outputs (filled lazily)
start_time_p = struct([]); %#ok<NASGU>
end_time_p   = struct([]); %#ok<NASGU>

% ---------------- de-dup state (optional, persists across cycles) ----------------
if enable_dedup
    dedup     = containers.Map('KeyType','char','ValueType','any');
    build_key = @(tile, sp, lon, lat, obs, yr) sprintf('%d_%d_%d_%d_%d_%d', ...
                   tile, sp, round(lon*LONLAT_Q), round(lat*LONLAT_Q), round(obs*OBS_Q), yr);
    cycle_id       = 0;
    KEEP_N_CYCLES  = 3;   % keep keys for last ~3 cycles
end

% ---------------- time loops ----------------
t0_assim = mod(t0_assim, dt_assim);
count = 0;

for imonth = 1:length(run_months)
  month = run_months(imonth);

  for day = 1:days_in_month(2014, month)
    if count < w_days, count = count + 1; else, count = w_days; end

    for seconds_in_day = t0_assim:dt_assim:(86400-1)
      hour    = floor( seconds_in_day/3600);
      minute  = floor((seconds_in_day-hour*3600)/60);
      seconds = seconds_in_day-hour*3600-minute*60;
      if (seconds ~= 0); input('something is wrong! Ctrl-c now'); end

      % New cycle: bump counter and evict old dedup keys
      if enable_dedup
        cycle_id = cycle_id + 1;
        if ~isempty(dedup)
          ks = keys(dedup);
          to_rm = false(size(ks));
          for ii = 1:numel(ks)
            v = dedup(ks{ii});
            if v.last_cycle_id < cycle_id - KEEP_N_CYCLES
              to_rm(ii) = true;
            end
          end
          if any(to_rm), remove(dedup, ks(to_rm)); end
        end
      end

      for year = start_year(imonth):end_year(imonth)
        YYYYMMDD = [num2str(year, '%4.4d'), num2str(month, '%2.2d'), num2str(day, '%2.2d')];
        HHMM     = [num2str(hour, '%2.2d'), num2str(minute, '%2.2d')];

        % ObsFcstAna file
        fname = [inpath, '/ana/ens_avg/', 'Y', YYYYMMDD(1:4), '/', 'M', YYYYMMDD(5:6), '/', ...
                 exp_run, '.ens_avg.ldas_ObsFcstAna.', YYYYMMDD, '_', HHMM, 'z.bin'];
        ifp = fopen(fname, 'r', 'l');

        if (ifp > 0)
          fclose(ifp);
          [date_time, obs_assim, obs_species, obs_tilenum, obs_lon, obs_lat, ...
           obs_obs, obs_obsvar, obs_fcst, obs_fcstvar, obs_ana, obs_anavar] = read_ObsFcstAna(fname); %#ok<ASGLU>

          % Drop rows with no forecast
          idx = isnan(obs_fcst);
          obs_assim(  idx) = [];
          obs_species(idx) = [];
          obs_tilenum(idx) = [];
          obs_lon(    idx) = [];
          obs_lat(    idx) = [];
          obs_obs(    idx) = [];
          obs_obsvar( idx) = [];
          obs_fcst(   idx) = [];
          obs_fcstvar(idx) = [];
          obs_ana(    idx) = [];
          obs_anavar( idx) = [];
          date_time(  idx) = [];

          % Species selection
          for scnt = 1:N_species
            if combine_species_stats
              ind = find(ismember(obs_species, species));
            else
              this_species = species(scnt);
              ind = find(obs_species == this_species);
            end

            if ~isempty(ind)
              obs_tilenum_i = obs_tilenum(ind);
              obs_species_i = obs_species(ind);
              obs_obs_i     = obs_obs(    ind);
              obs_fcst_i    = obs_fcst(   ind);
              obs_lon_i     = obs_lon(    ind);
              obs_lat_i     = obs_lat(    ind);

              % Map to grid
              i_idx = floor((obs_lon_i - ll_lon) / d_lon) + 1;
              j_idx = floor((obs_lat_i - ll_lat) / d_lat) + 1;
              [~, obs_idx] = ismember([i_idx, j_idx], [i_out, j_out], 'rows');
              valid = obs_idx > 0;
              obs_idx       = obs_idx(valid);
              obs_i         = obsnum(obs_idx);
              obs_tilenum_i = obs_tilenum_i(valid);
              obs_species_i = obs_species_i(valid);
              obs_obs_i     = obs_obs_i(valid);
              obs_fcst_i    = obs_fcst_i(valid);
              obs_lon_i     = obs_lon_i(valid);
              obs_lat_i     = obs_lat_i(valid);

              % Accumulate (with optional de-dup keeping the later occurrence)
              for r = 1:numel(obs_i)
                oi   = obs_i(r);
                tile = obs_tilenum_i(r);
                spid = obs_species_i(r);
                lonr = obs_lon_i(r);
                latr = obs_lat_i(r);
                oval_raw = obs_obs_i(r);
                mval     = double(obs_fcst_i(r));

                if isnan(oval_raw)
                  nval  = 0; oval = 0; o2val = 0;
                else
                  nval  = 1; oval = double(oval_raw); o2val = oval^2;
                end
                m2val = mval^2;

                if enable_dedup
                  key = build_key(tile, spid, lonr, latr, double(oval_raw), year);

                  % If seen recently, remove earlier contribution (keep later)
                  if isKey(dedup, key)
                    prev = dedup(key);
                    val = o_data_sum(prev.scnt, prev.oi, prev.count);  if isnan(val), val = 0; end
                    o_data_sum(prev.scnt, prev.oi, prev.count)  = val - prev.o;
                    val = m_data_sum(prev.scnt, prev.oi, prev.count);  if isnan(val), val = 0; end
                    m_data_sum(prev.scnt, prev.oi, prev.count)  = val - prev.m;
                    val = o_data_sum2(prev.scnt, prev.oi, prev.count); if isnan(val), val = 0; end
                    o_data_sum2(prev.scnt, prev.oi, prev.count) = val - prev.o2;
                    val = m_data_sum2(prev.scnt, prev.oi, prev.count); if isnan(val), val = 0; end
                    m_data_sum2(prev.scnt, prev.oi, prev.count) = val - prev.m2;
                    val = N_data(prev.scnt, prev.oi, prev.count);      if isnan(val), val = 0; end
                    N_data(prev.scnt, prev.oi, prev.count)      = val - prev.N;
                  end
                end

                % Add current contribution
                val = o_data_sum(scnt, oi, count);  if isnan(val), val = 0; end
                o_data_sum(scnt, oi, count)  = val + oval;
                val = m_data_sum(scnt, oi, count);  if isnan(val), val = 0; end
                m_data_sum(scnt, oi, count)  = val + mval;
                val = o_data_sum2(scnt, oi, count); if isnan(val), val = 0; end
                o_data_sum2(scnt, oi, count) = val + o2val;
                val = m_data_sum2(scnt, oi, count); if isnan(val), val = 0; end
                m_data_sum2(scnt, oi, count) = val + m2val;
                val = N_data(scnt, oi, count);      if isnan(val), val = 0; end
                N_data(scnt, oi, count)      = val + nval;

                % Update min/max with later forecast value (approximate)
                curr_min = m_data_min(scnt, oi, count);
                curr_max = m_data_max(scnt, oi, count);
                if isnan(curr_min), m_data_min(scnt, oi, count) = mval; else, m_data_min(scnt, oi, count) = min(curr_min, mval); end
                if isnan(curr_max), m_data_max(scnt, oi, count) = mval; else, m_data_max(scnt, oi, count) = max(curr_max, mval); end

                % Remember latest kept copy for this key
                if enable_dedup
                  dedup(key) = struct('scnt', scnt, 'oi', oi, 'count', count, ...
                                      'o', oval, 'm', mval, 'o2', o2val, 'm2', m2val, 'N', nval, ...
                                      'last_cycle_id', cycle_id);
                end
              end
            end
          end
        end
      end
    end

    % Write daily DOY output once the initial window is filled
    if count >= w_days
      end_time.year  = 2014; end_time.month = month; end_time.day = day;
      end_time.hour  = hour;  end_time.min   = minute; end_time.sec = seconds;
      start_time     = augment_date_time( -floor(w_days*(24*60*60)), end_time );

      % Clean nodata sentinels
      o_data_sum(abs(o_data_sum - nodata) <= nodata_tol) = NaN;
      m_data_sum(abs(m_data_sum - nodata) <= nodata_tol) = NaN;

      for i = 1:N_species
        N_window    =      sum(N_data(     i,:,1:w_days),   3,"omitnan");
        data2D(1,:) =      sum(o_data_sum( i,:,1:w_days),   3,"omitnan")./N_window;
        data2D(2,:) = sqrt(sum(o_data_sum2(i,:,1:w_days),   3,"omitnan")./N_window - data2D(1,:).^2);
        data2D(3,:) =      sum(m_data_sum( i,:,1:w_days),   3,"omitnan")./N_window;
        data2D(4,:) = sqrt(sum(m_data_sum2(i,:,1:w_days),   3,"omitnan")./N_window - data2D(3,:).^2);
        data2D(5,:) = N_window;
        data2D(6,:) =      min(m_data_min( i,:,1:w_days),[],3);
        data2D(7,:) =      max(m_data_max( i,:,1:w_days),[],3);

        data2D([1:Nf], N_window < Ndata_min) = NaN;

        DOY = augment_date_time(-floor(w_days*(24*60*60)/2.0), end_time).dofyr;
        if(is_leap_year(end_time.year) && DOY>=59)
          error('This code should never hit a leap year');
        end

        if print_each_DOY
          pentad = floor((DOY + 2)/5);
          if combine_species_stats
            fname_out = [fname_out_base_d, '_sp_ALL_DOY', num2str(DOY,'%3.3d'), '.nc4'];
          else
            fname_out = [fname_out_base_d,'_sp_', char(species_names(i)),'_DOY', num2str(DOY,'%3.3d'), '.nc4'];
          end
          if (exist(fname_out)==2 && overwrite)
            disp(['Output file exists. overwriting ', fname_out])
          elseif (exist(fname_out)==2 && ~overwrite)
            disp(['Output file exists. not overwriting. returning'])
            disp(['Writing ', fname_out])
            return
          else
            disp(['Creating ', fname_out])
          end

          write_netcdf_latlon_grid( fname_out, i_out, j_out, ll_lons, ll_lats, data2D, pentad, ...
                                    start_time, end_time, overwrite, Nf, ll_lon, ll_lat, d_lon, d_lat)
        end

        if mod((DOY + 2),5) == 0
          pentad                 = (DOY + 2)/5;
          data_out(i,:,:,pentad) = data2D;
          start_time_p(pentad)   = start_time; %#ok<AGROW>
          end_time_p(pentad)     = end_time;   %#ok<AGROW>
          if print_each_pentad
            if combine_species_stats
              fname_out = [fname_out_base_p, '_sp_ALL_p', num2str(pentad,'%2.2d'), '.nc4'];
            else
              fname_out = [fname_out_base_p, '_sp_', char(species_names(i)),'_p', num2str(pentad,'%2.2d'), '.nc4'];
            end
            if (exist(fname_out)==2 && overwrite)
              disp(['Output file exists. overwriting ', fname_out])
            elseif (exist(fname_out)==2 && ~overwrite)
              disp(['Output file exists. not overwriting. returning'])
              disp(['Writing ', fname_out])
              return
            else
              disp(['Creating ', fname_out])
            end

            write_netcdf_latlon_grid( fname_out, i_out, j_out, ll_lons, ll_lats, data2D, pentad, ...
                                      start_time, end_time, overwrite, Nf, ll_lon, ll_lat, d_lon, d_lat )
          end
        end

        % Roll the window
        o_data_sum( i,:,1:w_days-1) = o_data_sum( i,:,2:w_days);
        m_data_sum( i,:,1:w_days-1) = m_data_sum( i,:,2:w_days);
        o_data_sum2(i,:,1:w_days-1) = o_data_sum2(i,:,2:w_days);
        m_data_sum2(i,:,1:w_days-1) = m_data_sum2(i,:,2:w_days);
        m_data_min( i,:,1:w_days-1) = m_data_min( i,:,2:w_days);
        m_data_max( i,:,1:w_days-1) = m_data_max( i,:,2:w_days);
        N_data(     i,:,1:w_days-1) = N_data(     i,:,2:w_days);
        o_data_sum( i,:,w_days    ) = NaN;
        m_data_sum( i,:,w_days    ) = NaN;
        o_data_sum2(i,:,w_days    ) = NaN;
        m_data_sum2(i,:,w_days    ) = NaN;
        m_data_min( i,:,w_days    ) = NaN;
        m_data_max( i,:,w_days    ) = NaN;
        N_data(     i,:,w_days    ) = NaN;

        data2D = NaN+0.0.*data2D;
      end
    end
  end
end

% Optional: dump all pentads at end
if print_all_pentads
  for i = 1:N_species
    data_o = squeeze(data_out(i,:,:,:));
    if combine_species_stats
      fname_out = [fname_out_base_d, '_sp_ALL_all_pentads.nc4'];
    else
      fname_out = [fname_out_base_d,'_sp_', char(species_names(i)),'_all_pentads.nc4'];
    end
    if (exist(fname_out)==2 && overwrite)
      disp(['Output file exists. overwriting ', fname_out])
    elseif (exist(fname_out)==2 && ~overwrite)
      disp(['Output file exists. not overwriting. returning'])
      disp(['Writing ', fname_out])
      return
    else
      disp(['Creating ', fname_out])
    end

    write_netcdf_latlon_grid( fname_out, i_out, j_out, ll_lons, ll_lats, data_o, 1:73, ...
                              start_time_p, end_time_p, overwrite, Nf, ll_lon, ll_lat, d_lon, d_lat )
  end
end

% ==================== EOF ==============================================