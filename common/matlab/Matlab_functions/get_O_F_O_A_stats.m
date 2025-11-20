
function [] = get_O_F_O_A_stats(varname,                    ...
    run_months, exp_path, exp_run, domain, start_year, end_year,   ...
    dt_assim, species, hscale, N_data_min, sm_scaling, count_time, convert_grid, time_of_day_in_hours )

% Updated from
% get_model_and_obs_stats.m
%
% Compute mean, var, skewness, min, max, and histogram of O_A
% and innovations from tile-based "innov" files.
%
% The main purpose of this function is to aggregate the information
% from the "innov" files so that it can be plotted quickly afterwards
%
% Stats output file contains tile ID, number of
% data, mean, variance, skewness, min, max, and histogram
% for each catchment (see write statement)
%
% GDL, 11 feb 2011, add correlation as output
% GDL, 04 aug 2011, pass N_data_min as argument
% GDL, 11 oct 2012, allow sm-scaling of forecast values
%
% GDL, aug 2013: added convert_grid (= EASEv2_M36, EASE_M36, ...)
%                M09 obs are administered by tiles (0) that could be anywhere
%                around the center of the observed pixel (M36)
%                -----------
%                | X X X X |
%                | X O O X |
%                | X O O X |
%                | X X X X |
%                -----------
% -------------------------------------------------------------------
% begin user-defined inputs
% -------------------------------------------------------------------

% obs species to be processed (see ens_upd_inputs.nml for a list)
%
% (only observation species that represent observations of the same
%  model prognostic or diagnostic can be processed together!)

%species = [1 2];
%species = [1];
%species = [2];

% edges for histogram

old_output_style = 0;   % switch for old output file and dir conventions

nodata     = -9999;
nodata_tol = 1e-4;

if(strcmp(varname,'''sfmc''')||strcmp(varname,'sfmc'))         % surface soil moisture
    
    % NOTE: edge_dx should be chosen ridiculously small. Then the
    % resulting polynomial fit in get_cdf_match.m will approximate
    % most closely what would be obtained by fitting the ranked
    % scatter diagram.  (The latter cannot be obtained easily for
    % large domains and many observations because it is unlikely
    % to fit into memory, and looping over each tile would require
    % reading the entire data set of "innov" files once for each
    % tile.) - reichle, 22 Aug 2005
    
    edge_min = -0.4;
    edge_max = 0.4;
    edge_dx  = 0.001;    % increase edge_dx if memory (RAM) is a problem
    spec_tag = 'SM';
    
    % coarse histogram for SMMR vs. AMSR-E stats - CANNOT be used
    % with get_cdf_match.m!!!
    % also increase edge_max for stats of scaled data
    % (reichle, 27 Jan 2006)
    
    if 0
        
        edge_min = 0;
        edge_max = 0.6;
        edge_dx  = 0.03;    % increase edge_dx if memory (RAM) is a problem
        
        edge_min = -.18;
        edge_max = 0.72;
        edge_dx  = 0.03;    % increase edge_dx if memory (RAM) is a problem
        
    end
    
elseif(strcmp(varname,'''Tb''')||strcmp(varname,'Tb'))        %Tbh/v
    
    edge_min = -50.0;
    edge_max = 50.0;
    %  edge_dx  = 0.2;
    edge_dx  = 10; %14Oct15 to speed up and limit file size
    spec_tag = 'Tb';
    
elseif strcmp(varname,'Tbn')        %Tbh/v
    
    edge_min = -3.0;
    edge_max = 3.0;
    edge_dx  = 0.05;
    spec_tag = 'Tb';
    
    
else
    
    input('Unknown species requested. Ctrl-C now.')
    
end

for s=1:length(species)
    spec_tag = [spec_tag,'_',num2str(species(s))];
end

sc_tag = {'','_sc'};
sc_tag = sc_tag{sm_scaling+1};

% minimum number of data points to include in statistics

% N_data_min = 1;    % leave decision about "good" stats for later

% no-data-value for points that don't have good statistics

no_data_stats = -9999.;

% -------------------------------------------------------------------
% end user-defined inputs
% -------------------------------------------------------------------

edges = edge_min:edge_dx:edge_max;

N_edges = length(edges);

% -------------------------------------------------------------

% assemble input and output paths

inpath  = [ exp_path, '/', exp_run, '/output/', domain ];

if ~exist(inpath, 'dir')
    inpath  = [ exp_path, '/output/', domain ];
end

outpath = [ inpath, '/stats/'   ];

% create outpath if it doesn't exist

if exist(outpath)~=2
    eval(['!mkdir -p ', outpath]);
end

% -------------------------------------------------------------

% assemble output file name

month_string = {'Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; ...
    'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'};

if strcmp(varname, 'Tb')
    fname_out = [ outpath, '/', exp_run, '.O_F_O_Astats.',                    ...
        'hscale_', num2str(hscale,'%2.2f'), '_',             ...
        num2str(start_year), '-', num2str(end_year), '.' ];
elseif strcmp(varname, 'Tbn')
    fname_out = [ outpath, '/', exp_run, '.O_Fn_O_Astats.',                    ...
        'hscale_', num2str(hscale,'%2.2f'), '_',             ...
        num2str(start_year), '-', num2str(end_year), '.' ];
end

if (iscell(run_months))
    
    if ((end_year-start_year+1)~=length(run_months))
        error('cell entries should match number of years, i.e. need months for each year');
    else
        for i=1:length(run_months) %for all years
            fname_out = [fname_out, month_string{min([run_months{i}])},'_',month_string{max([run_months{i}])},'.' ];
        end
    end
    
else
    
    if length(run_months)==12
        fname_out = [fname_out, 'all_months' ];
    else
        for i=1:length(run_months)
            fname_out = [fname_out, month_string{run_months(i)} ];
        end
    end
    
end


fname_out = [fname_out, spec_tag];

if exist( 'time_of_day_in_hours', 'var')
    
    fname_out = [fname_out, '_', num2str(time_of_day_in_hours,'%2.2d'), 'z'];
    
end

if (exist('convert_grid','var'))
    fname_out = [fname_out,convert_grid];
end

if (count_time)
    fname_out = [fname_out,'_Nt'];
end

fname_out = [fname_out, sc_tag, '.dat'];

% check whether output file exists

if (exist(fname_out)==2)
    
    disp(['output file exists. not overwriting. returning'])
    disp(fname_out)
    return
    
else
    
    disp(['creating ', fname_out])
    
end

% -------------------------------------------------------------

% load catchment coordinates

if old_output_style
    
    %  fname = [inpath, '/', exp_run, '_tile_coord.dat'];
    fname = [inpath, '/rc_out/', exp_run, '.ldas_tilecoord.txt'];
    
else
    
    fname = [inpath, '/rc_out/', exp_run, '.ldas_tilecoord.bin'];
    
end

[ tile_coord ] = read_tilecoord_GEOS( fname );

N_tile = length(tile_coord.tile_id);

% -------------------------------------------------------------

% determine tiles to whose statistics the current obs will contribute

disp('pre-computing index for regional averaging')

central_lat        = tile_coord.com_lat;
central_lon        = tile_coord.com_lon;
tile_coord_tile_id = tile_coord.tile_id;

if (exist('convert_grid'))
    
    %1) convert to M35 EASE indices
    %2) convert back to lat/lon at center of obs
    if (~isempty(strfind(convert_grid, 'M36')) && ~isempty(strfind(convert_grid, 'EASEv2')))
        gridid = 'M36';
        [central_row,central_col] = smapeasev2_latlon2ind(central_lat,central_lon,gridid);
        [central_lat,central_lon] = smapeasev2_ind2latlon(central_row,central_col,gridid);
    elseif (~isempty(strfind(convert_grid, 'M36')) && ~isempty(strfind(convert_grid, 'EASEv1')))
        gridid = 'M36';
        [central_row,central_col] = smapeasev1_latlon2ind(central_lat,central_lon,gridid);
        [central_lat,central_lon] = smapeasev1_ind2latlon(central_row,central_col,gridid);
    else
        error(['Unable to convert to ',convert_grid])
    end
    
    row_col_tmp         = [central_row central_col];
    [unique_rc, ia, ic] = unique(row_col_tmp,'rows');
    N_tile_obs = length(ia);
    
    max_Hx_c = length(find(mode(ic)==ic));
    
    tile_coord_tile_id = zeros(N_tile_obs,max_Hx_c);
    
    disp(['centralizing ',num2str(N_tile_obs),' obs on ',convert_grid,...
        ' grid before doing stats: max ',num2str(max_Hx_c), 'tiles per obs cell'])
    
    for i=1:N_tile_obs
        
        tmp_ind = find(row_col_tmp(:,1) == unique_rc(i,1)  & row_col_tmp(:,2) == unique_rc(i,2));
        tile_coord_tile_id(i,1:length(tmp_ind)) = tile_coord.tile_id(tmp_ind);
        
    end
    
else
    
    max_Hx_c = 1;
    
    N_tile_obs = N_tile;
    ia = 1:N_tile;
    ic = 1:N_tile;
    
end

% -------------------------------------------------------------

% get parameters to rescale sm, if wanted

if (sm_scaling == 1)
    
    catfname = [inpath,'/rc_out/','/Y2010/M01/',...
        exp_run,'.ldas_catparam.','20100101_0000','z.bin'];
    
    mwparname = [inpath,'/rc_out/','/Y2010/M01/',...
        exp_run,'.ldas_mwRTMparam.','20100101_0000','z.bin'];
    
    cat_param = read_catparam( catfname, tile_coord.N_tile );
    
    [ mwRTM_param, tile_check, mwRTMtileid ] = read_mwRTMparam_18_20( mwparname );
    for t=1:tile_check
        if (mwRTMtileid(t) ~= tile_coord.tile_id(t))
            error('tileid of tile_coord.txt-file does not match that in RTMparamfile');
        end
    end
    
end

% -------------------------------------------------------------

% initialize output statistics

N_data    = zeros(N_tile_obs,1);               % same for obs and obsana
N_dat_time= zeros(N_tile_obs,1);
count     = zeros(N_tile_obs,1);

O_A_mean = zeros(N_tile_obs,1);
O_A_var  = zeros(N_tile_obs,1);
O_F_mean = zeros(N_tile_obs,1);
O_F_var  = zeros(N_tile_obs,1);
O_F_skew  = zeros(N_tile_obs,1);
O_A_corr    = zeros(N_tile_obs,1);
O_F_O_A_var= zeros(N_tile_obs,1);
O_ensvar_mean = zeros(N_tile_obs,1);
F_ensvar_mean = zeros(N_tile_obs,1);

O_A_min  =  inf*ones(N_tile_obs,1);
O_A_max  = -inf*ones(N_tile_obs,1);

O_A_hist = zeros(N_tile_obs,N_edges-1);

obs_mean = zeros(N_tile_obs,1);
obs_var  = zeros(N_tile_obs,1);
obsana_mean = zeros(N_tile_obs,1);
obsana_var  = zeros(N_tile_obs,1);
obsana_obs_mean = zeros(N_tile_obs,1);

% tmp debugging

O_F_all = [];

% determine tiles to whose statistics the current obs will contribute

disp('pre-computing index for regional averaging')

if hscale>0
    
    for i=1:N_tile_obs
        this_lat = central_lat(ia(i));
        this_lon = central_lon(ia(i));
        
        tmp_sq_distance =                           ...
            (central_lon - this_lon).^2 +    ...
            (central_lat - this_lat).^2;
        
        hscale_ind{i} = find( tmp_sq_distance <= hscale^2 );
    end
    
else
    
    hscale_ind =num2cell(ia);
    
end


% -------------------------------------------------------------

tmp_months = run_months;

for year = start_year:end_year
    
    if (iscell(tmp_months))
        
        run_months = [tmp_months{year-start_year+1}];
        
    end
    
    for month = run_months
        
        for day = 1:days_in_month( year, month)
            
            for seconds_in_day = 0:dt_assim:(86400-1)
                
                hour    = floor(seconds_in_day/3600);
                
                % check if diurnal stats are needed
                
                if exist('time_of_day_in_hours','var')
                    tmp_hour = time_of_day_in_hours;
                else
                    tmp_hour = hour;       % all hours of day will be included
                end
                
                if hour==tmp_hour
                    
                    minute  = floor( (seconds_in_day-hour*3600)/60 );
                    
                    seconds = seconds_in_day-hour*3600-minute*60;
                    
                    if (seconds~=0)
                        input('something is wrong! Ctrl-c now')
                    end
                    
                    YYYYMMDD = [ num2str(year,   '%4.4d'),     ...
                        num2str(month,  '%2.2d'),     ...
                        num2str(day,    '%2.2d')  ];
                    
                    HHMM     = [ num2str(hour,   '%2.2d'),     ...
                        num2str(minute, '%2.2d')  ];
                    
                    % read innov files
                    
                    fname = [ inpath, '/ana/ens_avg/',                  ...
                        'Y', YYYYMMDD(1:4), '/',                  ...
                        'M', YYYYMMDD(5:6), '/',                  ...
                        exp_run, '.ens_avg.ldas_ObsFcstAna.',        ...
                        YYYYMMDD, '_', HHMM, 'z.bin' ];

                    if ~exist(fname,'file')
                              fname = strrep(fname, 'Vv6032', 'Vv6030');
                              fname = strrep(fname, 'Vv7032', 'Vv7030');
                    end
                    
                    ifp = fopen( fname, 'r', 'b' );
                    
                    if (ifp > 0)           %Proceed only if file exists (e.g. irregular SMOS swaths!)
                        
                        fclose(ifp);
                        
                        [date_time,              ...
                            obs_assim,              ...
                            obs_species,            ...
                            obs_tilenum,            ...
                            obs_lon,                ...
                            obs_lat,                ...
                            obs_obs,                ...
                            obs_obsvar,             ...
                            obs_fcst,               ...
                            obs_fcstvar,            ...
                            obs_ana,                ...
                            obs_anavar              ...
                            ] =                      ...
                            read_ObsFcstAna_GEOS( fname );

                        % the obs in the product OL4001 is unscaled, so the obs 
                        % need to read from the DA product
                        if ~isempty(strfind(exp_run, 'OL4001'))
                           fname_DA = strrep(fname,'OL4001', 'Vv4030');

                           [date_time,              ...
                            obs_assim,              ...
                            obs_species_DA,            ...
                            obs_tilenum_DA,            ...
                            obs_lon,                ...
                            obs_lat,                ...
                            obs_obs,                ...
                            obs_obsvar,             ...
                            obs_fcst_tmp,               ...
                            obs_fcstvar_tmp,            ...
                            obs_ana_tmp,                ...
                            obs_anavar_tmp              ...
                            ] =                      ...
                            read_ObsFcstAna_GEOS( fname_DA );
                        
                           % need to reorder the obs_fcst read from the
                           % open loop file to match the tilenum from the
                           % DA file. non-existing obs_fcst are filled with
                           % NaNs.
                                                 
                               
                           [Lia, Locb] = ismember([obs_tilenum_DA obs_species_DA], ...
                                                  [obs_tilenum obs_species],'rows');
                           
                           % LocB == 0 means no tilenum in DA found in OL
                           tmpidx = find(Locb >= 1);
                           tmp = nodata *ones(size(obs_obs)) ;
                           tmp(tmpidx) = obs_fcst(Locb(tmpidx)); 
                           clear obs_fcst
                           obs_fcst = tmp; clear tmp
                           
                           tmp = nodata*ones(size(obs_obs));
                           tmp(tmpidx) = obs_fcstvar(Locb(tmpidx)); 
                           clear obs_fcstvar
                           obs_fcstvar = tmp; clear tmp
                           
                           tmp = nodata*ones(size(obs_obs));
                           tmp(tmpidx) = obs_ana(Locb(tmpidx)); 
                           clear obs_ana
                           obs_ana = tmp; clear tmp
                           
                           tmp = nodata*ones(size(obs_obs));
                           tmp(tmpidx) = obs_anavar(Locb(tmpidx)); 
                           clear obs_anavar
                           obs_anavar = tmp; clear tmp
                           
                           obs_tilenum = obs_tilenum_DA;     clear obs_tilenum_DA        
                           obs_species = obs_species_DA;     clear obs_species_DA                       

                           
                           clear Locb Lia tmpidx
                           
                           clear obs_fcst_tmp obs_fcstvar_tmp obs_ana_tmp obs_anavar_tmp
                        end

                        
                        % extract species of interest
                        
                        if length(obs_assim) >= 1
                            
                            ind = [];
                            
                            for this_species = species
                                
                                %ind = [ ind;  find( obs_species == this_species) ];
                                
                                %modified by Q. Liu to exclude non-assimilated obs points
                                ind = [ ind;  find( obs_species == this_species & obs_assim==1) ];
                                
                            end
                            
                            obs_species = obs_species(ind);
                            obs_tilenum = obs_tilenum(ind);
                            obs_value   = obs_obs(ind);
                            obs_obsvar   = obs_obsvar(ind);
                            obs_fcst    = obs_fcst(ind);
                            obs_fcstvar    = obs_fcstvar(ind);
                            obs_ana     = obs_ana(ind);
                            
                            if (sm_scaling == 1)
                                obs_fcst = obs_fcst.*mwRTM_param.poros(obs_tilenum)./cat_param.poros(obs_tilenum);
                            end
                            
                            obs_assim   = obs_assim(ind);
                            
                            if strcmp(varname,'Tbn')
                                innov       = (obs_value-obs_fcst)./sqrt(obs_obsvar+obs_fcstvar);
                            else
                                innov       = obs_value-obs_fcst;
                            end
                            
                            O_A         = obs_value - obs_ana;
                            
                            if (count_time)
                                count(:) = 0;
                            end
                            
                            for i=1:length(obs_value)
                                
                                if (  abs(obs_value(i) - nodata) > nodata_tol && ...
                                        abs(obs_ana(i) - nodata) > nodata_tol   &&...
                                        abs(obs_fcst(i)- nodata) > nodata_tol   &&...
                                        abs(O_A(i)  - nodata) > nodata_tol      &&...
                                        abs(innov(i)- nodata) > nodata_tol )
                                    
                                    %ind = hscale_ind{obs_tilenum(i)};
                                    %model tiles
                                    
                                    ind = hscale_ind{ic(obs_tilenum(i))};
                                    %obs tiles
                                    ind = unique(ic(ind));
                                    
                                    % add up data
                                    N_data(ind) = N_data(ind) + 1;  % same for obs and obsana
                                    if (count_time)
                                        if (count(ind) == 0)
                                            N_dat_time(ind) = N_dat_time(ind) + 1;
                                            count(ind) = 1;
                                        end
                                    end
                                    
                                    % O_F, O_A
                                    
                                    O_A_mean(ind) = O_A_mean(ind) + O_A(i);
                                    O_F_mean(ind) = O_F_mean(ind) + innov(i);
                                    
                                    O_A_var(ind)  = O_A_var(ind)  + O_A(i)^2;   %O_A var
                                    O_F_var(ind)  = O_F_var(ind)  + innov(i)^2; %O_F var
                                    O_F_skew(ind)  = O_F_skew(ind)  + innov(i)^3; %O_F skewness 
                                    O_ensvar_mean(ind) = O_ensvar_mean(ind) + obs_obsvar(i);
                                    F_ensvar_mean(ind) = F_ensvar_mean(ind) + obs_fcstvar(i);
                                    
                                    O_F_O_A_var(ind) = O_F_O_A_var(ind) + innov(i).*O_A(i);
                                    
                                    O_A_min(ind)  = min(O_A_min(ind), O_A(i));
                                    O_A_max(ind)  = max(O_A_max(ind), O_A(i));
                                    
                                    %HISTOGRAM
                                    binindex = find( innov(i) > edges );
                                    if (isempty(binindex))
                                        j = 1;          %Storing values below min edge in first bin
                                    elseif (max(binindex) >= N_edges-1)
                                        j = N_edges-1;  %Storing values above max edge in last bin
                                    else
                                        j = max(binindex);
                                    end
                                    
                                    obsana_value = obs_ana(i);
                                    
                                    obs_mean(ind) = obs_mean(ind) + obs_value(i);
                                    obs_var(ind)  = obs_var(ind)  + obs_value(i)^2;
                                    
                                    obsana_mean(ind) = obsana_mean(ind) + obsana_value;
                                    obsana_var(ind)  = obsana_var(ind)  + obsana_value^2;
                                    
                                    obsana_obs_mean(ind) = obsana_obs_mean(ind) + obsana_value*obs_value(i);
                                    
                                    O_A_hist(ind,j) = O_A_hist(ind,j) + 1;
                                    
                                end % condition on valid data
                                
                            end % loop through observations/innovations
                            
                        end % if file present
                    end
                end    % time_of_day_in_hours
                
            end      % seconds_in_day
        end        % day
    end          % month
end            % year

% normalize sums

N_data_all = N_data;
if (count_time)
    N_data = N_dat_time;
end

% pick catchments with at least N_data_min data points

disp(' ')
disp(['total number of catchments = ', num2str(N_tile)])

ind_incl = find(N_data >= N_data_min);
ind_excl = find(N_data <  N_data_min);

disp(['number of catchments with at least ', num2str(N_data_min), ...
    ' data points = ', num2str(length(ind_incl))])

disp(['number of catchments without any data points = ', ...
    num2str(length(find(N_data==0)))])

disp(['number of catchments with 1 to ', num2str(N_data_min-1), ...
    ' data points = ', num2str(length(ind_excl)-length(find(N_data==0)))])


% throw out catchments with less than N_data_min data points

NN = N_data_all(ind_incl);

O_A_mean(ind_excl)     = [];
O_A_var(ind_excl)      = [];
O_F_mean(ind_excl)     = [];
O_F_var(ind_excl)      = [];
O_F_skew(ind_excl)      = [];
O_ensvar_mean(ind_excl)      = [];
F_ensvar_mean(ind_excl)      = [];
O_F_O_A_var(ind_excl)      = [];

% normalize

O_A_mean = O_A_mean./NN;
O_F_mean = O_F_mean./NN;

O_ensvar_mean = O_ensvar_mean./NN;
F_ensvar_mean = F_ensvar_mean./NN;

% NOTE normalization of skewness is *not* from textbook
%      (result differs slightly from matlab function)

O_A_var  = (O_A_var-NN.*O_A_mean.^2)./(NN-1);
O_F_skew  = (O_F_skew./(NN-1) - 3*O_F_var./NN.*O_F_mean+2*O_F_mean.^3);
O_F_var  = (O_F_var-NN.*O_F_mean.^2)./(NN-1);
O_F_skew = O_F_skew./(O_F_var.^(1.5));

O_F_O_A_var = (O_F_O_A_var-NN.*(O_F_mean.*O_A_mean))./(NN-1);

O_A_corr = (obsana_obs_mean-(obsana_mean.*obs_mean))./...
    (sqrt(obs_var).*sqrt(obsana_var));

% expand stats back to arrays for all catchments,
% insert no-data-values for points without statistics

O_A_mean(ind_incl) = O_A_mean;
O_A_mean(ind_excl) = no_data_stats;
O_A_mean(isnan(O_A_mean)) = no_data_stats;

O_A_var(ind_incl) = O_A_var;
O_A_var(ind_excl) = no_data_stats;
O_A_var(isnan(O_A_var)) = no_data_stats;

O_F_O_A_var(ind_incl) = O_F_O_A_var;
O_F_O_A_var(ind_excl) = no_data_stats;
O_F_O_A_var(isnan(O_F_O_A_var)) = no_data_stats;

O_A_corr(ind_incl) = O_A_corr;
O_A_corr(ind_excl) = no_data_stats;
O_A_corr(isnan(O_A_corr)) = no_data_stats;

% set no-data for min/max

O_A_min( isinf(O_A_min)) = no_data_stats;
O_A_max( isinf(O_A_max)) = no_data_stats;

% additional quality control
% (variance of moisture content may be close to zero and become negative
%  due to roundoff error)

ind = find( O_A_var(ind_incl)<0 );
O_A_var(ind_incl(ind)) = no_data_stats;
O_A_corr(ind_incl(ind)) = no_data_stats;

ind = find( obs_var(ind_incl)<0 );
O_A_corr(ind_incl(ind)) = no_data_stats;

ind = find( obsana_var(ind_incl)<0 );
O_A_corr(ind_incl(ind)) = no_data_stats;

ind = find( O_A_corr(ind_incl)<-1);
O_A_corr(ind_incl(ind)) = no_data_stats;

ind = find( O_A_corr(ind_incl)>1);
O_A_corr(ind_incl(ind)) = no_data_stats;

if length(ind)>0
    disp(['found negative obs variance ', ...
        ' for catchments ', num2str(ind_incl(ind)') ])
    disp(['setting respective components of O_A_var to no-data-value'])
end

clear ind

% write output file

disp(' ')
disp(['writing ', fname_out])

ofp = fopen(fname_out, 'w');

format_string = [''];
for i=1:max_Hx_c
    format_string = [ format_string, ' %10d'];
end
format_string = [format_string, '%7d%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%13.3f%13.3f'];

for i=1:N_edges-1
    format_string = [ format_string, ' %d'];
end
format_string = [ format_string, '\n'];

for k=1:N_tile_obs
    
    fprintf(ofp, format_string,                       ...
        tile_coord_tile_id(k,:), N_data(k),  O_F_mean(k), O_F_var(k),     ...
        O_A_mean(k), O_A_var(k), O_F_skew(k),  ...
        F_ensvar_mean(k),  O_ensvar_mean(k), O_A_corr(k),     ...
        edge_min,    edge_max,                    ...
        O_A_hist(k,:) );
    
end


fclose(ofp);

% ================================================================


% -------------------------------------------------------------------

% echo some diagnostics

% figure out what fraction of data was thrown to satisfy N_data_min
%  and other quality control

N_data_elim  = sum( N_data( find(isnan(O_A_mean)) ) );

N_data_total = sum( N_data );

disp(' ')
disp(['N_data_total = ', num2str(N_data_total)])
disp(['N_data_elim  = ', num2str(N_data_elim)])
disp(['fraction of lost data points = ', num2str(N_data_elim/N_data_total)])

%disp([num2str(NN)])
%disp([num2str(N_tile)])

% ==================== EOF ==============================================
