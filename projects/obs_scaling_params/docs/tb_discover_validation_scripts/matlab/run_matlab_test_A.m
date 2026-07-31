% Temporary MATLAB parity-check wrapper for the SMAP Tb tile-space port.
% Calls the UNMODIFIED production reference get_model_and_obs_clim_stats.m
% and UNMODIFIED matlab2python/shared/matlab/read_obsparam.m against a real
% GEOS-LDAS OL experiment (.bin ObsFcstAna archive), descending orbit.

clear

addpath('/discover/nobackup/projects/land_da/geosldas-analysis/projects/obs_scaling_params/matlab_reference');
addpath('/discover/nobackup/projects/land_da/geosldas-analysis/projects/matlab2python/shared/matlab');

exp_path = '/discover/nobackup/projects/land_da/obs_scaling_test';
exp_run  = {'LS_OLv8_M36'};
domain   = 'SMAP_EASEv2_M36_GLOBAL';

run_months = [1 2 3];
start_year = [2016 2016 2016];
end_year   = [2016 2016 2016];

orbit    = [ 1 ]; % 1=A, 2=D
pol      = [ 1 2 ]; % 1=H, 2=V
inc_ang  = [ 40.0 ];

prefix_out = 'MATLAB_TEST_SMAP_Tb_ASC_';

dt_assim = 3*60*60;
t0_assim = 0;

obs_param_fname = [exp_path, '/', exp_run{1}, '/output/', domain, '/rc_out/', ...
    '/Y2016/M01/', exp_run{1}, '.ldas_obsparam.20160101_0000z.txt'];

var_name = {'Tb'};
descr = 'SMAP_L1C';

hscale = 0.0;
w_days = 75;
Ndata_min = 20;
convert_grid = 'EASEv2_M36';

[N_obs_param, obs_param] = read_obsparam(obs_param_fname);

species = [];
for oo=1:length(orbit)
    for pp=1:length(pol)
        for aa=1:length(inc_ang)
            add_species = obs_param(strcmp(var_name,{obs_param.varname}) & ...
                orbit(oo) == [obs_param.orbit] & ...
                inc_ang(aa) == [obs_param.ang] & ...
                pol(pp) == [obs_param.pol] & ...
                ~cellfun(@isempty, strfind({obs_param.descr},descr))).species;
            species = union(species,add_species);
        end
    end
end

disp('Selected species:');
disp(species);

if (orbit(1) == 1) int_Asc = 1; end
if (orbit(1) == 2) int_Asc = 0; end

tic;
get_model_and_obs_clim_stats( var_name, ...
    run_months, exp_path, exp_run{1}, domain, start_year, end_year, ...
    dt_assim, t0_assim, species, obs_param, ...
    hscale, inc_ang, int_Asc, w_days, Ndata_min, prefix_out, ...
    convert_grid );
elapsed = toc;
fprintf('MATLAB TEST RUN (descending) COMPLETE, elapsed=%.1f s\n', elapsed);
