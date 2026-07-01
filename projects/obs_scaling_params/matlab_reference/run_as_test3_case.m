% Run MATLAB reference scaling for LS_DAv8_M36_as_test3 (all 14 species,
% May 2020, 3-hourly assim, 3-day window, 1-deg grid).
% Output goes to scratch; input data read through symlinks.

clear

repo_root = getenv('OBS_SCALING_REPO_ROOT');
if isempty(repo_root)
    repo_root = char(java.io.File(pwd).getCanonicalPath());
end

scratch_root = getenv('OBS_SCALING_SCRATCH_ROOT');
if isempty(scratch_root)
    scratch_root = fullfile(getenv('HOME'), 'obs_scaling_as_test3');
end

addpath(fullfile(repo_root, 'projects', 'obs_scaling_params', 'matlab_reference'));
addpath(fullfile(repo_root, 'projects', 'matlab2python', 'shared', 'matlab'));

exp_path = scratch_root;
exp_run  = 'LS_DAv8_M36_as_test3';
domain   = 'SMAP_EASEv2_M36_GLOBAL';

species_names = {'SMOS_fit_Tbh_A', 'SMOS_fit_Tbh_D', ...
                 'SMOS_fit_Tbv_A', 'SMOS_fit_Tbv_D', ...
                 'SMAP_L1C_Tbh_A', 'SMAP_L1C_Tbh_D', ...
                 'SMAP_L1C_Tbv_A', 'SMAP_L1C_Tbv_D', ...
                 'ASCAT_META_SM',  'ASCAT_METB_SM',   'ASCAT_METC_SM', ...
                 'MYD10C1',        'MOD10C1',          'CYGNSS_SM_6hr'};
species              = 1:14;
combine_species_stats = 1;

run_months    = 5;
start_year    = 2020;
end_year      = 2020;
dt_assim      = 3 * 60 * 60;
t0_assim      = 0;
grid_resolution = 1.0;
w_days        = 3;
Ndata_min     = 1;
prefix_out    = 'as_test3_matlab_zscore_stats_';
print_each_DOY    = 1;
print_each_pentad = 0;
print_all_pentads = 0;
out_dir       = 'matlab_as_test3_zscore_1deg';
enable_dedup  = true;

get_model_and_obs_clim_stats_latlon_grid( ...
    species_names, run_months, exp_path, exp_run, domain, start_year, ...
    end_year, dt_assim, t0_assim, species, combine_species_stats, ...
    grid_resolution, w_days, Ndata_min, prefix_out, print_each_DOY, ...
    print_each_pentad, print_all_pentads, out_dir, enable_dedup);

disp('MATLAB as_test3 run complete');
