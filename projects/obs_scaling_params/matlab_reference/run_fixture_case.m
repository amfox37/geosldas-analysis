% Run a one-timestamp MATLAB reference case for Python parity checks.
%
% Expected input layout is created by scripts/prepare_fixture_case.py.

clear

repo_root = getenv('OBS_SCALING_REPO_ROOT');
if isempty(repo_root)
    repo_root = char(java.io.File(pwd).getCanonicalPath());
end

fixture_root = getenv('OBS_SCALING_FIXTURE_ROOT');
if isempty(fixture_root)
    fixture_root = '/tmp/obs_scaling_fixture';
end

addpath(fullfile(repo_root, 'projects', 'obs_scaling_params', 'matlab_reference'));
addpath(fullfile(repo_root, 'projects', 'matlab2python', 'shared', 'matlab'));

exp_path = fixture_root;
exp_run = 'LS_DAv8_M36_as_test2';
domain = 'SMAP_EASEv2_M36_GLOBAL';

species_names = {'ASCAT_META_SM', 'ASCAT_METB_SM', 'ASCAT_METC_SM'};
species = [9, 10, 11];
combine_species_stats = 1;

run_months = 5;
start_year = 2020;
end_year = 2020;
dt_assim = 24 * 60 * 60;
t0_assim = 9 * 60 * 60;
grid_resolution = 1.0;
w_days = 1;
Ndata_min = 1;
prefix_out = 'fixture_matlab_zscore_stats_';
print_each_DOY = 1;
print_each_pentad = 0;
print_all_pentads = 0;
out_dir = 'matlab_fixture_zscore_1deg';
enable_dedup = true;

get_model_and_obs_clim_stats_latlon_grid( ...
    species_names, run_months, exp_path, exp_run, domain, start_year, ...
    end_year, dt_assim, t0_assim, species, combine_species_stats, ...
    grid_resolution, w_days, Ndata_min, prefix_out, print_each_DOY, ...
    print_each_pentad, print_all_pentads, out_dir, enable_dedup);

disp('MATLAB fixture run complete');
