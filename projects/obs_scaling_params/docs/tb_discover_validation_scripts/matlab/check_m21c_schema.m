addpath('/discover/nobackup/projects/land_da/geosldas-analysis/projects/matlab2python/shared/matlab/');
fname = '/discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_OLv8_M36/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/rc_out/Y2016/M01/LS_OLv8_M36.ldas_obsparam.20160101_0000z.txt';
[N, obs_param] = read_obsparam(fname);
for i=1:N
  fprintf('i=%d descr=%s species=%g orbit=%g pol=%g ang=%g varname=%s units=%s\n', ...
    i, obs_param(i).descr, obs_param(i).species, obs_param(i).orbit, obs_param(i).pol, ...
    obs_param(i).ang(1), obs_param(i).varname, obs_param(i).units);
end
