clear

addpath('/discover/nobackup/projects/land_da/ldas_code/GEOSldas/src/Components/@GEOSldas_GridComp/GEOSldas_App/util/shared/matlab/')
% addpath('/gpfsm/dnb34/amfox/GEOSldas_diagnostics/Jupyter/Matlab_functions')
% addpath /discover/nobackup/qliu/gdelanno_RTM/MATLAB_LDASSA

data_path = '/discover/nobackup/projects/land_da/Evaluation/IVs/output/';

% choose 2 data for comparison
D1_version = 'OLv8_M36_cd';
D2_version = 'DAv8_M36_cd'; 

% preprocessed skill files
f_D1 =  [data_path,'ASCL4_',D1_version,'_IVD_IVS_stats_lag2day_201808_202405.mat'];
f_D2 =  [data_path,'ASCL4_',D2_version,'_IVD_IVS_stats_lag2day_201808_202405.mat'];

% load skill values, note it's R^2
D1 = load(f_D1,'R2_ivs_obs','R2_ivs_mod');
D2 = load(f_D2,'R2_ivs_obs','R2_ivs_mod');

% get EASE grid coord info
[lat,lon] = EASEv2_ind2latlon([0:405],[0:963],'M36');
lon_EASE = repmat(lon',[1,length(lat)]);
lat_EASE = repmat(lat,[length(lon),1]);

% Compute R and exclude points where estimated R's are too poo
R_D1 = sqrt(D1.R2_ivs_mod);  R_D1(R_D1 < 0.1) = NaN;
R_D2 = sqrt(D2.R2_ivs_mod);  R_D2(R_D2 < 0.1) = NaN;
R_OBS = sqrt(D2.R2_ivs_obs); R_OBS(R_OBS < 0.1) = NaN;

R_D1(isnan(R_OBS)) = NaN;
R_D2(isnan(R_OBS)) = NaN;

% Usually plot the differences only as the R values tend to be overestimated with IVs method
% Also, data are stored in 1-D and need to reshape in 2-D
Rdiff_D2minusD1 = reshape(R_D2 - R_D1, 964, 406);

% Plot with own preferences
% figure, pcolor(lon_EASE,lat_EASE, Rdiff_D2minusD1); shading('flat') 

lons = reshape(lon_EASE, [], 1);
lats = reshape(lat_EASE, [], 1);
Rdiff_vector = R_D2 - R_D1;

% Save the data for later use by Python
outname = ['Rdiff_',D2_version,'_minus_',D1_version,'_ASCL4.mat'];
save(outname,'Rdiff_vector','lats','lons');
