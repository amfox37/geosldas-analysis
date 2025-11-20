% Modified triple collocation analysis that reads in data from two separate observation files
% as produced by the IVs/scripts/step2 and step3 scripts
% AMF 2025-10-09

clear

addpath /discover/nobackup/qliu/gdelanno_RTM/MATLAB_LDASSA

data_path = '/discover/nobackup/projects/land_da/Evaluation/IVs/output/';

MOD_version = 'OLv8_M36_cd';
resolution = 'M36';

do_year = 1;
do_summer = 0;

Nmin = 20 ;

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

if end_time.month == 1
    time_tag = [num2str(start_time.year),num2str(start_time.month,'%2.2d'),'_',...
        num2str(end_time.year-1),'12'];
else
    time_tag = [num2str(start_time.year),num2str(start_time.month,'%2.2d'),'_',...
        num2str(end_time.year),num2str(end_time.month-1,'%2.2d')];
end

% get EASEv2 coord
if strcmp(resolution,'M36')
    [lat,lon] = smapeasev2_ind2latlon([0:405],[0:963],'M36');
    lon_L4 = repmat(lon',[1,length(lat)]);
    lat_L4 = repmat(lat,[length(lon),1]);
    clear lat lon
else
    f_EASE = '/css/smapl4/public/L4_Products/L4_SM/Vv5030/lmc/SMAP_L4_SM_lmc_00000000T000000_Vv5030_001.h5';
    lon_L4 = double(h5read(f_EASE,'/cell_lon'));
    lat_L4 = double(h5read(f_EASE,'/cell_lat'));
end

M09_Nlon = size(lon_L4,1);
M09_Nlat = size(lon_L4,2);

L3_sm_sum =  zeros(M09_Nlon*M09_Nlat,1);
L3_sm2_sum =  zeros(M09_Nlon*M09_Nlat,1);

mod_sm_sum =  zeros(M09_Nlon*M09_Nlat,1);
mod_sm2_sum =  zeros(M09_Nlon*M09_Nlat,1);

ASC_sm_sum =  zeros(M09_Nlon*M09_Nlat,1);
ASC_sm2_sum =  zeros(M09_Nlon*M09_Nlat,1);

modxASC_sm_sum =  zeros(M09_Nlon*M09_Nlat,1);
L3xASC_sm_sum =  zeros(M09_Nlon*M09_Nlat,1);
L3xmod_sm_sum =  zeros(M09_Nlon*M09_Nlat,1);

N_sm =  zeros(M09_Nlon*M09_Nlat,1);

% load ASCATL4 climatology file

f_clim= [data_path,'ASCL4_',MOD_version,'_clim_pentad_',time_tag,'_w31.mat'];
disp(f_clim)
tmp = load(f_clim,'mod_sm_clim','obs_sm_clim');
mod_clim = tmp.mod_sm_clim; %[N_tile, N_pentad]
ASC_clim = tmp.obs_sm_clim/200.;
clear tmp

% load SMAPL3 climatology file
f_clim= [data_path,'SMPL3_',MOD_version,'_clim_pentad_',time_tag,'_w31.mat'];
disp(f_clim)
tmp = load(f_clim,'mod_sm_clim','obs_sm_clim');
L3_clim = tmp.obs_sm_clim;
clear tmp

date_time = start_time;
date_time = get_dofyr_pentad(date_time);
id = 0;

while 1
    
    if (date_time.year ==end_time.year   && ...
            date_time.month==end_time.month  && ...
            date_time.day  >=end_time.day    )
        break
        
    end
    
    date_string = [num2str(date_time.year,'%4.4d'), ...
        num2str(date_time.month,'%2.2d'),num2str(date_time.day,'%2.2d')];
    
    fname_asc = [data_path,'ASCL4_H119_SMSF_L4_',MOD_version,'_QC_1_',date_string,'.mat'];
    fname_smp = [data_path,'SMPL3_R19_SMSF_MOD_',MOD_version,'_QC_1_',date_string,'.mat'];
    
    if (exist(fname_asc,'file') && exist(fname_smp,'file'))
        
        disp([fname_asc, ' ',fname_smp])
        
        A = load(fname_asc,'sm_mod','sm_obs','idx_EASEv2_lonxlat');
        S = load(fname_smp,'sm_mod','sm_obs','idx_EASEv2_lonxlat');

        % Alias to avoid collisions and clarify roles
        idx_smap  = double(S.idx_EASEv2_lonxlat(:));
        idx_asc   = double(A.idx_EASEv2_lonxlat(:));

        mod_smap  = S.sm_mod(:);    % model from SMAP file
        obs_smap  = S.sm_obs(:);    % SMAP observation

        mod_asc   = A.sm_mod(:);    % model from ASCAT file
        obs_asc   = A.sm_obs(:);    % ASCAT observation

        % Intersect tiles available in BOTH files today and align order
        [idx_common, ia, ib] = intersect(idx_smap, idx_asc);
        if isempty(idx_common)
            % Nothing to do for this day
            date_time = augment_date_time(86400, date_time);
            continue
        end

        % Extract aligned vectors
        sm_mod_S   = mod_smap(ia);
        sm_mod_A   = mod_asc(ib);
        sm_SMAP    = obs_smap(ia);
        sm_ASCAT   = obs_asc(ib);

        % Prefer a single model field on the overlap (they should be identical)
        % If tiny diffs, just take SMAP-side model
        % (Optional sanity)
        % if any(isfinite(sm_mod_S) & isfinite(sm_mod_A))
        %     assert(max(abs(sm_mod_S - sm_mod_A),[],'omitnan') < 1e-6);
        % end
        sm_mod = sm_mod_S;

        % ASCAT unit normalization (only if your TC used ÷200; remove if already normalized in IV files)
        sm_ASCAT = sm_ASCAT / 200;
        
        idx = idx_common;
        
        if ~do_year
            if  do_summer
                if date_time.month <= 4 || date_time.month >=10
                    idx = [];
                end
            else
                if date_time.month >=5 &&  date_time.month <= 9
                    idx = [];
                end
            end
        end
        
        if ~isempty(idx)
            % compute anomaly
            tmp_clim = L3_clim(:, date_time.pentad);
            sm_L3 = sm_SMAP - tmp_clim(idx); clear tmp_clim
            
            tmp_clim = mod_clim(:, date_time.pentad);
            sm_mod = sm_mod - tmp_clim(idx); clear tmp_clim
            
            tmp_clim = ASC_clim(:, date_time.pentad);
            sm_ASC = sm_ASCAT - tmp_clim(idx); clear tmp_clim
            
            % require ALL THREE finite after anomaly subtraction
            iv = isfinite(sm_mod) & isfinite(sm_L3) & isfinite(sm_ASC);
            if any(iv)
                ii = idx(iv);   % full-grid indices for accumulation

                mod_sm_sum(ii)    = mod_sm_sum(ii)    + sm_mod(iv);
                mod_sm2_sum(ii)   = mod_sm2_sum(ii)   + sm_mod(iv).^2;

                L3_sm_sum(ii)     = L3_sm_sum(ii)     + sm_L3(iv);
                L3_sm2_sum(ii)    = L3_sm2_sum(ii)    + sm_L3(iv).^2;

                ASC_sm_sum(ii)    = ASC_sm_sum(ii)    + sm_ASC(iv);
                ASC_sm2_sum(ii)   = ASC_sm2_sum(ii)   + sm_ASC(iv).^2;

                L3xASC_sm_sum(ii) = L3xASC_sm_sum(ii) + sm_L3(iv).*sm_ASC(iv);
                L3xmod_sm_sum(ii) = L3xmod_sm_sum(ii) + sm_L3(iv).*sm_mod(iv);
                modxASC_sm_sum(ii)= modxASC_sm_sum(ii)+ sm_mod(iv).*sm_ASC(iv);

                N_sm(ii) = N_sm(ii) + 1;
            end

            clear sm_L3 sm_mod sm_ASC tmp iv
        end
        
        clear idx fname_asc fname_smp A S sm_mod_S sm_mod_A sm_SMAP sm_ASCAT idx_common ia ib 
    else
        error(['no file found ',fname])
    end
    
    date_time = augment_date_time(86400, date_time);
    
end

NN = N_sm;
NN(NN < Nmin) = NaN;

L3_sm_mean = L3_sm_sum ./NN;
mod_sm_mean = mod_sm_sum ./NN;
ASC_sm_mean = ASC_sm_sum ./NN;

L3_sm_mean(isinf(L3_sm_mean)) = NaN;
mod_sm_mean(isinf(mod_sm_mean)) = NaN;
ASC_sm_mean(isinf(ASC_sm_mean)) = NaN;

C_L3_L3 = L3_sm2_sum ./NN - L3_sm_mean .^2; 
C_L3_L3(C_L3_L3 < 0.) = NaN;
C_mod_mod = mod_sm2_sum ./NN - mod_sm_mean .^2; 
C_mod_mod(C_mod_mod < 0.) = NaN;
C_ASC_ASC = ASC_sm2_sum ./NN - ASC_sm_mean .^2; 
C_ASC_ASC(C_ASC_ASC < 0.) = NaN;

C_mod_ASC = modxASC_sm_sum ./NN  - mod_sm_mean .* ASC_sm_mean; 
% C_mod_ASC(C_mod_ASC < 0.) = NaN; 
C_L3_mod = L3xmod_sm_sum ./NN  - L3_sm_mean .* mod_sm_mean;    
% C_L3_mod(C_L3_mod < 0.) = NaN;
C_L3_ASC = L3xASC_sm_sum ./NN  - L3_sm_mean .* ASC_sm_mean; 
% C_L3_ASC(C_L3_ASC < 0.) = NaN;

R_mod_L3 = C_L3_mod ./sqrt(C_mod_mod .* C_L3_L3);
R_mod_ASC = C_mod_ASC ./sqrt(C_mod_mod .* C_ASC_ASC);
R_ASC_L3 = C_L3_ASC ./sqrt(C_L3_L3 .* C_ASC_ASC);

R2_TC_L3 = C_L3_mod .* C_L3_ASC ./ C_mod_ASC ./C_L3_L3; %R2_TC_L2(R2_TC_L2 <0 |R2_TC_L2>1) = NaN;
R2_TC_L3(R2_TC_L3 < 0.001) = NaN; R2_TC_L3(R2_TC_L3 > 1.) = 1.;
R2_TC_mod = C_L3_mod .* C_mod_ASC ./ C_L3_ASC ./C_mod_mod; %R2_TC_mod(R2_TC_mod <0 |R2_TC_mod>1) = NaN;
R2_TC_mod(R2_TC_mod < 0.001) = NaN; R2_TC_mod(R2_TC_mod > 1.) = 1.;
R2_TC_ASC = C_mod_ASC .* C_L3_ASC ./ C_L3_mod ./C_ASC_ASC; %R2_TC_ASC(R2_TC_ASC <0 |R2_TC_ASC>1) = NaN;
R2_TC_ASC(R2_TC_ASC < 0.001) = NaN; R2_TC_ASC(R2_TC_ASC > 1.) = 1.;

sigma2_L3 = C_L3_L3 - C_L3_ASC.*C_L3_mod./C_mod_ASC;
sigma2_mod = C_mod_mod - C_mod_ASC.*C_L3_mod./C_L3_ASC;
sigma2_ASC = C_ASC_ASC - C_L3_ASC.*C_mod_ASC./C_L3_mod;

lons = reshape(lon_L4, [], 1);
lats = reshape(lat_L4, [], 1);

fname_out = [data_path,'ASCL4_SMPL3_',MOD_version,'_TC_stats_',time_tag];

if ~do_year
    if do_summer
        fname_out = [fname_out,'_summer'];
    else
        fname_out = [fname_out,'_winter'];
    end
end

fname_out = [fname_out,'.mat'];

save(fname_out,'lons','lats','N_sm','Nmin','R2_TC_L3','R2_TC_ASC','R2_TC_mod','sigma2_L3','sigma2_mod',...
    'sigma2_ASC','R_mod_L3','R_mod_ASC','R_ASC_L3','C_L3_mod','C_mod_ASC','C_L3_ASC')
