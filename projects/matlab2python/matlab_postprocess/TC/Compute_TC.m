clear

data_path = '/discover/nobackup/qliu/merra_land/matlab/SMAP/L3L4ASCAT/';

MOD_version = 'NRv11.4exp_M36';
resolution = 'M36';
SMAPL3_version = 'R18290';
ASCAT_version = 'H119';

do_year = 1;
do_summer = 0;

Nmin = 20 ;

start_time.year = 2015;
start_time.month = 4;
start_time.day = 1;
start_time.hour = 0;
start_time.min = 0;
start_time.sec = 0;

end_time.year = 2023;
end_time.month = 4;
end_time.day = 1;
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

% load climatology file

f_clim= [data_path,'ASC_L3_',MOD_version,'_clim_pentad_',time_tag,'_w31.mat'];

disp(f_clim)

tmp = load(f_clim,'mod_sm_clim','L3_sm_clim','ASC_sm_clim');
mod_clim = tmp.mod_sm_clim; %[N_tile, N_pentad]
L3_clim = tmp.L3_sm_clim;
ASC_clim = tmp.ASC_sm_clim/200.;
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
    
    fname = [data_path,'ASCAT_',ASCAT_version(1:4), ...
        '_L3_',SMAPL3_version,'_',MOD_version,'_SMSF_',date_string,'.mat'];
    
    if exist(fname,'file')
        
        disp(fname)
        
        tmp = load(fname,'sm_mod','sm_L3','sm_ASC','idx_lonxlat');
        tmp.sm_ASC = tmp.sm_ASC/200.;
        
        idx = tmp.idx_lonxlat;
        
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
            sm_L3 = tmp.sm_L3 - tmp_clim(idx); clear tmp_clim
            
            tmp_clim = mod_clim(:, date_time.pentad);
            sm_mod = tmp.sm_mod - tmp_clim(idx); clear tmp_clim
            
            tmp_clim = ASC_clim(:, date_time.pentad);
            sm_ASC = tmp.sm_ASC - tmp_clim(idx); clear tmp_clim
            
            % exclude point where no climatology values
            iv = find(~isnan(sm_mod));
            % convert compressed data back on global
            mod_sm_sum(idx(iv)) = mod_sm_sum(idx(iv)) + sm_mod(iv);
            mod_sm2_sum(idx(iv)) = mod_sm2_sum(idx(iv)) + sm_mod(iv).^2;
            
            L3_sm_sum(idx(iv)) = L3_sm_sum(idx(iv)) + sm_L3(iv);
            L3_sm2_sum(idx(iv)) = L3_sm2_sum(idx(iv)) + sm_L3(iv).^2;
            
            ASC_sm_sum(idx(iv)) = ASC_sm_sum(idx(iv)) + sm_ASC(iv);
            ASC_sm2_sum(idx(iv)) = ASC_sm2_sum(idx(iv)) + sm_ASC(iv).^2;
            
            L3xASC_sm_sum(idx(iv)) = L3xASC_sm_sum(idx(iv)) + sm_L3(iv) .* sm_ASC(iv);
            L3xmod_sm_sum(idx(iv)) = L3xmod_sm_sum(idx(iv)) + sm_L3(iv) .* sm_mod(iv);
            modxASC_sm_sum(idx(iv)) = modxASC_sm_sum(idx(iv)) + sm_mod(iv) .* sm_ASC(iv);
            
            N_sm(idx(iv)) = N_sm(idx(iv)) + 1;
            
            clear sm_L3 sm_mod sm_ASC tmp iv
        end
        
        clear idx fname
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
C_mod_ASC(C_mod_ASC < 0.) = NaN;
C_L3_mod = L3xmod_sm_sum ./NN  - L3_sm_mean .* mod_sm_mean;    
C_L3_mod(C_L3_mod < 0.) = NaN;
C_L3_ASC = L3xASC_sm_sum ./NN  - L3_sm_mean .* ASC_sm_mean; 
C_L3_ASC(C_L3_ASC < 0.) = NaN;

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

fname_out = [data_path,'ASC_L3_',MOD_version,'_TC_stats_',time_tag];

if ~do_year
    if do_summer
        fname_out = [fname_out,'_summer'];
    else
        fname_out = [fname_out,'_winter'];
    end
end

fname_out = [fname_out,'.mat'];

save(fname_out,'N_sm','Nmin','R2_TC_L3','R2_TC_ASC','R2_TC_mod','sigma2_L3','sigma2_mod',...
    'sigma2_ASC','R_mod_L3','R_mod_ASC','R_ASC_L3')

