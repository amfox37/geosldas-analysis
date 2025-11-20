clear all
addpath /discover/nobackup/qliu/gdelanno_RTM/MATLAB_LDASSA

data_path = '/discover/nobackup/projects/land_da/Evaluation/IVs/output/';

mod_version = 'OLv8_M36_cd'; 

mod_resolution = 'M36'; 

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

% get EASEv2 coord
if strcmp(mod_resolution,'M36')
[lat,lon] = smapeasev2_ind2latlon([0:405],[0:963],'M36');
lon_L4 = repmat(lon',[1,length(lat)]);
lat_L4 = repmat(lat,[length(lon),1]);
clear lat lon
else
%f_EASE = '/css/smapl4/public/L4_Products/L4_SM/Vv5030/lmc/SMAP_L4_SM_lmc_00000000T000000_Vv5030_001.h5';
f_EASE = '/home/qliu/smap/SMAP_L4/L4_SM/Vv5030/lmc/SMAP_L4_SM_lmc_00000000T000000_Vv5030_001.h5';
lon_L4 = double(h5read(f_EASE,'/cell_lon'));
lat_L4 = double(h5read(f_EASE,'/cell_lat'));
end

Nlon = size(lon_L4,1);
Nlat = size(lon_L4,2);
mod_sm_sum =  zeros(Nlon*Nlat,365);
obs_sm_sum =  zeros(Nlon*Nlat,365);
N_sm_clim =  zeros(Nlon*Nlat,365);

% minimum number requirement for computing the climatology of  day
Nday_min = 4 * (end_time.year - start_time.year);

clear dofyr_list
for dofyr = 1:365
    if dofyr <= 15
        dofyr_list(:,dofyr) = [1:dofyr+15 365-(15-dofyr):365];
    elseif dofyr >= 351
        dofyr_list(:,dofyr) = [dofyr-15:365 1:15-(365-dofyr)];
    else
        dofyr_list(:,dofyr) = [dofyr-15 : dofyr+15];
    end
end

date_time = start_time;
date_time = get_dofyr_pentad(date_time);

while 1
    
    if (date_time.year ==end_time.year   && ...
            date_time.month==end_time.month  && ...
            date_time.day  ==end_time.day    )
        break
        
    end
    
    % leap year 2/29 will contribute to the climatology of 2/28 of a normal
    % year
    
    % get the dofyr in a normal year of this date_time 
    tmp_date_time = date_time;
    % any non-leap year
    tmp_date_time.year = 2017;
    
    if tmp_date_time.month == 2 && tmp_date_time.day == 29
        tmp_date_time.day = 28;
    end
    tmp_date_time = get_dofyr_pentad(tmp_date_time);
    this_doy = tmp_date_time.dofyr;
    clear tmp_date_time
    
    date_string = [num2str(date_time.year,'%4.4d'), ...
        num2str(date_time.month,'%2.2d'),num2str(date_time.day,'%2.2d')]
    
        fname = [data_path,'ASCL4_H119_SMSF_L4_',mod_version,'_QC_1_',date_string,'.mat'];
        
        if exist(fname,'file')
            disp(fname)
            tmp = load(fname,'sm_mod','sm_obs','idx_EASEv2_lonxlat');
            
            idx = tmp.idx_EASEv2_lonxlat; 
            
            if ~isempty(idx)
                idofyr = dofyr_list(:,this_doy);
                for i = 1:length(idofyr)
                    mod_sm_sum(idx,idofyr(i)) = mod_sm_sum(idx,idofyr(i)) + tmp.sm_mod;
                    obs_sm_sum(idx,idofyr(i)) = obs_sm_sum(idx,idofyr(i)) + tmp.sm_obs;
                    N_sm_clim(idx,idofyr(i)) = N_sm_clim(idx,idofyr(i)) + 1;
                end
            end
            
            clear idx tmp
            
        end
    
    date_time = augment_date_time(86400, date_time);
    
end

% compute the climatology
N_sm_clim_to_save =  N_sm_clim;
N_sm_clim(N_sm_clim < Nday_min) = NaN;
mod_sm_clim = single(mod_sm_sum ./N_sm_clim); 
obs_sm_clim = single(obs_sm_sum ./N_sm_clim); 
N_sm_clim = int32(N_sm_clim_to_save); clear N_sm_clim_to_save

% save pentad clim to reduce size
mod_sm_clim = mod_sm_clim(:,3:5:365);
obs_sm_clim = obs_sm_clim(:,3:5:365);
N_sm_clim = N_sm_clim(:,3:5:365);

% if end_time.month == 1
% fout= [data_path,'SMPL3_',mod_version,'_clim_pentad_',num2str(start_time.year), ...
%     num2str(start_time.month,'%2.2d'),'_', num2str(end_time.year-1), ...
%     '12_w31',orbit,'.mat'];
% else
% fout= [data_path,'SMPL3_',mod_version,'_clim_pentad_',num2str(start_time.year), ...
%     num2str(start_time.month,'%2.2d'),'_', num2str(end_time.year), ...
%     num2str(end_time.month-1,'%2.2d'), '_w31.mat'];
% end

if end_time.month == 1
fout= [data_path,'ASCL4_',mod_version,'_clim_pentad_',num2str(start_time.year), ...
    num2str(start_time.month,'%2.2d'),'_', num2str(end_time.year-1), ...
    '12_w31',orbit,'.mat'];
else
fout= [data_path,'ASCL4_',mod_version,'_clim_pentad_',num2str(start_time.year), ...
    num2str(start_time.month,'%2.2d'),'_', num2str(end_time.year), ...
    num2str(end_time.month-1,'%2.2d'), '_w31.mat'];
end

save(fout,'-v7.3','mod_sm_clim','obs_sm_clim','N_sm_clim','Nday_min');
