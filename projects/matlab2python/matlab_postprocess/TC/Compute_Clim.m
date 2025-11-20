clear

data_path = '/discover/nobackup/qliu/merra_land/matlab/SMAP/L3L4ASCAT/';

MOD_version = 'NRv11.4exp_M36';
resolution = 'M36';
SMAPL3_version = 'R18290';
ASCAT_version = 'H119';

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
    f_EASE = '/css/smapl4/public/L4_Products/L4_SM/Vv7030/lmc/SMAP_L4_SM_lmc_00000000T000000_Vv7030_001.h5';
    lon_L4 = double(h5read(f_EASE,'/cell_lon'));
    lat_L4 = double(h5read(f_EASE,'/cell_lat'));
end
M09_Nlon = size(lon_L4,1);
M09_Nlat = size(lon_L4,2);

mod_sm_sum =  zeros(M09_Nlon*M09_Nlat,365);
L3_sm_sum =  zeros(M09_Nlon*M09_Nlat,365);
ASC_sm_sum =  zeros(M09_Nlon*M09_Nlat,365);

N_sm_clim =  zeros(M09_Nlon*M09_Nlat,365);

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
    
    fname = [data_path,'ASCAT_',ASCAT_version(1:4), ...
        '_L3_',SMAPL3_version,'_',MOD_version,'_SMSF_',date_string,'.mat'];
    
    if exist(fname,'file')
        disp(['loading ',fname])
        tmp = load(fname,'sm_L3','sm_mod','sm_ASC','idx_lonxlat');
        
        idx = tmp.idx_lonxlat;
        
        if ~isempty(idx)
            idofyr = dofyr_list(:,this_doy);
            for i = 1:length(idofyr)
                mod_sm_sum(idx,idofyr(i)) = mod_sm_sum(idx,idofyr(i)) + tmp.sm_mod;
                L3_sm_sum(idx,idofyr(i)) = L3_sm_sum(idx,idofyr(i)) + tmp.sm_L3;
                ASC_sm_sum(idx,idofyr(i)) = ASC_sm_sum(idx,idofyr(i)) + tmp.sm_ASC;
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
L3_sm_clim = single(L3_sm_sum ./N_sm_clim);
ASC_sm_clim = single(ASC_sm_sum ./N_sm_clim);
N_sm_clim = int32(N_sm_clim_to_save); clear N_sm_clim_to_save

% save pentad clim to reduce size
mod_sm_clim = mod_sm_clim(:,3:5:365);
L3_sm_clim = L3_sm_clim(:,3:5:365);
ASC_sm_clim = ASC_sm_clim(:,3:5:365);
N_sm_clim = N_sm_clim(:,3:5:365);

fout= [data_path,'ASC_L3_',MOD_version,'_clim_pentad_',time_tag,'_w31.mat'];

save(fout,'-v7.3','mod_sm_clim','L3_sm_clim','ASC_sm_clim','N_sm_clim','Nday_min');
