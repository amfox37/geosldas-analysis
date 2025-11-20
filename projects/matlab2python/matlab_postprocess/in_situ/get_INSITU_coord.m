function [INSITU_lat, INSITU_lon, INSITU_id] = get_INSITU_coord(INSITU_tag);

dirbase = '/home/qliu/smap/SMAP_Nature/';%'/Volumes/smap-tb.jpl.nasa.gov';%

addpath([dirbase, '/validation/tools'])

calvalglobalvar(dirbase);

if strcmp(INSITU_tag, 'SCAN')
   % read SCAN coord file
% ============================================================================================
SCAN_coord_file = '/discover/nobackup/qliu/merra_land/DATA/SCAN/data/coord/ALL_SCAN_coord.txt';

fid = fopen(SCAN_coord_file);
C = textscan(fid,'%s%d%f%f%d%d');
fclose(fid);

tmp.SCAN_id = C{2};
tmp.lat = C{3};
tmp.lon = C{4};

clear C

% exclude non-conus sites
ind1 = union(find(tmp.lat > 45), find(tmp.lat <=25));
ind2 = union(find(tmp.lon >= -60), find(tmp.lon < -130));
ind = union(ind1, ind2);
tmp.lat(ind) = [];
tmp.lon(ind) = [];
tmp.SCAN_id(ind) = [];

clear ind

SCAN_lat = tmp.lat;
SCAN_lon = tmp.lon;
for i=1:length(tmp.SCAN_id)
    SCAN_id{i} = num2str(tmp.SCAN_id(i));
end

clear tmp


INSITU_lat = SCAN_lat;
INSITU_lon = SCAN_lon;
INSITU_id  = SCAN_id;

return

elseif strcmp(INSITU_tag, 'USCRN')
% read USCRN coord file
% ==============================================================================================
USCRN_path = '/discover/nobackup/qliu/merra_land/DATA/USCRN/data/'; 
USCRN_id_fname = [USCRN_path, 'list_of_USCRNsites.txt'];

fid = fopen(USCRN_id_fname);
tmp = textscan(fid,'%s');
fclose(fid);

sites = tmp{1};
clear tmp


USCRN_finfo = [USCRN_path, 'USCRN_Station_info.txt'];

fid = fopen(USCRN_finfo);
C = textscan(fid,'%s%f%f%s%s','headerlines',1);
fclose(fid);

tmp.USCRN_id = C{1};
tmp.lat = C{2};
tmp.lon = C{3};

clear C

for i = 1:length(sites)
    s_ind(i) = find(strcmp(tmp.USCRN_id, sites{i}));
    USCRN_id{i} = tmp.USCRN_id{s_ind(i)};
    USCRN_lat(i)= tmp.lat(s_ind(i));
    USCRN_lon(i)= tmp.lon(s_ind(i));
end

clear tmp

INSITU_lat = USCRN_lat;
INSITU_lon = USCRN_lon;
INSITU_id  = USCRN_id;

return
elseif strcmp(INSITU_tag, 'Oznet')

% read Oznet coord info
% =================================================================================================
Oznet_id_fname = [ '/discover/nobackup/qliu/merra_land/DATA/Murrumbidgee/QLiu_201904/data/coord/list_MURRUM_sites.txt'];

fid = fopen(Oznet_id_fname);
tmp = textscan(fid,'%s');
fclose(fid);

sites = tmp{1};
clear tmp


Oznet_finfo = [ '/discover/nobackup/qliu/merra_land/DATA/Murrumbidgee/QLiu_201904/data/coord/MURRUM_coord.txt'];

fid = fopen(Oznet_finfo);
C = textscan(fid,'%s%s%f%f%f%f');
fclose(fid);

tmp.Oznet_id = C{2};
tmp.lat = C{3};
tmp.lon = C{4};

clear C

for i = 1:length(sites)
    s_ind(i) = find(strcmp(tmp.Oznet_id, sites{i}));
    Oznet_id{i} = tmp.Oznet_id{s_ind(i)};
    Oznet_lat(i)= tmp.lat(s_ind(i));
    Oznet_lon(i)= tmp.lon(s_ind(i));
end

clear tmp

INSITU_lat = Oznet_lat;
INSITU_lon = Oznet_lon;
INSITU_id  = lower(Oznet_id);

return
elseif strcmp(INSITU_tag, 'Msnet')

% read Mesonet coord info
% ======================================================================================================
Station_file = '/discover/nobackup/qliu/merra_land/DATA/Oklahoma_Mesonet/geoinfo/Mesonet_geoinfo.csv';
Sinfo = readtable(Station_file);
sites = Sinfo.stnm;

Msnet_lat = Sinfo.nlat;
Msnet_lon = Sinfo.elon;
for i = 1:length(Sinfo.stnm)
   Msnet_id{i} = num2str(Sinfo.stnm(i));
end

clear Sinfo

INSITU_lat = Msnet_lat;
INSITU_lon = Msnet_lon;
INSITU_id  = Msnet_id;

return

elseif strcmp(INSITU_tag, 'SMOSM')

% read Mesonet coord info
% ======================================================================================================
Station_file = '/discover/nobackup/qliu/merra_land/DATA/ISMN/data/coord/list_SMOSMANIA_sites.txt';

fid = fopen(Station_file);
tmp = textscan(fid,'%s');
fclose(fid);

sites = tmp{1};
clear tmp

SMOSM_finfo = [ '/discover/nobackup/qliu/merra_land/DATA/ISMN/data/coord/coord_SMOSMANIA_sites.txt'];

fid = fopen(SMOSM_finfo);
C = textscan(fid,'%s%s%f%f%f%f');
fclose(fid);

tmp.SMOSM_id = C{2};
tmp.lat = C{3};
tmp.lon = C{4};

clear C

for i = 1:length(sites)
    s_ind(i) = find(strcmp(tmp.SMOSM_id, sites{i}));
    SMOSM_id{i} = tmp.SMOSM_id{s_ind(i)};
    SMOSM_lat(i)= tmp.lat(s_ind(i));
    SMOSM_lon(i)= tmp.lon(s_ind(i));
end

clear tmp

INSITU_lat = SMOSM_lat;
INSITU_lon = SMOSM_lon;
INSITU_id  = SMOSM_id;

return

elseif strcmp(INSITU_tag, 'CalVal_M33')

% read CalVal coord info
% ======================================================================================================

refpixlist = { '03013302', ... % '03010903', '03010908', ...
           '04013302', ... % '04010907', '04010910', ...
           '07013301', ... %'07010902', '07010916', ...
           '09013301', ... %'09010906', ... 
           '12033301', ...
           '16013302', ... %'16010906', '16010907', '16010913', ...
           '16023302', ... %'16020907', ...
           '16033302', ... %'16030911', '16030916', ...
           '16043302', ... %'16040901', ...
           '16063302', ...
           '16073302', ... %'16070909', '16070911', ...
           '19023301', ... %'19020902', ...
           '25013301', ... %'25010911', ...
           '27013301', ... %'27010910', '27010911', ...
           ...%'41010906',  ...  
           '45013301', ... %'45010902', ...
           '45023301', ... %'45020902', ...            
           '48013301', ... %'48010902',  '48010911', ...
           ...%'53013301', ... %
           '67013301'};

[col, row, scale, active, site, refpixii, description, refpixids, rootzone, sitealias, lat,lon] = get_refpix_master_info(refpixlist, 'SM', 'active');

INSITU_lat = lat; 
INSITU_lon = lon; 
INSITU_id = refpixids;
 
return

elseif strcmp(INSITU_tag, 'CalVal_M09')

% read CalVal coord info
% ======================================================================================================

refpixlist = { '03010903', '03010908', ...
           '04010907', '04010910', ...
           '07010902', '07010916', ...
           '09010906', ... 
           '16010906', '16010907', '16010913', ...
           '16020905', '16020906','16020907',...
           '16030911', '16030916', ...
           '16040901', '16040906',...
           '16060907', ...
           '16070909', '16070910','16070911', ...
           '19020902', ...
           '25010911', ...
           '27010910', '27010911', ...
           '41010906',  ...  
           '45010902', ...
           '45020902', ...            
           '48010902',  '48010911', ...
           '67010901' ...
           };

[col, row, scale, active, site, refpixii, description, refpixids, rootzone, sitealias, lat,lon] = get_refpix_master_info(refpixlist, 'SM', 'active');

INSITU_lat = lat; 
INSITU_lon = lon; 
INSITU_id = refpixids;
 
return

else

   error(['no coord information available for ', INSITU_tag]);

end


