clear all

addpath /home/qliu/projects/matlab_code/SMAP/

exp_path = '/discover/nobackup/qliu/merra_land/temp/SMAP_runs/SMAP_Nature_v03/output/';

exp_run = 'SMAP_Nature_v03';

domain = 'SMAP_EASEv2_M09_GLOBAL';

fout_path = '/discover/nobackup/qliu/temp/';

fname_tilecoord = [exp_path, exp_run,'/', domain,'/rc_out/', ...
    exp_run, '.ldas_tilecoord.bin'];

tile_coord = read_tilecoord(fname_tilecoord);
tc = tile_coord;


fname_tilegrids = [exp_path, exp_run,'/', domain,'/rc_out/', ...
    exp_run, '.ldas_tilegrids.bin'];

[tile_grid_g, tile_grid_d] = read_tilegrids(fname_tilegrids);
tg = tile_grid_g;

if ~isempty(strfind(tile_grid_g.gridtype, 'EASE'))
    i_ind = [tg.ind_base:tg.N_lon+(tg.ind_base-1)];
    j_ind = [tg.ind_base:tg.N_lat+(tg.ind_base-1)];
    %if tg.i_dir == -1, i_ind = [i_ind(end):-1:i_ind(1)]; end
    %if tg.j_dir == -1, j_ind = [j_ind(end):-1:j_ind(1)]; end
    
    [lat, lon] = smapeasev2_ind2latlon(j_ind, i_ind, tg.gridtype(end-2:end));
    
else
    
    lon = [(tg.ll_lon+tg.dlon/2.):tg.dlon:(tg.ur_lon-tg.dlon/2.)];
    lat = [(tg.ll_lat+tg.dlat/2.):tg.dlat:(tg.ur_lat-tg.dlat/2.)];
    
end



for ifile = 1:2
    
    if ifile == 1, file_out = 'std_force_pert.nc4'; end
    if ifile == 2, file_out = 'std_progn_pert.nc4'; end
    
    fname_out = [fout_path, file_out];
    
    netcdf.setDefaultFormat('FORMAT_NETCDF4');
    fout_id = netcdf.create(fname_out, 'NETCDF4');
    
    if fout_id < 0, error(['Creating ' fname_out 'failed']); end
    
    % Setup global attributes
    NC_GLOBAL = netcdf.getConstant('GLOBAL');
    
    if ifile == 1
        netcdf.putAtt(fout_id, NC_GLOBAL, 'Title', 'LDAS standard deviation for perturbations of forcing variables');        
    else
        netcdf.putAtt(fout_id, NC_GLOBAL, 'Title', 'LDAS standard deviation for perturbations of prognostics');
    end
    netcdf.putAtt(fout_id, NC_GLOBAL, 'Filename', file_out);
    netcdf.putAtt(fout_id, NC_GLOBAL, 'Institution', 'NASA GMAO');
    netcdf.putAtt(fout_id, NC_GLOBAL, 'History', 'File written by matlab-R2012a');
    netcdf.putAtt(fout_id, NC_GLOBAL, 'Contact', 'NASA/GMAO Rolf Reichle/Qing Liu');
    netcdf.putAtt(fout_id, NC_GLOBAL, 'Comments', 'NETCDF-4');
    
    % Define dimensions:
    dimidY = netcdf.defDim(fout_id,'lat',tg.N_lat);
    dimidX = netcdf.defDim(fout_id,'lon',tg.N_lon);
    
    % Define axis:
     
    varid = netcdf.defVar(fout_id,'lon','double',dimidX);
    netcdf.putAtt(fout_id,varid,'standard_name','longitude');
    netcdf.putAtt(fout_id,varid,'long_name','longitude');
    netcdf.putAtt(fout_id,varid,'units','degrees_east');
    netcdf.putVar(fout_id,varid,lon);
    
    varid = netcdf.defVar(fout_id,'lat','double',dimidY);
    netcdf.putAtt(fout_id,varid,'standard_name','latitude');
    netcdf.putAtt(fout_id,varid,'long_name','latitude');
    netcdf.putAtt(fout_id,varid,'units','degrees_north');
    netcdf.putVar(fout_id,varid,lat);
    
    % Synchronize global 
    netcdf.sync(fout_id)
    
    
    % Define groups and variables in each group
    
    tg_group_id = netcdf.defGrp(fout_id, 'tilegrids');
    netcdf.putAtt(tg_group_id, NC_GLOBAL, 'gridtype', 'EASEv2_M09');   
 
    if ifile ==1
        
        da_group_id = netcdf.defGrp(fout_id, 'std_force_pert');
        vars = {'pcp', 'sw', 'lw', 'tmp2m', 'dpt2m', 'wnd'};
        units = {'mm', 'W m-2', 'W m-2', 'K', 'm', 'm s-1'};
        tmp_value = [0.5 0.3 50 0.1 0.1 0.1];
    else
        
        da_group_id = netcdf.defGrp(fout_id, 'std_progn_pert');
        vars = {'catdef', 'rzexc', 'srfexc', 'snow', 'tc', ...
            'ght1', 'ght2', 'ght3', 'ght4', 'ght5', 'ght6'};
        units = {'mm', 'mm', 'mm', 'mm', 'K', 'W m-2', 'W m-2', ...
            'W m-2', 'W m-2','W m-2', 'W m-2'};
        tmp_value = [0.05 0.1 0.02 0.02 0.2 50000 100000 20000 400000 800000 5000000];
    end
    
    fillValue = single(1.e15);
    
    %% Put data in the tilegrid group:
    
    varid = netcdf.defVar(tg_group_id, 'ind_base', 'int', []);
    netcdf.putVar(tg_group_id, varid,tg.ind_base);
    
    varid = netcdf.defVar(tg_group_id, 'i_dir', 'int', []);
    netcdf.putVar(tg_group_id, varid,tg.i_dir);
    
    varid = netcdf.defVar(tg_group_id, 'j_dir', 'int', []);
    netcdf.putVar(tg_group_id, varid,tg.j_dir);
    
    varid = netcdf.defVar(tg_group_id, 'N_lon', 'int', []);
    netcdf.putVar(tg_group_id, varid,tg.N_lon);
    
    varid = netcdf.defVar(tg_group_id, 'N_lat', 'int', []);
    netcdf.putVar(tg_group_id, varid,tg.N_lat);
    
    varid = netcdf.defVar(tg_group_id, 'i_offg', 'int',[]);
    netcdf.putVar(tg_group_id, varid,tg.i_offg);
    
    varid = netcdf.defVar(tg_group_id, 'j_offg', 'int',[]);
    netcdf.putVar(tg_group_id, varid,tg.j_offg);
    
    varid = netcdf.defVar(tg_group_id, 'll_lon', 'float',[]);
    netcdf.putVar(tg_group_id, varid,tg.ll_lon);
    
    varid = netcdf.defVar(tg_group_id, 'll_lat', 'float',[]);
    netcdf.putVar(tg_group_id, varid,tg.ll_lat);
    
    varid = netcdf.defVar(tg_group_id, 'ur_lon', 'float', []);
    netcdf.putVar(tg_group_id, varid,tg.ur_lon);
    
    varid = netcdf.defVar(tg_group_id, 'ur_lat', 'float',[]);
    netcdf.putVar(tg_group_id, varid,tg.ur_lat);
    
    varid = netcdf.defVar(tg_group_id, 'dlon', 'float', []);
    netcdf.defVarFill(tg_group_id,varid,false,fillValue);
    if isempty(strfind(tg.gridtype, 'EASE'))
        netcdf.putVar(tg_group_id, varid,tg.dlon);
    else
        netcdf.putVar(tg_group_id, varid,fillValue);
    end
    
    varid = netcdf.defVar(tg_group_id, 'dlat', 'float', []);
    netcdf.defVarFill(tg_group_id,varid,false,fillValue);
    if isempty(strfind(tg.gridtype, 'EASE'))
        netcdf.putVar(tg_group_id, varid,tg.dlat);
    else
        netcdf.putVar(tg_group_id, varid,fillValue);
    end
    
    % Synchronize:
    netcdf.sync(tg_group_id)
    
    % Put data into the data group:
    
    
    % Insert data:
    
    for i=1:length(vars)
        
        var_tmp = tmp_value(i);
        tmp = repmat(var_tmp, tc.N_tile, 1);
        tile_data = tmp .* (tc.com_lat + 90)/90.;
        
        data = tile2grid( tile_data', tc, tg);
        data(data<-9998) = fillValue;
        data(isnan(data)) = fillValue;
        
        varid = netcdf.defVar(da_group_id,vars{i},'float',[dimidX dimidY]);
        netcdf.putAtt(da_group_id,varid,'name', ['std_error_pert_',vars{i}]);
        netcdf.putAtt(da_group_id,varid,'units',units{i});
        netcdf.defVarFill(da_group_id,varid,false,fillValue);
        
        
        netcdf.putVar(da_group_id,varid,data');
        clear data
        
    end
    
    % Synchronize:
    netcdf.sync(da_group_id)
    
    netcdf.close(fout_id)
    
end




