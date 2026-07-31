function [ N_obs_param, obs_param ] = read_obsparam_newschema( fname )
%
% TEMPORARY, test-only reader for the current GEOSldas ldas_obsparam schema:
%  - adds fcstvarname/fcstunits after varname/units, relative to the
%    checked-in matlab2python/shared/matlab/read_obsparam.m (last updated
%    2017-06-08 for flistpath/flistname). Mirrors
%    obs_scaling/io.py:read_obs_param field order exactly.
%  - uses a quote-aware tokenizer (tokenize_obsparam_file.m) so quoted
%    fields with embedded spaces (e.g. fcstunits='m3 m-3') are not
%    mis-split, unlike fscanf('%s').
% Does not modify or replace the production reference reader.
%
tokens = tokenize_obsparam_file(fname);
disp(['Reading (new-schema, quote-aware) ',fname]);
idx = 1;
    function tok = next_raw()
        tok = tokens{idx};
        idx = idx + 1;
    end
    function s = next_str()
        tok = next_raw();
        s = tok(2:end-1);
    end
    function v = next_num()
        v = str2double(next_raw());
    end

N_obs_param = round(next_num());
for i=1:N_obs_param
    obs_param(i).descr           = next_str();
    obs_param(i).species         = next_num();
    obs_param(i).orbit           = next_num();
    obs_param(i).pol             = next_num();
    obs_param(i).N_ang           = next_num();
    ang = zeros(1, obs_param(i).N_ang);
    for a=1:obs_param(i).N_ang
        ang(a) = next_num();
    end
    obs_param(i).ang             = ang;
    obs_param(i).freq            = next_num();
    obs_param(i).FOV             = next_num();
    obs_param(i).FOV_units       = next_str();
    obs_param(i).assim           = next_str();
    obs_param(i).scale           = next_str();
    obs_param(i).getinnov        = next_str();
    obs_param(i).RTM_ID          = next_num();
    obs_param(i).bias_Npar       = next_num();
    obs_param(i).bias_trel       = next_num();
    obs_param(i).bias_tcut       = next_num();
    obs_param(i).nodata          = next_num();
    obs_param(i).varname         = next_str();
    obs_param(i).units           = next_str();
    obs_param(i).fcstvarname     = next_str();   % new field
    obs_param(i).fcstunits       = next_str();   % new field
    obs_param(i).path            = next_str();
    obs_param(i).name            = next_str();
    obs_param(i).maskpath        = next_str();
    obs_param(i).maskname        = next_str();
    obs_param(i).scalepath       = next_str();
    obs_param(i).scalename       = next_str();
    obs_param(i).flistpath       = next_str();
    obs_param(i).flistname       = next_str();
    obs_param(i).errstd          = next_num();
    obs_param(i).std_normal_max  = next_num();
    obs_param(i).zeromean        = next_str();
    obs_param(i).coarsen_pert    = next_str();
    obs_param(i).xcorr           = next_num();
    obs_param(i).ycorr           = next_num();
    obs_param(i).adapt           = next_num();
end
disp(['Done reading obs_param for ',num2str(N_obs_param),' species']);
end
