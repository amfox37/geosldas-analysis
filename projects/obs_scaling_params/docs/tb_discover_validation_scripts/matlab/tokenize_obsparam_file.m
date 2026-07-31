function tokens = tokenize_obsparam_file(fname)
%
% TEMPORARY, test-only quote-aware tokenizer for ldas_obsparam files.
% MATLAB's fscanf('%s') is not quote-aware and mis-splits any quoted
% field containing an embedded space (e.g. fcstunits='m3 m-3' for ASCAT
% soil-moisture species), unlike Python's shlex-based reader in
% obs_scaling/io.py. This mirrors shlex whitespace_split=True with quote
% handling. Does not modify the production reference reader.
%
txt = fileread(fname);
n = length(txt);
tokens = {};
i = 1;
while i <= n
    c = txt(i);
    if isspace(c)
        i = i + 1;
        continue;
    end
    if c == ''''
        j = i + 1;
        while j <= n && txt(j) ~= ''''
            j = j + 1;
        end
        tokens{end+1} = txt(i:j); %#ok<AGROW>
        i = j + 1;
    else
        j = i;
        while j <= n && ~isspace(txt(j))
            j = j + 1;
        end
        tokens{end+1} = txt(i:j-1); %#ok<AGROW>
        i = j;
    end
end
end
