#!/usr/bin/env python3
# download daily snow cover CMG files from NSIDC, renamed for the GEOSldas readers
# usage: get_snow_cmg_obs.py <outdir> <first_day YYYY-MM-DD> <last_day YYYY-MM-DD> <product> [...]
import json, os, subprocess, sys, urllib.request
versions = {'VNP10C1': '2', 'VJ110C1': '2', 'VJ210C1': '2', 'MYD10C1': '61', 'MOD10C1': '61'}
outdir, first_day, last_day = sys.argv[1], sys.argv[2], sys.argv[3]
cookies = os.path.join(outdir, '.edl_cookies')
os.makedirs(outdir, exist_ok=True)
for product in sys.argv[4:]:
    url = ('https://cmr.earthdata.nasa.gov/search/granules.json?short_name=' + product +
           '&version=' + versions[product] +
           '&temporal=' + first_day + 'T00:00:00Z,' + last_day + 'T23:59:59Z&page_size=2000')
    entries = json.load(urllib.request.urlopen(url))['feed']['entry']
    print(product, ':', len(entries), 'granules found')
    for entry in entries:
        parts = entry['title'].split('.')          # product, AYYYYDDD, version, production time, ext
        fname = '.'.join([parts[0], parts[1], parts[2], parts[-1]])
        data_url = [l['href'] for l in entry['links'] if 'data#' in l.get('rel', '')
                    and l['href'].startswith('https') and l['href'].endswith(parts[-1])][0]
        fpath = os.path.join(outdir, product, parts[1][1:5], fname)
        if os.path.exists(fpath):
            continue
        os.makedirs(os.path.dirname(fpath), exist_ok=True)
        subprocess.run(['curl', '-sS', '-f', '-n', '-L', '-c', cookies, '-b', cookies, '-o', fpath, data_url], check=True)
        print('  ', entry['title'], '->', fpath)
