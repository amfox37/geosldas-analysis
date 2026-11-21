#!/usr/bin/env python3
# checks 3-4 of the VIIRS SCF reader test: parse the ldas_log (logit) files of one experiment
# usage: check_jpss_logs.py <exp_dir>
# per analysis cycle: files read by read_JPSS_SCF_hdf / read_MODIS_SCF_hdf, N_good_data, "cannot find file" lines
import collections, glob, os, re, sys

exp_dir = sys.argv[1]
logs = sorted(glob.glob(os.path.join(exp_dir, 'output/*/rc_out/Y*/M*/*.ldas_log.*.txt')))
print('log files:', [os.path.basename(f) for f in logs])

cycles = collections.OrderedDict()      # anal time -> list of (reader, product, fname, N_good_data)
not_found = []
cycle = None
pending = None                          # (reader, product, fname) waiting for its N_good_data line
for fn in logs:
    for line in open(fn, errors='replace'):
        m = re.search(r'get_enkf_increments\(\): enter at anal time (\d{8}_\d{6})z', line)
        if m:
            cycle = m.group(1)
            cycles.setdefault(cycle, [])
            continue
        m = re.search(r'(read_\w+_SCF_hdf)\(\): reading \w+ SCF obs from\s+(\S+)', line)
        if m:
            fname = m.group(2)
            pending = (m.group(1), os.path.basename(fname)[:7], fname)
            continue
        m = re.search(r'N_good_data\s*=\s*(\d+)', line)
        if m and pending is not None:
            cycles[cycle].append(pending + (int(m.group(1)),))
            pending = None
            continue
        if 'cannot find file' in line:
            not_found.append((cycle, line.strip()))

products = sorted({r[1] for c in cycles.values() for r in c})
print('analysis cycles found:', len(cycles), ' products read:', products)
for prod in products:
    print('\n==', prod)
    n_two = collections.Counter()
    for cyc, reads in cycles.items():
        r = [x for x in reads if x[1] == prod]
        if len(r) == 2:
            n_two[cyc[:8]] += 1
        print(' ', cyc, ' n_files=%d' % len(r), ' N_good_data =', [x[3] for x in r],
              ' files =', [os.path.basename(x[2])[8:16] for x in r])
    days = sorted({c[:8] for c in cycles})
    print('  cycles per day with two files:', {d: n_two[d] for d in days})
print('\ncannot find file lines:', len(not_found))
for c, l in not_found[:20]:
    print(' ', c, l)
