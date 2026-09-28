#!/usr/bin/env python3
"""
Month x arm table of monitor O-F std change vs OL, % (DA-OL)/OL.

For each arm and each month in [--start, --end] (plus the full period), run
score_cygl1_arm.py (or reuse an existing score log in --log-dir), then reduce
the per-species "% vs OL" lines to species-group means (SMOS, SMAP, ASCAT, L3)
plus the CygL1 own-O-F change and N_L1 - the same reduction as the overnight
wait_score_month.sh awk.

Usage:
  build_month_arm_table.py --start 202001 --end 202012 \\
      --arms DA_L1_full_xc015_fixedop:full_xc015_fixedop DA_SMAP_fixedop:smap_fixedop \\
      --log-dir ../../output/month_arm_table_2020 --out-prefix ../../output/month_arm_table_2020/table
"""
import argparse
import os
import re
import subprocess
import sys
from datetime import datetime
from dateutil.relativedelta import relativedelta

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
GROUPS = ['SMOS', 'SMAP', 'ASCAT', 'L3']


def grp(species):
    for g in ('SMOS', 'SMAP', 'ASCAT'):
        if species.startswith(g):
            return g
    return 'L3'


def parse(log, tag):
    txt = open(log).read()
    if '=' * 20 not in txt:
        return None
    txt = txt[txt.index('=' * 20):]
    out = {}
    pct = {g: [] for g in GROUPS}
    for line in txt.splitlines():
        f = line.split()
        if f and f[0] == tag:
            out['L1'] = float(f[-2])
            out['N_L1'] = int(f[-1])
        m = re.match(r'\s*(\S+)\s+OmF_stdv=.*?(-?[\d.]+)%\s*$', line)
        if m:
            pct[grp(m.group(1))].append(float(m.group(2)))
    for g in GROUPS:
        out[g] = sum(pct[g]) / len(pct[g]) if pct[g] else float('nan')
    return out


def md_table(t):
    lines = ['| arm | ' + ' | '.join(t.columns) + ' |',
             '|---|' + '---:|' * len(t.columns)]
    for arm, r in t.iterrows():
        lines.append(f'| {arm} | ' + ' | '.join(f'{v:+.2f}' for v in r) + ' |')
    return '\n'.join(lines)


def score(expid, tag, ol, beg, end_excl, log, rescore):
    if rescore or not os.path.exists(log) or parse(log, tag) is None:
        cmd = [PY, 'score_cygl1_arm.py', '--da-expid', expid, '--arm-tag', tag,
               '--ol-expid', ol, '--start', beg, '--cap-end', end_excl]
        with open(log, 'w') as fh:
            subprocess.run(cmd, cwd=HERE, stdout=fh, stderr=subprocess.STDOUT)
    return parse(log, tag)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--arms', nargs='+', required=True, help='EXPID:TAG ...')
    p.add_argument('--ol-expid', default='OLv8_M36_AZ_fixedop')
    p.add_argument('--start', required=True, help='yyyymm')
    p.add_argument('--end', required=True, help='yyyymm (inclusive)')
    p.add_argument('--log-dir', required=True)
    p.add_argument('--reuse-dirs', nargs='*', default=[],
                   help='other dirs with score_<tag>_<yyyymm|yyyy>.log to reuse')
    p.add_argument('--out-prefix', required=True)
    p.add_argument('--rescore', action='store_true')
    a = p.parse_args()
    os.makedirs(a.log_dir, exist_ok=True)

    first = datetime.strptime(a.start, '%Y%m')
    last = datetime.strptime(a.end, '%Y%m')
    periods = []
    m = first
    while m <= last:
        periods.append((f'{m:%Y%m}', f'{m:%Y%m%d}', f'{m + relativedelta(months=1):%Y%m%d}'))
        m += relativedelta(months=1)
    full_lbl = 'ALL' if (a.start, a.end) != (f'{first:%Y}01', f'{first:%Y}12') else f'{first:%Y}'
    periods.append((full_lbl, f'{first:%Y%m%d}', f'{last + relativedelta(months=1):%Y%m%d}'))

    rows = []
    for arm in a.arms:
        expid, tag = arm.split(':')
        for lbl, beg, end in periods:
            name = f'score_{tag}_{lbl}.log'
            log = os.path.join(a.log_dir, name)
            for d in a.reuse_dirs:
                src = os.path.join(d, name)
                if not os.path.exists(log) and os.path.exists(src) and parse(src, tag):
                    log = src
            r = score(expid, tag, a.ol_expid, beg, end, log, a.rescore)
            if r is None:
                print(f'!! no score for {tag} {lbl} (see {log})', flush=True)
                continue
            rows.append(dict(arm=tag, period=lbl, **r))
            print(f'{tag:40s} {lbl:6s} ' + ' '.join(f'{r[g]:6.2f}' for g in GROUPS)
                  + f' | {r["L1"]:6.2f} {r["N_L1"]}', flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(a.out_prefix + '.csv', index=False)
    order = [lbl for lbl, _, _ in periods]
    with open(a.out_prefix + '.md', 'w') as fh:
        fh.write(f'# Monitor O-F std change vs {a.ol_expid}, % (DA-OL)/OL (negative = better)\n\n')
        for g in GROUPS + ['L1']:
            t = df.pivot(index='arm', columns='period', values=g)
            t = t[[c for c in order if c in t.columns]].reindex([x.split(':')[1] for x in a.arms])
            fh.write(f'## {g}\n\n' + md_table(t) + '\n\n')
    print('wrote', a.out_prefix + '.csv', a.out_prefix + '.md')


if __name__ == '__main__':
    main()
