#!/usr/bin/env python3
"""
Part B4: does averaging J1 and J2 fit MODIS Aqua and IMS better than J1 or J2 alone? (class 1 dropped, P3-P7)

On the cells where J1 and J2 are both kept (and, for the Aqua comparison, Aqua too, so all three are scored on
the same cells): the cell-level average 0.5 * (J1 + J2) against J1 and J2 alone.
  vs Aqua : RMSD and bias, all matched cells and snowy ones (either > 0)
  vs IMS  : accuracy, false snow, missed snow at SCF >= SCF_ANA_BETA = 0.55 (Northern Hemisphere; nearest IMS pixel
            and the 'decisive' 3x3 window of ims_check.py)
GEOSldas averages per tile (the species' tile means, itself a mean of the kept cells of the tile); the cell-level
average is the same idea one step earlier. This tests the obs, not the update.

Prints the table to stdout (markdown); nothing is written to disk.

usage: jpss_average_check.py
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
from lib import scf_qc  # noqa: E402
from periods import cmg_file, days  # noqa: E402
from jpss_compare import AQUA, J1, J2, PERIODS_B, load_viirs  # noqa: E402
from ims_check import NH_ROWS, cell_pixels, ims_state  # noqa: E402
from get_ims import ims_file  # noqa: E402

NAMES = ('J1', 'J2', 'avg')
THR = np.float32(0.55)


def main():
    r, c = cell_pixels()
    aq = {k: np.zeros(3) for k in NAMES}            # n, sum sq error, sum error (all matched cells)
    aqs = {k: np.zeros(3) for k in NAMES}           # the same, snowy cells
    ims = {d: {k: np.zeros(4) for k in NAMES} for d in ('nn', 'dec')}   # both snow, obs only, IMS only, neither
    for p in PERIODS_B:
        for d in days(p):
            s1, k1, _ = load_viirs(J1, d)
            s2, k2, _ = load_viirs(J2, d)
            sa, ka = scf_qc.modis_scf(scf_qc.read_modis_cmg(cmg_file(AQUA, d)))
            both = k1 & k2
            obs = {'J1': s1, 'J2': s2, 'avg': np.where(both, 0.5 * (s1 + s2), np.nan)}
            m = both & ka
            for name, x in obs.items():
                e = (x[m] - sa[m]).astype(np.float64)
                sn = (x[m] > 0) | (sa[m] > 0)
                aq[name] += [e.size, (e ** 2).sum(), e.sum()]
                aqs[name] += [sn.sum(), (e[sn] ** 2).sum(), e[sn].sum()]
            st = ims_state(ims_file(d), r, c)
            for dn in ('nn', 'dec'):
                valid, isnow = st[dn]
                mm = both[:NH_ROWS] & valid
                for name, x in obs.items():
                    v, sn = x[:NH_ROWS][mm] >= THR, isnow[mm]
                    ims[dn][name] += [(v & sn).sum(), (v & ~sn).sum(), (~v & sn).sum(), (~v & ~sn).sum()]
        print(p, file=sys.stderr, flush=True)

    print('### B4. Average of J1 and J2 vs Aqua (cells where J1, J2 and Aqua are kept), P3-P7\n')
    print('| obs | cells (M) | RMSD | bias | snowy RMSD | snowy bias |\n|---|---|---|---|---|---|')
    for k in NAMES:
        n, sse, se = aq[k]
        ns, ssee, see = aqs[k]
        print(f'| {k} | {n / 1e6:.1f} | {np.sqrt(sse / n):.4f} | {se / n:+.4f} | {np.sqrt(ssee / ns):.4f} | {see / ns:+.4f} |')
    for dn, nm in (('nn', 'nearest IMS pixel'), ('dec', 'decisive 3x3 IMS window')):
        print(f'\n### B4. Average of J1 and J2 vs IMS, SCF >= 0.55, {nm} (Northern Hemisphere cells where J1 and J2 are kept)\n')
        print('| obs | cells (M) | accuracy % | false snow % of obs snow | missed snow % of IMS snow |\n|---|---|---|---|---|')
        for k in NAMES:
            b, o, i, nn = ims[dn][k]
            print(f'| {k} | {(b + o + i + nn) / 1e6:.1f} | {100 * (b + nn) / (b + o + i + nn):.2f} | {100 * o / (b + o):.2f} | {100 * i / (b + i):.2f} |')


if __name__ == '__main__':
    main()
