#!/usr/bin/env python3
"""
Build the dense075 CYGNSS L1 thinning tier for a given date range: a nested
superset of DA-intermediate (min_sep_deg=2.40), further greedily relaxed at
min_sep_deg=0.75 -- lands between DA-intermediate and DA-dense (unthinned) on
the density spectrum, targeting ~1-2 M36 grid cells (~0.36-0.9deg) median
nearest-neighbor spacing among kept obs. Calibrated 2026-09-06 following up on
the DA-intermediate/coh05 obs-spacing diagnosis (see
[[cygl1-paired-thinning-density-experiment]]): sweeping min_sep_deg 1.50 down
to 0.40 gave 0.75 -> N=26,340 (5.0x intermediate's 5,275 over the original
Jan-Jun 2020 range), median NN dist=0.824deg (~1.8-2.3 grid cells),
local-interaction median=12/p90=18 -- the value the user picked.

Nested superset of DA-intermediate (itself a nested superset of DA-sparse):
reproduces DA-intermediate's own kept set deterministically via the same
build_sparse()+build_intermediate_candidate(min_sep=2.40) calls over
whatever date range is passed (greedy admission only depends on same-window
candidates, so this exactly reproduces DA-intermediate's kept obs for that
range with no need to re-derive from already-thinned output files), then
relaxes further at min_sep_deg=0.75. Built from the ORIGINAL full/unthinned
CYGNSS_L1 stream (not the already-thinned intermediate files) so
tile_start/tile_count/support remapping in write_thinned_files() stays
correct.

write_thinned_files() only writes files for the dates passed in, so this can
be (and has been) rerun for successive non-overlapping date ranges to extend
an existing CYGNSS_L1_thinned_dense075_6mo tree without touching earlier
dates already written -- confirmed pattern, same as
thin_cygl1_nested_density_6mo.py's --beg-date/--end-date extension usage.

Consolidates three prior near-duplicate scripts that differed only in
BEG_DATE/END_DATE (build_cygl1_dense075_6mo.py [Jan-Jun 2020],
_extend_jul_dec.py [Jul-Dec 2020], _extend_2021.py [Jan-Dec 2021]) --
merged 2026-09-08, now driven by CLI args instead of hardcoded dates.

Usage: build_cygl1_dense075_thinning.py --beg-date 20200101 --end-date 20200630
"""
import argparse
import sys

sys.path.insert(0, "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/scripts")
import thin_cygl1_nested_density_6mo as thin  # noqa: E402

DST_ROOT = "/discover/nobackup/projects/land_da/cygl1_operator_test/CYGNSS_L1_thinned_dense075_6mo"
MIN_SEP_DEG = 0.75
LABEL = f"nested-superset-of-intermediate, min_sep_deg={MIN_SEP_DEG}, xcompact=ycompact=1.25deg"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beg-date", required=True, help="YYYYMMDD")
    ap.add_argument("--end-date", required=True, help="YYYYMMDD, inclusive")
    args = ap.parse_args()

    thin.BEG_DATE = args.beg_date
    thin.END_DATE = args.end_date

    df, src_paths, dates = thin.load_all_obs()
    sparse_mask = thin.build_sparse(df)
    intermediate_mask = thin.build_intermediate_candidate(df, sparse_mask, min_sep_deg=2.40)
    dense075_mask = thin.build_intermediate_candidate(df, intermediate_mask, min_sep_deg=MIN_SEP_DEG)

    print(f"sparse={sparse_mask.sum()}, intermediate={intermediate_mask.sum()}, "
          f"dense075={dense075_mask.sum()} ({dense075_mask.sum()/intermediate_mask.sum():.2f}x intermediate)")

    assert (intermediate_mask & ~dense075_mask).sum() == 0, "intermediate obs missing from dense075 -- nesting broken!"

    thin.write_thinned_files(df, dense075_mask, src_paths, dates, DST_ROOT, LABEL)
    print(f"Done. Output: {DST_ROOT}")


if __name__ == "__main__":
    main()
