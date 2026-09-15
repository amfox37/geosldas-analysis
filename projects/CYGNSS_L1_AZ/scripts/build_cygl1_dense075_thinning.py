#!/usr/bin/env python3
"""
Build a nested CYGNSS L1 thinning tier for a given date range and min_sep_deg:
a nested superset of DA-intermediate (min_sep_deg=2.40), further greedily
relaxed at the given --min-sep-deg -- lands between DA-intermediate and
DA-dense (unthinned) on the density spectrum. Calibration sweep (2026-09-06,
see [[cygl1-paired-thinning-density-experiment]], sweeping min_sep_deg 1.50
down to 0.40) over the original Jan-Jun 2020 range:

  min_sep_deg | N      | vs intermediate | median NN dist | local-interaction med/p90
  1.25        | 12,756 | 2.4x            | 1.336deg        | 4/7
  1.00        | 17,842 | 3.4x            | 1.091deg        | 7/11
  0.75        | 26,340 | 5.0x            | 0.824deg        | 12/18  (the "dense075" tier)
  0.60        | 34,284 | 6.5x            | 0.681deg        | 16/25
  0.50        | 42,044 | 8.0x            | 0.578deg        | 21/32  (the "dense050" tier)
  0.40        | 53,613 | 10.2x           | 0.473deg        | 27/42

Nested superset of DA-intermediate (itself a nested superset of DA-sparse):
reproduces DA-intermediate's own kept set deterministically via the same
build_sparse()+build_intermediate_candidate(min_sep=2.40) calls over
whatever date range is passed (greedy admission only depends on same-window
candidates, so this exactly reproduces DA-intermediate's kept obs for that
range with no need to re-derive from already-thinned output files), then
relaxes further at --min-sep-deg. Built from the ORIGINAL full/unthinned
CYGNSS_L1 stream (not the already-thinned intermediate files) so
tile_start/tile_count/support remapping in write_thinned_files() stays
correct.

write_thinned_files() only writes files for the dates passed in, so this can
be (and has been) rerun for successive non-overlapping date ranges to extend
an existing CYGNSS_L1_thinned_<tier>_6mo tree without touching earlier
dates already written -- confirmed pattern, same as
thin_cygl1_nested_density_6mo.py's --beg-date/--end-date extension usage.

Generalized 2026-09-09 (per [[feedback_avoid_per_run_script_copies]]) from the
original dense075-only version to take --min-sep-deg/--tier-name as CLI args,
so a new density tier (e.g. dense050) doesn't need a new script file. Output
dir is derived as CYGNSS_L1_thinned_<tier-name>_6mo; --tier-name defaults to
a slug of --min-sep-deg (e.g. 0.50 -> "dense050") when omitted.

Usage: build_cygl1_dense075_thinning.py --min-sep-deg 0.75 \\
           --beg-date 20200101 --end-date 20200630 [--tier-name dense075]
"""
import argparse
import sys

sys.path.insert(0, "/gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ/scripts")
import thin_cygl1_nested_density_6mo as thin  # noqa: E402

BASE = "/discover/nobackup/projects/land_da/cygl1_operator_test/"


def default_tier_name(min_sep_deg):
    return f"dense{int(round(min_sep_deg * 100)):03d}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-sep-deg", type=float, required=True)
    ap.add_argument("--tier-name", default=None, help="defaults to a slug of --min-sep-deg, e.g. 0.50 -> dense050")
    ap.add_argument("--beg-date", required=True, help="YYYYMMDD")
    ap.add_argument("--end-date", required=True, help="YYYYMMDD, inclusive")
    args = ap.parse_args()

    tier_name = args.tier_name or default_tier_name(args.min_sep_deg)
    dst_root = f"{BASE}CYGNSS_L1_thinned_{tier_name}_6mo"
    label = f"nested-superset-of-intermediate, min_sep_deg={args.min_sep_deg}, xcompact=ycompact=1.25deg"

    thin.BEG_DATE = args.beg_date
    thin.END_DATE = args.end_date

    df, src_paths, dates = thin.load_all_obs()
    sparse_mask = thin.build_sparse(df)
    intermediate_mask = thin.build_intermediate_candidate(df, sparse_mask, min_sep_deg=2.40)
    tier_mask = thin.build_intermediate_candidate(df, intermediate_mask, min_sep_deg=args.min_sep_deg)

    print(f"sparse={sparse_mask.sum()}, intermediate={intermediate_mask.sum()}, "
          f"{tier_name}={tier_mask.sum()} ({tier_mask.sum()/intermediate_mask.sum():.2f}x intermediate)")

    assert (intermediate_mask & ~tier_mask).sum() == 0, f"intermediate obs missing from {tier_name} -- nesting broken!"

    thin.write_thinned_files(df, tier_mask, src_paths, dates, dst_root, label)
    print(f"Done. Output: {dst_root}")


if __name__ == "__main__":
    main()
