# Analysis B: key findings (hand-written from `B_summary.md`; run 2026-11-21, SLURM 59364204)

Scripts: `scripts/viirs_transition/incr_stats.py` (catch_progn_incr), `incr_innov.py` (ObsFcstAna), `incr_stats.j`, `incr_summary.py`.
Raw per-arm NetCDFs: `output/incr_stats/` (git-ignored). Numbers are N of 30N, 33 gap/boundary days dropped from every arm, one member.

1. **Double increments are a MODIS Terra+Aqua feature only.** In TA arms with MODIS (TA_to_J1J2 SY1-SY3, TA_control SY4-SY6)
   51-60% of tile-days with an increment have two (never three); in every Terra-only or VIIRS configuration it is 0% (J1 and J2
   share one step and are averaged into one obs per tile). In 97% of doubles the two increments are in consecutive 3-hourly steps; the other 3% are steps 00z and 21z of the same UTC file (180W-135W sector, where the Terra and Aqua passes straddle the UTC date change, so they are really Aqua of one local day and Terra of the next, not a pair).
   TA_control's double share drifts down (55.5% -> 51.3%), TA_to_J1J2's single-step VIIRS arm has 100% single.
2. **Pair behaviour.** Of double tile-days 69-75% have the same sign (reinforcing); pair net/abs is 61-64%, i.e. about a third of
   the absolute adjustment cancels between Terra and Aqua. The second increment is ~10% smaller than the first (0.93-0.98 vs 1.05-1.10
   kg m-2 mean abs pre-switch). Same-sign share falls in TA_control over time (73 -> 69%).
3. **Increment volume.** Annual abs increment per tile (kg m-2 yr-1): T-only ~65-76, TA (MODIS) ~120-138, TA_to_J1J2 post-switch
   79-85 (drops ~35% at the switch; tiles with >=1 increment/day fall 10.4k -> 8.7k), T_to_J2 post-switch 65-70 (no change). Net/abs is
   64-71% in all arms. So switching TA -> J1+J2 removes the second update and with it about a third of the total adjustment; the VIIRS arm is
   closer to a Terra-only arm in volume (80 vs 65-76) than to TA MODIS (120-138).
4. **Cap.** 1.7-3.3% of increments hit the 5 kg m-2 cap (Terra-only 2.0-3.3%, TA 1.7-2.3%, TA_to_J1J2 post-switch 2.1-2.3%); no arm is cap-dominated.
5. **Second-sensor innovations** (B4; tiles with both sensors, RMS OmF in SCF): the Terra update already moves Aqua's forecast toward Aqua obs:
   Aqua RMS 0.32 (OL) -> 0.244 (T_to_J2, Terra assimilated, Aqua not) -> 0.237 (TA, SY1); Terra 0.312 -> 0.264 (assimilated). Same for VIIRS:
   J1 RMS 0.320 (OL) -> 0.2725 (T_to_J2 SY4, J2 only) -> 0.2689 (TA_to_J1J2). So most of the innovation reduction comes from the first
   update, and assimilating the second adds only about 3% RMS. Mean OmF stays positive (+0.06 to +0.10), i.e. obs have more snow than the model in all arms.
6. **UTC pattern** (`B_2`): increments follow local time as expected (Terra step at 10:30 local -> 09z-03z diagonal; Aqua one step later), doubles
   concentrated 45E-135E and 135W-45W.

## Caveats
- Counting is per UTC-day file: a Terra (21z day N) + Aqua (00z day N+1) pair in the 180W-135W sector is split over two files and counted as two singles, so the double share is slightly underestimated there and the 3% above is mis-paired.
- OL writes no `catch_progn_incr` (getinnov only); it only enters B4.
- MODIS innovation rows (Aqua, Terra) are limited to days up to 2025-12-31 in every arm (OL/VIIRS arms have no monitor MODIS in 2026), so
  SY6 MODIS rows cover Aug-Dec 2025 only; use SY6 MODIS means only between arms, not against SY4/SY5. SY6 Aqua/Terra OmF is much larger
  in every arm (different season sample).
- An increment is |sum of 3 layer increments| >= 0.01 kg m-2; the share of "tiny but non-zero" values was not separately assessed. Cap hit = >= 4.995.
- One member, three snow years per contrast: spread across years is shown per SY, but no decision rule has been applied (still unconfirmed).
- Day 2023-08-01 is excluded (segment boundary); the controls start 2023-08-01 and have no increment file that day. Gap-day exclusion is the same for all arms.
- Part of B not done: tile-level maps of double fraction / annual net (data are in `tile_*` variables of `incr_stats_<arm>.nc`).
- Fractional "frac assimilated" in B4 is a day-weighted mean; 1.00 = assimilated all days of that SY.
