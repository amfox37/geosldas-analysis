# scripts/viirs_transition

Scripts for `docs/viirs_transition_plan.md`. All of them read `config/arms.json` (arms, dates, species, paths), run on
Discover with the GEOSpyD py3.13 python (`/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3`, no
`g5_modules`), and write compact files under `output/` (git-ignored, on `land_da`).

## Analysis A: observation availability (done 2026-11-21)

| File | What it does |
|---|---|
| `obs_availability.py` | One arm: reads every 3-hourly `ana/ens_avg/*ldas_ObsFcstAna*` file and writes `output/obs_availability/obs_availability_<arm>.nc` |
| `obs_availability_discover.sh` | Runs it for all arms (or `ARMS="..."`; test range with `BEG=`/`END=`/`OUT_DIR=`) |
| `obs_availability.j` | SLURM array, one task per arm of `arms.json` (account `g0610`, as in `main_run/bat_inp_main.txt`). Submit from this directory: `sbatch obs_availability.j`; delete the `obs_avail.*.out/err` afterwards |
| `obs_availability_summary.py` | Reads the five NetCDFs; writes the tables, checks and figures below |

Run: `sbatch obs_availability.j` (about 12 min), then
`python3 obs_availability_summary.py` (a few seconds; set `MPLCONFIGDIR` to a writable directory under `land_da`).

Timing (job 59362332): one arm-month takes about 8 s; the three 6-year arms took 11 min 11 s to 11 min 24 s, the two
3-year controls 6 min. The array was submitted with `--time=00:12:00` (about 1.2 x the estimate) and finished only
about 40 s inside it. Use 15 minutes if you rerun, because the run time varies with file-system load.

### Output NetCDF (`obs_availability_<arm>.nc`)

Species are mapped through `obsparam_descr` of each file (never the namelist indices): `species_name` is
Aqua `MYD10C1`, Terra `MOD10C1`, J1 `VJ110C1`, J2 `VJ210C1` (local ids 1-4 inside `ObsFcstAna`).
Days are the UTC date of the file time stamp; bands are `S_hemisphere`, `N_0_30`, `N_30_45`, `N_45_60`, `N_60_90`
(by tile latitude).

| Variable | Dimensions | Meaning |
|---|---|---|
| `n_obs` | day, species, band | obs in all 3-hourly files of the day (monitored or assimilated) |
| `n_assim` | day, species, band | of these, obs with `assim_flag == 1` |
| `n_tiles` | day, species, band | unique tiles with at least one obs (equals `n_obs`: a tile gets one obs of a species per day) |
| `n_steps` | day | 3-hourly files found (8 normally; 7 on 2020-08-01, which has no 00z file, and 1 on the last day) |
| `tiles_both_terra_aqua`, `tiles_both_j1_j2` | day | tiles with both sensors that day |
| `tiles_any_modis`, `tiles_any_viirs` | day | tiles with either sensor |
| `tiles_terra_only`, `tiles_aqua_only`, `tiles_j1_only`, `tiles_j2_only` | day | tiles with one of the pair only |
| `assim_any` | day, species | 1 if any obs of the species was assimilated that day |
| `days_with_obs` | month, species, tile | days in the month with at least one obs (int16; tile index = `tilenum - 1`) |
| `tile_lat`, `tile_lon` | tile | position of tiles that ever had an obs |

### Summary tables, checks, figures

| File | Where | Content |
|---|---|---|
| `summary_obs_per_day.csv`, `summary_overlap.csv`, `zero_obs_days.csv` | `output/obs_availability/` | mean obs, assimilated obs and overlap per day by arm, snow year and species; runs of days with no obs of a species |
| `obs_availability_summary.md`, `obs_availability_checks.md` | `docs/viirs_transition_report_figures/` | the same tables, as markdown, plus the check results |
| `figA_1_tiles_per_day_by_species.png`, `figA_2_constellation_overlap.png`, `figA_3_tiles_by_latitude_band.png` | `docs/viirs_transition_report_figures/` | monthly tiles/day by species and arm; same-day overlap of Terra+Aqua and J1+J2; tiles/day by latitude band |

Note: `.gitignore` ignores `*.png` and un-ignores only specific figure directories, and `docs/viirs_transition_report_figures/`
is not one of them. The three PNGs (and the markdown, which is not ignored) are therefore not picked up by `git add`
as it stands. I did not edit `.gitignore`.

### Checks run (all pass)

- A single-day count for `TA_to_J1J2` 2024-01-15 12z matches an independent `bincount` of the files and the numbers seen
  earlier (2741 / 510 / 3964 / 3765 for Aqua / Terra / J1 / J2).
- `n_obs == n_tiles` on every day, in all five arms (no tile has two obs of one species in a day).
- Assimilated species have all their obs assimilated, and monitor-only species have none, on every full day in every
  arm (the segment pattern of `arms.json`; 2023-08-01 is skipped, see below).

### Things to know when using the output

- **2023-08-01 is a mixed day.** The 00z file of that day is the analysis at the end of segment 1 and still assimilates
  the old sensors (Terra, or Terra+Aqua); the other steps use the new configuration. For example `T_to_J2` on that day
  has 1221 assimilated Terra and 4144 assimilated J2 obs.
- **Partial days.** 2020-08-01 (7 steps) and 2026-08-01 (1 step) are not full days; the summary skips days with fewer
  than 8 steps.
- **Data gaps (missing input days, not model problems):** Aqua 2022-04-02 to 2022-04-16 (15 days), Terra 2022-10-11 to
  2022-10-22 (12 days), Terra 2025-05-11 to 2025-05-14 (4 days) and J2 2024-02-28 (1 day). The affected sensor is
  assimilated in the arms that use it, so those arms had no assimilation from it on those days (`zero_obs_days.csv`).
- **2026 MODIS monitoring gap.** The OL, `T_to_J2` and `TA_to_J1J2` have no Aqua or Terra obs from 2026-01-01 (Terra) /
  2026-01-02 (Aqua) to the end (`gmao_ops` ends 2025-12-31); the controls do. This is recorded as zeros, not as an error.
  It makes the SY6 rows of Aqua and Terra in the summary tables low for those three arms.
