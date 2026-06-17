# ASCAT/SMAP Paper Figure Provenance

This note maps the older ASCAT/SMAP paper figure machinery into the current
project layout. The workflow predates the `projects/ascat_da` reorganization,
so the final paper figure notebook lives under `projects/utils`, while earlier
ASCAT-specific diagnostics now live here.

## Final figure directory

Local draft/revision figures are in:

`/Users/amfox/Documents/Publications/ASCAT_SMAP_paper/Revisions/Draft_figures`

That directory currently contains 69 PNG figure files. File timestamps suggest
three production/revision batches:

- 2024-06-26 around 12:28-12:30: main map/time-series figure batch.
- 2024-07-01 around 15:45-16:04: revised in-situ/Rdiff and relative-stddev maps.
- 2025-03-13 to 2025-03-31: revised/300dpi time-series and in-situ updates.

## Experiment names

The final publication-style notebook uses the `fp_scaled` experiment set:

| Label | Experiment id |
| --- | --- |
| `CNTL` | `OLv7_M36_MULTI_type_13_comb_fp_scaled` |
| `ASC_DA` | `DAv7_M36_ASCAT_type_13_comb_fp_scaled` |
| `SMP_DA` | `DAv7_M36_SMAP_type_13_comb_fp_scaled` |
| `MLT_DA` | `DAv7_M36_MULTI_type_13_comb_fp_scaled` |

The older ASCAT-project precursor used the unscaled/earlier `comb_fp` set:

| Label | Experiment id |
| --- | --- |
| `OL` | `OLv7_M36_ASCAT_type_13_comb_fp` |
| `ASCAT` | `DAv7_M36_ASCAT_type_13_comb_fp` |
| `SMAP` | `DAv7_M36_SMAP_type_13_comb_fp` |
| `MULTI` | `DAv7_M36_MULTI_type_13_comb_fp` |

## Main notebooks and scripts

| Path | Role |
| --- | --- |
| `projects/utils/notebooks/paper_figures_083024.ipynb` | Main paper figure notebook. Produces the final `fp_scaled` map, time-series, in-situ, and Rdiff-style figures. |
| `projects/utils/notebooks/preliminary_paper_figures_081424.ipynb` | Earlier figure assembly notebook with overlapping cells and older output styling. |
| `projects/ascat_da/notebooks/compare_comb_fp_043024.ipynb` | ASCAT-project precursor for combined-footprint experiment diagnostics. Uses `comb_fp` rather than `fp_scaled`. |
| `projects/utils/scripts/mapper_functions.py` | Plot helper that auto-generates many PNG filenames from plot titles. |
| `projects/matlab2python/matlab_postprocess/IVs/Plot_script_sample.m` | Older MATLAB IV workflow that writes `Rdiff_*.mat` products used by the paper notebook. |
| `projects/matlab2python/matlab_postprocess/in_situ/SMskill_vs_INSITU_single_expt.m` | Older MATLAB in-situ workflow that writes `*_stats.mat` products used by the in-situ bar figures. |

## Filename generation

Many `Draft_figures` filenames are not literal `savefig("...")` calls. They
come from `plot_global_tight()` or `plot_global_tight_pcm()` in
`projects/utils/scripts/mapper_functions.py`, which saves figures with:

```python
savename = re.sub('[^0-9a-zA-Z]+', '_', plot_title) + '.png'
```

For example:

| Plot title | Saved filename |
| --- | --- |
| `SMAP OmF (Tb_V): 2020701_0300z` | `SMAP_OmF_Tb_V_2020701_0300z.png` |
| `ASCAT OmF (SFDS): 2020701_0300z` | `ASCAT_OmF_SFDS_2020701_0300z.png` |
| `MLT_DA - CNTL: $\Delta$ StdDev of OmF (Tb)` | `MLT_DA_CNTL_Delta_StdDev_of_OmF_Tb_.png` |
| `ASC_DA: Number of surface SM increments` | `ASC_DA_Number_of_surface_SM_increments.png` |

## Figure groups

### One-day snapshot maps

Likely source: `projects/utils/notebooks/paper_figures_083024.ipynb`, with
overlapping earlier cells in `preliminary_paper_figures_081424.ipynb`.

These use local daily files under `../test_data/fp_scaled`, including
`DAv7_M36_MULTI_type_13_comb_fp_scaled.inst3_1d_lndfcstana_Nt.20200701.nc4`
and `DAv7_M36_MULTI_type_13_comb_fp_scaled.catch_progn_incr.20200701.nc4`.

Representative outputs:

- `SMAP_OmF_Tb_H_2020701_0300z.png`
- `SMAP_OmF_Tb_V_2020701_0300z.png`
- `ASCAT_OmF_SFDS_2020701_0300z.png`
- `srfexc_Increment_2020701_0300z.png`
- `rzexc_Increment_2020701_0300z.png`
- `SFMC_Analysis_2020701_0300z.png`
- `RZMC_Analysis_2020701_0300z.png`
- `Surface_Wetness_Analysis_2020701_0300z.png`
- `Rootzone_Wetness_Analysis_2020701_0300z.png`

### OmF time series

Likely source: `projects/utils/notebooks/paper_figures_083024.ipynb`.

Inputs are daily `*_OmF_ts.npz` products under `../test_data/fp_scaled`.
The notebook contains explicit `savefig()` calls for:

- `OL_ASCAT_SMAP_MULTI_fp_scaled_SFDS_OmF_ts_stddev_only.png`
- `OL_ASCAT_SMAP_MULTI_fp_scaled_Tb_OmF_ts_stddev_only.png`
- `OL_ASCAT_SMAP_MULTI_fp_scaled_SFDS_OmF_ts_stddev_only_revised_300dpi.png`
- `OL_ASCAT_SMAP_MULTI_fp_scaled_Tb_OmF_ts_stddev_only_revised_300dpi.png`

The non-300dpi `_revised.png` copies in `Draft_figures` were not found as
literal save names in the notebook and may have been manual exports or renames.

### OmF stddev and relative-stddev maps

Likely source: `projects/utils/notebooks/paper_figures_083024.ipynb`.

Inputs are monthly/stat arrays from the `fp_scaled` diagnostics. Outputs are
auto-named from plot titles, including:

- `CNTL_StdDev_of_OmF_SFDS_.png`
- `CNTL_StdDev_of_OmF_Tb_.png`
- `ASC_DA_CNTL_Delta_StdDev_of_OmF_SFDS_.png`
- `ASC_DA_CNTL_Delta_StdDev_of_OmF_Tb_.png`
- `SMP_DA_CNTL_Delta_StdDev_of_OmF_SFDS_.png`
- `SMP_DA_CNTL_Delta_StdDev_of_OmF_Tb_.png`
- `MLT_DA_CNTL_Delta_StdDev_of_OmF_SFDS_.png`
- `MLT_DA_CNTL_Delta_StdDev_of_OmF_Tb_.png`
- `MLT_DA_ASC_DA_Delta_StdDev_of_OmF_SFDS_.png`
- `MLT_DA_SMP_DA_Delta_StdDev_of_OmF_Tb_.png`
- Relative-stddev variants with `Delta_Relative_StdDev_of_OmF_*`.

### Increment count and increment stddev maps

Likely source: `projects/utils/notebooks/paper_figures_083024.ipynb`.

Inputs are `*_increment_stats.npz` products under `../test_data/fp_scaled`.
Outputs include:

- `ASC_DA_Number_of_increments_0_4.png`
- `SMP_DA_Number_of_increments_0_2.png`
- `MLT_DA_Number_of_increments_0_5.png`
- `ASC_DA_Number_of_surface_SM_increments.png`
- `SMP_DA_Number_of_surface_SM_increments.png`
- `MLT_DA_Number_of_surface_SM_increments.png`
- `ASC_DA_StdDev_of_srfexc_increment.png`
- `ASC_DA_StdDev_of_rzexc_increment.png`
- `SMP_DA_StdDev_of_srfexc_increment.png`
- `SMP_DA_StdDev_of_rzexc_increment.png`
- `MLT_DA_StdDev_of_srfexc_increment.png`
- `MLT_DA_StdDev_of_rzexc_increment.png`
- surface/root-zone soil-moisture increment stddev and delta maps.

### In-situ soil-moisture skill bar figures

Likely source: `projects/utils/notebooks/paper_figures_083024.ipynb`.

Inputs are MATLAB stats files under `../test_data/fp_scaled`, with names like:

`{experiment_id}_CalVal_M33_SM_3h__6yr_stats.mat`

The upstream MATLAB workflow appears to be
`projects/matlab2python/matlab_postprocess/in_situ/SMskill_vs_INSITU_single_expt.m`.

Outputs:

- `fp_scaled_surf_rz_stats.png`
- `fp_scaled_surf_rz_stats_revised.png`
- `fp_scaled_surf_rz_stats_revised_300dpi.png`

### Surface soil-moisture anomaly-R maps

Likely source: `projects/utils/notebooks/paper_figures_083024.ipynb`, using
Rdiff `.mat` products under `../test_data/fp_scaled`.

The upstream Rdiff products are generated by the older MATLAB/Python IV
workflow under `projects/matlab2python/matlab_postprocess/IVs/` and
`projects/matlab2python/scripts/plot_script_sample.py`.

Inputs include:

- `Rdiff_DAv7_M36_ASCAT_type_13_comb_fp_scaled_minus_OLv7_M36_MULTI_type_13_comb_fp_scaled.mat`
- `Rdiff_DAv7_M36_SMAP_type_13_comb_fp_scaled_minus_OLv7_M36_MULTI_type_13_comb_fp_scaled.mat`

Outputs include:

- `ASC_DA_minus_CNTL_Surface_Soil_Moisture_Skill_anomaly_R_RB.png`
- `ASC_DA_minus_CNTL_Surface_Soil_Moisture_Skill_anomaly_R_RBr.png`
- `SMP_DA_minus_CNTL_Surface_Soil_Moisture_Skill_anomaly_R_RB.png`
- `SMP_DA_minus_CNTL_Surface_Soil_Moisture_Skill_anomaly_R_RBr.png`

The `RB` versus `RBr` suffix difference likely reflects color-map/style
revision rather than a different metric, but this has not been fully verified.

## ASCAT project precursor

`projects/ascat_da/notebooks/compare_comb_fp_043024.ipynb` is related to this
paper era but does not appear to be the final source of `Draft_figures`.

It uses `comp_dir = 'comb_fp'` and the earlier experiment names without the
`_scaled` suffix. It includes:

- ASCAT SFDS OmF time-series diagnostics.
- SMAP Tb OmF time-series diagnostics.
- ASCAT/SMAP observation-count time series.
- Cross-masked observation-count comparisons.
- Increment and stddev diagnostics.

This notebook is best treated as a precursor to the later `fp_scaled` paper
figure notebook.

## Known gaps and caveats

- The final paper figure workflow is notebook-first and order-dependent. It is
  not currently packaged as a reproducible script.
- Many inputs are referenced by relative paths such as `../test_data/fp_scaled`;
  moving notebooks would likely break reruns.
- Several figure files in `Draft_figures` may be manual exports, copies, or
  renames rather than direct notebook save outputs.
- The exact distinction between `*_RB.png` and `*_RBr.png` Rdiff outputs needs
  visual or notebook-output verification.
- The `projects/ascat_da` project appears to have been carved out after the
  2024 paper work. Git history shows this project area being added later than
  the Draft_figures timestamps, which explains why the main figure notebook
  still lives in `projects/utils`.
