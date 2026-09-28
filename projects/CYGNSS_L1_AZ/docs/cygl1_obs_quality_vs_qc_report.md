# CYGNSS L1 observation quality vs raw-granule QC fields (fixed operator build)

*2026-09-26. Arizona limited domain (lon -118..-106, lat 29..40), EASEv2 M36, Jan–Jun 2020.*

## Summary

Every CYGNSS L1 observation assimilated in the full-stream DA arm was matched to its raw-granule QC record. Observation quality was then
stratified by `coherency_ratio` and other QC fields, using three independent measures:

1. innovation spread against an independent open loop,
2. the Desroziers estimate of observation-error variance, and
3. how strongly each L1 innovation co-varies with collocated SMAP Tb and CYGNSS L3 soil-moisture innovations.

**Quality is U-shaped in `coherency_ratio`.** The bottom decile (< 0.40) and the top decile (> 2.16):

- have 2–3× the Desroziers error variance of the middle 80%,
- have 1.3–1.8× its innovation standard deviation, and
- carry the least soil-moisture signal (weakest correlation with the SMAP Tb and L3 innovations).

The high tail is also strongly positively biased (+1.35 dB). After per-tile z-score scaling it receives the *smallest* assumed error variance, so
the least informative observations get the most weight. The tails are not confined to a few tiles, so a per-observation screen is needed rather
than a tile mask. DDM SNR is a second, partly independent quality discriminator. Removing both coherency tails (keeping 0.40–2.16, 80% of
observations):

- lowers the Desroziers error variance from 7.31 to 5.74,
- strengthens the SMAP Tb correlation from −0.187 to −0.224 and the L3 correlation from 0.235 to 0.251,
- but shifts the mean innovation from −0.25 to −0.44 dB, so the scaling climatology must be rebuilt on the screened population.

## Data and method

**Runs** (GEOSldas, fixed L1 operator build md5 6a70b496, commits cdfaaf4/11dfdb1 of `@GEOSldas_GridComp`,
branch `feature/amfox/cygnss-ascat-hsaf-v8`):

- DA: `DA_L1_full_xc015_fixedop`. Full-stream CYGNSS L1 (`CYGNSS_L1_DDM3X5_CROP_SCALAR`, obs_param index 56) assimilated with errstd 2.75 dB,
  xcorr = ycorr = 0.15°, z-score scaling (clim from `OLv8_M36_AZ_fixedop` 2020–2022), 24 members, 3-h analysis. All other species
  (SMOS/SMAP Tb, ASCAT A/B/C, CYGNSS L3 6-hr SM) are monitored only.
- OL: `OLv8_M36_AZ_fixedop` (unscaled, all species monitored, same restart, forcing and parameters).

**Observations:** the 86,456 L1 observations the DA run actually assimilated (one per owner tile per 3-h window, selected by the reader as the
specular point nearest the tile center).

**Join** (exact keys only, no value matching; ~100% match: 86,456 of 86,459):

1. DA ObsFcstAna L1 row (window T, lon, lat) → preprocessed obs-file row (window of its own timestamp, `sp_lon`/`sp_lat` as float32).
   This is exact because reader fix 11dfdb1 writes the selected observation's specular-point lon/lat into ObsFcstAna.
2. obs-file row (`year`, `day`, `sc_num`, `sample_id`, `ch_id`) → per-day, per-spacecraft QC-pass CSV
   (`CYGNSS_operator/artifacts/out_images/cygnss_qc_m36_window_counts_<yyyymmdd>_cyg<NN>/`). This supplies `coherency_state`,
   `coherency_ratio`, `ddm_snr`, `sp_rx_gain`, `srtm_slope`, `modis_land_cover`, `pekel_sp_water_percentage_5km` and `brcs_*`.
3. DA row → OL row on (T, species, tile, lon, lat) → OL forecast `F_OL`.

**Measures** (O = the DA run's scaled observation throughout):

| Quantity | Definition | What it measures |
|---|---|---|
| inn | O − F_OL | Innovation against a DA-independent forecast |
| inn std / mean | std and mean of inn | Total mismatch and bias |
| R_desroz | E[(O−A)(O−F)] from the DA run | Implied observation-error variance (Desroziers et al. 2005) |
| R_assumed | mean scaled `obsvar` | Error variance the filter actually used |
| HPH_desroz | E[(A−F)(O−F)] | Implied forecast-error variance in observation space |
| r_SMAPh | corr(inn, same-tile SMAP Tbh A/D innovation) | Shared soil-moisture signal (expected negative: wetter → higher L1 dB, lower Tb) |
| r_L3 | corr(inn, same-tile CYGNSS L3 innovation) | Shared soil-moisture signal (expected positive) |

The SMAP and L3 innovations use the DA run's scaled monitor observations minus the OL forecast, matched to the L1 observation on the same tile
at the nearest analysis time within ±12 h. With N ≈ 5–6k pairs per bin, the standard error of r is about 0.013.

## Results

### coherency_ratio deciles, Jan–Jun 2020

| coherency_ratio | N | inn mean | inn std | R_desroz | R_assumed | HPH_desroz | r_SMAPh | r_L3 |
|---|---|---|---|---|---|---|---|---|
| 0.07–0.40 | 8646 | −0.40 | **3.96** | **15.43** | 9.34 | 0.40 | **−0.158** | 0.201 |
| 0.40–0.53 | 8645 | −0.83 | 2.81 | 8.24 | 8.78 | 0.24 | −0.248 | 0.227 |
| 0.53–0.66 | 8646 | −0.78 | 2.45 | 6.29 | 8.30 | 0.24 | −0.277 | 0.237 |
| 0.66–0.81 | 8645 | −0.73 | 2.51 | 6.50 | 8.26 | 0.24 | −0.255 | 0.222 |
| 0.81–0.99 | 8646 | −0.64 | 2.27 | 5.22 | 7.62 | 0.21 | −0.252 | 0.251 |
| 0.99–1.22 | 8645 | −0.48 | 2.16 | 4.61 | 7.33 | 0.21 | −0.267 | 0.284 |
| 1.22–1.48 | 8645 | −0.32 | 2.20 | 4.64 | 7.00 | 0.20 | −0.226 | 0.273 |
| 1.48–1.78 | 8646 | −0.08 | 2.24 | 4.75 | 6.59 | 0.20 | −0.209 | 0.280 |
| 1.78–2.16 | 8646 | +0.36 | 2.43 | 5.70 | 6.07 | 0.28 | −0.172 | 0.253 |
| 2.16–10.17 | 8645 | **+1.35** | **3.26** | **11.74** | **5.88** | **0.69** | **−0.088** | 0.175 |

Key points:

- **Both tails are worse on every measure.** The middle 80% has an innovation std of 2.2–2.8 dB and a Desroziers R of 4.6–8.2. The bottom tail
  has 3.96 dB and 15.4, the top tail 3.26 dB and 11.7.
- **Soil-moisture signal is weakest in the tails.** The SMAP Tb correlation in the top decile (−0.088) is a third of the middle's peak
  (−0.277). That gap is about 10 standard errors of the difference.
- **Mean innovation rises monotonically with coherency ratio**, from −0.8 to +1.35 dB. Strongly coherent reflections read as "wet" compared
  with the model.
- **The assumed R falls monotonically with coherency ratio** (9.3 → 5.9). This comes from the per-tile z-score scaling of the error variance.
  So the top tail, which has the second-largest true error, gets the smallest assumed error. Its HPH_desroz (0.69, about 3× the middle)
  confirms that it drives disproportionately large analysis increments.
- **Middle deciles are over-assumed.** In the 0.8–1.8 deciles the Desroziers R (4.6–5.2) is below the assumed R (6.6–7.6).

### Season

The same U-shape appears in both halves of the period (CSV tables `coherency_deciles_JFM` and `coherency_deciles_AMJ`):

| | Bottom decile: R_desroz / r_SMAPh | Middle deciles: r_SMAPh | Top decile: R_desroz / r_SMAPh |
|---|---|---|---|
| Jan–Mar | 13.1 / −0.22 | −0.26 to −0.35 | 9.3 / −0.13 |
| Apr–Jun | 17.6 / −0.08 | −0.13 to −0.19 | 13.4 / −0.07 |

The shared SMAP signal is weaker in spring in every decile. This matches the spring loss of monitor skill seen in the DA scores and the
noise-vs-gain split.

### coherency_state

| coherency_state | fraction | inn mean | inn std | R_desroz | R_assumed | r_SMAPh | r_L3 |
|---|---|---|---|---|---|---|---|
| 0 | 84.1% | −0.50 | 2.65 | 6.94 | 7.86 | −0.212 | 0.240 |
| 1 | 13.8% | +1.14 | 3.11 | 10.34 | 5.89 | −0.108 | 0.186 |
| 2 | 2.1% | +0.23 | 1.50 | 2.14 | 4.53 | −0.317 | 0.287 |

State 1 closely matches the poor high-`coherency_ratio` tail. State 2 is a small, very clean subset. The product definitions of the states still
need to be taken from the CYGNSS L1 documentation before citing them.

### Other QC fields (quintiles)

- **`ddm_snr`: monotonic, and the strongest single discriminator after coherency.**

  | ddm_snr (dB) | 2.0–3.3 | 3.3–5.1 | 5.1–7.5 | 7.5–10.9 | 10.9–34.5 |
  |---|---|---|---|---|---|
  | R_desroz | 12.2 | 8.5 | 6.9 | 4.9 | 4.0 |
  | r_SMAPh | −0.13 | −0.18 | −0.23 | −0.23 | −0.19 |
  | r_L3 | 0.24 | 0.24 | 0.24 | 0.28 | 0.27 |

- **`srtm_slope`:** R_desroz is higher on steeper terrain (5.7 → 8.6–8.9) but the soil-moisture signal is nearly flat (r_SMAPh −0.16 to −0.20).
  Slope is a weaker discriminator.
- **`sp_inc_angle`:** no clear pattern (R_desroz 6.7–8.0, r_SMAPh −0.17 to −0.21).
- **`pekel_sp_water_percentage_5km`:** zero for every observation in this domain, so it is uninformative here.

### coherency × SNR

corr(coherency_ratio, ddm_snr) = 0.43.

| coherency class | SNR tercile | N | inn std | R_desroz | r_SMAPh | r_L3 |
|---|---|---|---|---|---|---|
| low (< 0.40) | low | 5428 | 4.26 | 17.97 | −0.105 | 0.218 |
| low | mid | 3055 | 3.32 | 10.72 | −0.295 | 0.191 |
| low | high | 221 | 3.70 | 15.93 | −0.514 | 0.130 |
| mid (0.40–2.16) | low | 21231 | 2.87 | 8.62 | −0.176 | 0.257 |
| mid | mid | 24356 | 2.40 | 5.59 | −0.241 | 0.256 |
| mid | high | 23469 | 1.88 | 3.30 | −0.277 | 0.289 |
| high (> 2.16) | low | 2160 | 3.84 | 14.67 | −0.127 | 0.160 |
| high | mid | 1407 | 3.92 | 17.46 | −0.117 | 0.197 |
| high | high | 5129 | 2.71 | 8.88 | −0.078 | 0.185 |

- Within the middle coherency band, SNR still separates quality strongly (R_desroz 8.6 → 3.3).
- In the high-coherency tail, high SNR lowers the error variance but *not* the lack of soil-moisture signal (r_SMAPh −0.08). Those
  observations are bright and clean, but they are not measuring the model's soil moisture.
- The low-tail/high-SNR cell (N = 221) is too small to interpret.

### Spatial and land-cover structure

- **The tails are not tile artifacts.** Only 13 (low) and 12 (high) of 569 tiles have a majority of tail observations, and those tiles hold only
  7% and 11% of the tail observations. The 10% of tiles with the most tail observations hold 36% and 47% of them. The tails are somewhat
  concentrated but spread across the domain, so a per-observation screen is needed.
- **IGBP land cover** (share within each class):

  | Land cover | low tail | middle | high tail |
  |---|---|---|---|
  | 7 open shrubland | 64% | 71% | 60% |
  | 16 barren | 21% | 15% | 27% |
  | 12 cropland | 0.3% | 2.9% | 6.0% |
  | 10 grassland | 11% | 9% | 4% |

  The high tail is enriched in barren and cropland (smooth or flooded surfaces); the low tail is enriched in grassland.

### Two-sided screen (keep 0.40 ≤ coherency_ratio ≤ 2.16)

| subset | N | inn mean | inn std | R_desroz | R_assumed | r_SMAPh | r_L3 |
|---|---|---|---|---|---|---|---|
| all | 86,456 | −0.25 | 2.76 | 7.31 | 7.52 | −0.187 | 0.235 |
| keep 0.40–2.16 | 69,056 (80%) | **−0.44** | 2.42 | **5.74** | 7.49 | **−0.224** | **0.251** |

The screen improves every quality measure. The bias shift happens because the high tail was cancelling the negative bias of the rest; with the
original scaling climatology, the screen would add a drying bias. This is why the pre-fix screening experiment (see below) went one-sided. The
proper fix is to rebuild the z-score climatology on the screened population.

## Caveats

- One domain (arid and semi-arid SW US), one half-year, one DA configuration. Desroziers estimates depend on the DA run's gain, and the
  observation-error correlations are known to be misspecified (see `cygl1_obs_error_correlation_report.md`), so absolute R values are
  indicative. The *relative* ordering across bins is the robust result. The innovation std and the monitor correlations use the OL forecast
  only and do not depend on the DA configuration.
- The monitor correlations measure shared innovation (forecast-error) signal. They include forecast-error covariance by design, and up to 12 h
  of time mismatch and representativeness differences (L1 specular footprint vs 36-km tile).
- Only reader-selected observations are analysed. Observations that were never selected (not nearest to the owner-tile center) are not
  represented.
- `coherency_state` meanings are not yet confirmed from product documentation.

## Relation to earlier work

The 2026-08 coherency-screening experiment (one-sided, keep ≥ 0.523, gated binary, errstd 3.0) was run on the **pre-fix operator**, so its L1
numbers are invalid. The later operator fix (cdfaaf4/11dfdb1) changed which forecast is compared with which observation. The present analysis is
the first on the fixed build. It replaces the one-sided choice with a two-sided screen plus a rebuilt climatology.

## Follow-up

The DA test arm `DA_L1_full_xc015_coh040216_fixedop` (2020-01..12) is identical to `DA_L1_full_xc015_fixedop` except that it uses
coherency-screened obs files (`CYGNSS_L1_coh040_216/`, built by `scripts/cygl1_coherency_filter.py`) and a scaling climatology rebuilt on the
screened OL population (`scaling_params/cygnss_l1_z_score_clim_coh040216/`). errstd stays at 2.75 dB (one intervention). Note that the rebuilt
climatology's smaller obs std raises the effective scaled R somewhat even at unchanged errstd.

## Reproduce

```bash
PY=/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3
cd /gpfsm/dnb06/projects/p284/geosldas-analysis/projects/CYGNSS_L1_AZ
$PY scripts/cygl1_obs_quality_by_coherency.py --da-expid DA_L1_full_xc015_fixedop --start 202001 --end 202006   # add --rebuild to redo the join
```

Outputs (gitignored), under `output/obs_quality_by_coherency/`:

- joined per-observation table: `DA_L1_full_xc015_fixedop_202001_202006_l1_qc_joined.parquet`
- all tables above: `DA_L1_full_xc015_fixedop_202001_202006_<table>.csv`
  - `coherency_deciles_{all,JFM,AMJ}`, `coherency_state`
  - `{ddm_snr,srtm_slope,sp_inc_angle}_quintiles`
  - `tail_tile_clustering`, `landcover_share_by_class`
  - `snr_tercile_edges`, `coherency_class_x_snr_tercile`, `keep_band_summary`
