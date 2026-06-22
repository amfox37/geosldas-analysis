# Legacy vs H121 high-latitude bias investigation

## 1. The headline finding

H121 reads systematically drier than Legacy above roughly 50N, by up to ~20% sat in the
60-70N band (10-day global super-ob sample, 2020-06-01..10; full breakdown and figures in
`projects/ascat_da/notebooks/legacy_vs_h121_highlat_bias.ipynb`):

| Band | Bias (H121 - Legacy) |
|---|---:|
| <30N | +1.9% sat |
| 30-60N | -11.9% sat |
| >60N | -17.7% sat |

This is confirmed in both the OFA tile/cycle stream (`innov_newQC` run) and independently
reproduced from raw BUFR/H121 files with the current default Python QC -- the two agree
closely at every latitude band, so the bias is not an artifact of the GEOSldas reader or its
cycle assignment.

## 2. What was ruled out

- **Tile/super-ob binning**: rebinning the same raw obs onto a plain 0.25-deg grid instead of
  M36 GEOS tiles gives the same latitude-bias curve to within ~1% sat. Not a tile-lookup
  artifact.
- **The new H121 QC thresholds** (`sens_min`, `subsfc_max`, `bsflag` noise bit; see
  `legacy_vs_h121_qc_flags.md`): re-running with the original, pre-investigation QC (no
  sensitivity/subsurface-scattering/backscatter-noise screening) gives essentially the same
  high-latitude bias (-16.6% vs -17.7% above 60N). The bias predates these thresholds.
- **Residual snow/frozen-ground contamination**: layering an *additional* screen on H121's own
  `snow_cover_probability`/`frozen_soil_probability` fields (>=50% reject; these are themselves
  ERA5 climatological probabilities, not real-time flags) does not reduce the bias -- it gets
  slightly worse (-17.7% to -20.1% above 60N) while cutting the sample roughly in half. The
  obs that survive still show the same bias.

## 3. What the bias correlates with

H121's own per-obs `surface_soil_moisture_sensitivity` (mean backscatter sensitivity to soil
moisture, dB) drops from ~4-6 dB through the tropics/mid-latitudes to ~2.4-2.7 dB in the
50-70N band -- consistent with degraded retrieval performance over boreal-forest/tundra
vegetation. This is suggestive but not a full explanation on its own: there's a comparably low
sensitivity dip near the equator (rainforest) where the bias is much smaller, so latitude-band
zonal averaging conflates multiple biomes and sensitivity alone doesn't cleanly predict bias
magnitude. `subsurface_scattering_probability` is not informative at high latitude --
missing/defaulted-to-0 for ~59% of obs at 50-60N and ~99% at 60-70N (H121 doesn't appear to
compute it that far north).

## 4. Which product is actually wrong? (2026-06-22 revision)

Earlier framing in this investigation assumed H121 was "the biased one" because it disagrees
with Legacy. That assumption doesn't hold up: comparing each product's own innovation (O-F)
against the same Catchment background, by latitude band (full June, `innov_newQC` run, unscaled
monitor-mode O-F):

| Lat band | Legacy mean innov | H121 mean innov |
|---|---:|---:|
| 30-40N | 17.4 | 17.5 |
| 40-50N | 20.7 | 14.7 |
| 50-60N | **40.4** | 13.5 |
| 60-70N | **49.6** | 23.0 |
| 70-80N | **44.5** | 25.3 |

Legacy's innovation roughly doubles from its 30-50N baseline (~19) to 50-80N (~40-50), a jump
of +25 to +30. H121's innovation barely moves (+10 over its own baseline, after a dip at
50-60N). **Legacy is the one that diverges from the model background at high latitude, not
H121.** If Catchment is closer to truth, H121's high-latitude values are more consistent with
it than Legacy's are -- the opposite of the naive read of the raw Legacy-vs-H121 disagreement.
(Caveat: this doesn't prove H121 is "right" -- Catchment's own snow/frozen-soil physics could
share a similar high-latitude blind spot as whichever obs product it was tuned/evaluated
against.)

Notebook: `projects/ascat_da/notebooks/legacy_vs_h121_june_newqc.ipynb`, section 14.

## 5. A separate, related finding: why Legacy obs are sparser everywhere (not just high lat)

Single-window case study (Oklahoma/S Kansas box, Metop-A, 2020-06-01 cycle 1, raw BUFR/H121,
not specific to high latitude): in this window, GEOS-tile super-obbing gave Legacy 21 tiles
vs H121's 64 in the same box. Disabling Legacy's QC entirely recovered footprints across the
whole box (n=67 -> n=184, comparable in extent to H121) -- so the gap is not missing raw data,
it's QC rejection. Isolating each Legacy filter individually showed `smpf_ok`
(`soilMoistureProcessingFlag == 0`) alone reproduces the full-QC gap almost exactly (n=68 vs
n=67); the other three filters (`smcf_ok`, `tpcx_max`, `iwfr_max`) removed almost nothing on
their own.

This matches and concretely demonstrates (at the single-pass level) what `legacy_vs_h121_qc_flags.md`
Sect. 2-3 already found from the bit-level statistics: `SMPF` is a 16-bit `PROCESSING_FLAGS`
field, and ~60% of valid-range Legacy obs fail it, dominated by bits 16/32 ("slope Mid-Fore/Mid-Aft
beam out of range"). This is structural: Legacy estimates beam slope locally per swath pass
(inherently noisy, needs a per-pass sanity check that often fails); H121 derives slope/curvature
from an offline, kernel-smoothed multi-year climatology (Hahn et al. 2026 Sect. 3.3.3) and has
no equivalent per-pass failure mode. Legacy's QC is doing what it's designed to do here --
`legacy_vs_h121_qc_flags.md` already concludes "nothing needs to change" on the Legacy QC side.

Notebook: `projects/ascat_da/notebooks/legacy_vs_h121_june_newqc.ipynb`, sections 13 and 13b.

## 6. Open questions / not yet done

- No independent (ISMN/GLDAS/in-situ) validation of which product is more accurate at high
  latitude -- everything above compares the two products and the model background to each
  other, not to ground truth.
- Have not derived H121-specific z-score scaling parameters yet (Legacy-derived scaling
  parameters are known to distort H121's spatial innovation pattern -- see
  `ascat_h121_integration_test.md`). This matters more for high latitude given the bias found
  here.
- No decision yet on how to handle the high-latitude bias before assimilating H121
  operationally (options: derive latitude-aware/H121-specific scaling, mask high latitudes,
  or proceed and monitor -- not evaluated here).
