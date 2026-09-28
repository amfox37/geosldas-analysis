# CYGNSS L1 DA, fixed operator: 2020 arm comparison and noise-vs-gain split, extended through 2022

*2026-09-27. Arizona limited domain (lon −118..−106, lat 29..40), EASEv2 M36, 2020-01-01 → 2021-01-01, fixed-operator build
(md5 6a70b496). Every arm is scored against the single unscaled open loop `OLv8_M36_AZ_fixedop`, using the DA run's scaled obs and
cross-masking. Extended 2026-09-28: four arms continued through 2022-12-31, see §3.*

## Summary

**2021–2022 update (§3):**
- **Filter+err39 is still the best L1 arm**, and still beats the OL on every monitor except ASCAT, but by about half the 2020 margin.
  Pooled 2021–22: SMOS −0.71, SMAP −0.93, ASCAT +0.61, L3 −1.51. Whole period 2020–22: SMAP −1.20, L3 −2.08, ASCAT +0.55.
- **The spring degradation recurs every year, with a moving window:** May–Jun 2020, Apr–Jun 2021, Feb–May 2022.
- **It has two mechanisms.** In 2020 and 2022 the L1 → Tb signal collapses (corr < 0.1, α_opt ≈ 0). In 2021 the increments overshoot
  (corr normal, α_opt ≈ 0.45–0.75). A spring R increase over a wide window (Feb–Jun) is the next experiment. A fixed May–Aug gate would
  have missed 2022.
- **L3-only is still better than every L1 arm on every monitor; SMAP-only's ASCAT damage grew** (+12.7% in 2021–22).

**2020 results (§1–2):**

- **Best L1 configuration: full-stream L1, xcorr 0.15, coherency filter 0.40–2.16 with a rebuilt z-score climatology, errstd 3.9**
  (`DA_L1_full_xc015_coh040216_err39_fixedop`). Full-year monitor O-F std change vs OL: SMOS −1.35%, SMAP −1.52%, ASCAT +0.23%, CYGNSS L3
  −2.99%, and L1's own O-F −0.71%. It beats the OL on SMAP in 10 of 12 months and on L3 in 11 of 12 months, and has the smallest ASCAT
  penalty of the full-stream L1 arms.
- **The coherency filter provides most of the gain; errstd 3.9 adds robustness.** The filter alone reaches SMOS −1.29 / SMAP −1.51 /
  L3 −3.01 but ASCAT +1.99. Adding errstd 3.9 halves the ASCAT damage and damps the spring and fall degradations, at the cost of some
  Jan–Feb gain.
- **The best L1 arm is still behind the single-sensor baselines.** It comes close to L3-only on SMOS/SMAP Tb (−1.35/−1.52 vs −1.67/−1.75)
  but is well behind on L3 (−3.0 vs −7.6). SMAP-only is in a different class (Tb −12 to −14%) but damages ASCAT (+8.3%) and L3 (+1.8%).
- **May–June is bad in every L1 arm, and the noise-vs-gain split shows why.** In May–June the correlation between the forecast change and
  the OL innovation, corr(dF, O−F_OL), drops to 0.02–0.05 for SMAP/SMOS Tb in all L1 arms, and α_opt ≈ 0. The L1 increments carry
  essentially no information about Tb in those months. This is a **signal problem, not an R-amplitude problem**, so no global R value
  can fix it. L3-only keeps corr ≈ 0.15 in those months. L1's own α_opt stays ≈ 0.6–0.7 in May–June, so the L1 updates still fit the
  L1 obs but don't carry over to Tb.
- **The filter works by raising corr(dF, innov),** e.g. SMAP Jan .29 → .38 and Oct .22 → .35, not just by removing noise.
- **errstd 3.9 roughly halves the noise but gives up less gain,** so α_opt rises. In Jan–Apr the filter+err39 arm has SMAP α_opt of
  1.6–2.3, meaning its winter increments are now **too small**.
- **Implication: use a seasonal R, or seasonal gating, instead of one global R.** Use a smaller R in winter and fall (errstd ≈ 2.75 or
  lower with the filter), and a very large R or no L1 assimilation in about May–August. Separately, investigate why the L1 → Tb signal
  decouples in late spring (vegetation green-up, roughness, operator, or scaling).

## Arms

| tag | EXP_ID | L1 stream | xcorr | errstd | notes |
|---|---|---|---|---|---|
| full_xc015_fixedop | DA_L1_full_xc015_fixedop | full | 0.15 | 2.75 | **benchmark** |
| full_xc015_err39_fixedop | DA_L1_full_xc015_err39_fixedop | full | 0.15 | 3.9 | benchmark + larger R |
| full_xc015_coh040216_fixedop | DA_L1_full_xc015_coh040216_fixedop | full, coherency 0.40–2.16 | 0.15 | 2.75 | + rebuilt L1 clim on the screened population |
| full_xc015_coh040216_err39_fixedop | DA_L1_full_xc015_coh040216_err39_fixedop | full, coherency 0.40–2.16 | 0.15 | 3.9 | filter + clim + larger R |
| dense075_coh05_fixedop | DA_L1_dense075_coh05_fixedop | thinned | 0.625 | 2.75 | pre-fix "best" config rerun |
| smap_fixedop | DA_SMAP_fixedop | — (monitored) | — | — | SMAP Tb assimilated only |
| l3_fixedop | DA_L3_fixedop | — (monitored) | — | — | CYGNSS L3 SM assimilated only |

All 7 runs finished: cap_restart = 20210101, 0 `LDAS ERROR`/`forrtl` in the log, and the full set of ens_avg ObsFcstAna files in every
month (247/232/248/240/248/240/248/248/240/248/240/248).

## 1. Monitor O-F std change by month

Each value is the species-mean `% (DA−OL)/OL` of O-F stdv; negative = better than OL. SMOS and SMAP each average 4 species (H/V × A/D),
ASCAT averages 3 (A/B/C), and L3 is CYGNSS_SM_6hr. L1 is the arm's own CygL1 O-F change. The "2020" column is one pooled full-year
statistic, **not** the mean of the months: it also rewards reductions in month-to-month bias, which is why, for example, the benchmark's
full-year L3 is −1.65 even though most of its individual months are positive.

### SMOS

| arm | 202001 | 202002 | 202003 | 202004 | 202005 | 202006 | 202007 | 202008 | 202009 | 202010 | 202011 | 202012 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -5.48 | -1.81 | +0.18 | +1.86 | +2.66 | +1.62 | +0.74 | +2.71 | -0.47 | +0.18 | +1.55 | -0.13 | +0.12 |
| full_xc015_err39_fixedop | -4.67 | -1.24 | -0.03 | +0.89 | +1.12 | +0.30 | +0.19 | +1.54 | -1.41 | -1.58 | -0.01 | -0.48 | -0.42 |
| full_xc015_coh040216_fixedop | -7.29 | -2.49 | -0.58 | -1.27 | +1.40 | +1.47 | -0.12 | +0.95 | -1.85 | -4.57 | -0.07 | -1.39 | -1.29 |
| full_xc015_coh040216_err39_fixedop | -5.62 | -1.88 | -0.44 | -1.11 | +0.39 | +0.52 | -0.24 | -0.18 | -2.07 | -5.09 | -1.01 | -1.25 | -1.35 |
| dense075_coh05_fixedop | -0.41 | -0.66 | -0.18 | +1.02 | +1.05 | +1.23 | +0.08 | -0.28 | -1.19 | -1.27 | +1.12 | +0.51 | -0.06 |
| smap_fixedop | -22.93 | -15.96 | -10.41 | -14.73 | -3.44 | -7.53 | -10.19 | -13.98 | -21.61 | -20.83 | -18.49 | -15.01 | -14.41 |
| l3_fixedop | -5.99 | -2.11 | -0.38 | -3.67 | -1.15 | -0.49 | -0.64 | -1.65 | -3.13 | -5.54 | -3.55 | -0.74 | -1.67 |

### SMAP

| arm | 202001 | 202002 | 202003 | 202004 | 202005 | 202006 | 202007 | 202008 | 202009 | 202010 | 202011 | 202012 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -4.46 | -0.90 | -0.43 | +0.79 | +4.24 | +3.16 | +0.41 | +1.51 | -0.47 | +1.02 | +2.98 | -0.16 | +0.17 |
| full_xc015_err39_fixedop | -3.84 | -0.61 | -0.42 | -0.02 | +2.22 | +1.23 | -0.00 | +0.09 | -1.27 | -1.90 | +0.50 | -0.68 | -0.49 |
| full_xc015_coh040216_fixedop | -6.29 | -2.10 | -1.94 | -1.51 | +3.26 | +2.64 | -0.21 | -0.48 | -2.04 | -6.10 | +0.04 | -1.49 | -1.51 |
| full_xc015_coh040216_err39_fixedop | -4.88 | -1.48 | -1.44 | -1.64 | +1.88 | +1.24 | -0.36 | -0.69 | -2.10 | -6.77 | -1.50 | -1.28 | -1.52 |
| dense075_coh05_fixedop | -0.37 | -0.46 | -0.15 | +1.09 | +2.37 | +2.05 | -0.32 | -0.21 | -0.39 | -0.88 | +0.80 | +0.11 | +0.09 |
| smap_fixedop | -20.83 | -13.94 | -10.00 | -11.54 | -0.41 | -5.18 | -6.55 | -6.95 | -15.70 | -23.50 | -16.51 | -17.54 | -11.81 |
| l3_fixedop | -6.25 | -2.80 | -1.17 | -3.87 | -0.82 | -1.08 | -0.24 | -1.29 | -2.73 | -7.48 | -3.94 | -1.45 | -1.75 |

### ASCAT

| arm | 202001 | 202002 | 202003 | 202004 | 202005 | 202006 | 202007 | 202008 | 202009 | 202010 | 202011 | 202012 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -0.73 | +0.38 | +4.54 | +6.58 | +5.08 | +2.12 | +2.85 | +3.61 | +6.87 | +13.16 | +8.22 | +3.69 | +3.23 |
| full_xc015_err39_fixedop | -1.68 | -0.44 | +1.94 | +3.64 | +2.46 | +1.03 | +0.52 | +0.95 | +3.20 | +6.53 | +3.21 | +0.82 | +0.79 |
| full_xc015_coh040216_fixedop | -1.89 | -0.18 | +3.55 | +5.74 | +2.59 | -0.79 | +1.26 | +2.71 | +6.28 | +8.72 | +5.63 | +2.72 | +1.99 |
| full_xc015_coh040216_err39_fixedop | -2.36 | -0.72 | +1.52 | +3.15 | +1.16 | -0.80 | +0.11 | +0.77 | +2.82 | +4.25 | +1.87 | +0.79 | +0.23 |
| dense075_coh05_fixedop | +1.09 | +0.33 | +0.92 | +0.74 | +1.28 | -0.45 | +0.65 | +1.28 | +1.89 | +2.75 | +1.46 | +0.16 | +0.78 |
| smap_fixedop | -7.05 | -1.47 | +4.72 | +7.36 | +5.71 | +10.26 | +8.00 | +5.62 | +5.90 | +9.99 | +17.56 | +15.26 | +8.29 |
| l3_fixedop | -4.23 | -2.80 | -0.21 | +0.80 | -1.23 | -3.01 | -3.67 | -2.12 | -1.56 | -0.20 | -2.38 | -1.95 | -4.89 |

### L3

| arm | 202001 | 202002 | 202003 | 202004 | 202005 | 202006 | 202007 | 202008 | 202009 | 202010 | 202011 | 202012 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -1.17 | +0.45 | +0.60 | +1.71 | +3.24 | +0.83 | -0.09 | -0.12 | +1.09 | +1.69 | -0.78 | +1.62 | -1.65 |
| full_xc015_err39_fixedop | -1.47 | +0.01 | -0.02 | +0.29 | +1.11 | -0.81 | -1.36 | -1.09 | -0.49 | -0.26 | -2.26 | +0.09 | -2.29 |
| full_xc015_coh040216_fixedop | -3.18 | -1.27 | -0.96 | +0.48 | +1.99 | +0.08 | -1.21 | -1.77 | -1.24 | +0.69 | -4.04 | +0.74 | -3.01 |
| full_xc015_coh040216_err39_fixedop | -2.86 | -1.18 | -1.09 | -0.25 | +0.58 | -0.87 | -1.51 | -1.72 | -2.01 | -0.68 | -4.34 | -0.36 | -2.99 |
| dense075_coh05_fixedop | +0.14 | -0.40 | -0.27 | +0.72 | +1.95 | +1.05 | +0.09 | -0.08 | -0.05 | +0.84 | -0.54 | +0.33 | -0.34 |
| smap_fixedop | -7.20 | -1.72 | +2.54 | +2.14 | +6.34 | +6.94 | +0.83 | +2.92 | -1.57 | -1.02 | -4.04 | +0.84 | +1.80 |
| l3_fixedop | -5.31 | -3.26 | -1.63 | -3.08 | -1.87 | -4.46 | -3.35 | -3.20 | -3.86 | -2.85 | -6.67 | -3.64 | -7.61 |

### L1

| arm | 202001 | 202002 | 202003 | 202004 | 202005 | 202006 | 202007 | 202008 | 202009 | 202010 | 202011 | 202012 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -0.74 | -0.27 | -0.55 | -0.06 | -0.16 | -0.15 | -0.14 | -0.18 | -0.75 | +0.16 | -0.34 | -0.43 | -0.36 |
| full_xc015_err39_fixedop | -0.60 | -0.25 | -0.39 | -0.08 | -0.15 | -0.19 | -0.11 | -0.22 | -0.67 | -0.00 | -0.41 | -0.54 | -0.34 |
| full_xc015_coh040216_fixedop | -0.95 | -0.59 | -1.06 | -0.38 | -0.03 | -0.12 | -0.28 | -0.69 | -0.89 | -0.77 | -0.90 | -0.47 | -0.81 |
| full_xc015_coh040216_err39_fixedop | -0.71 | -0.47 | -0.70 | -0.32 | -0.07 | -0.17 | -0.29 | -0.62 | -0.82 | -0.83 | -0.99 | -0.46 | -0.71 |
| dense075_coh05_fixedop | -0.02 | -0.26 | -0.41 | -0.08 | -0.10 | -0.05 | -0.14 | -0.27 | -0.34 | -0.26 | -0.58 | -0.78 | -0.26 |
| smap_fixedop | -1.64 | -0.26 | -2.03 | -0.56 | +0.14 | +0.00 | -0.28 | -0.52 | -1.23 | -0.38 | -0.38 | -0.58 | -0.67 |
| l3_fixedop | -0.55 | -0.57 | -0.50 | -0.40 | -0.35 | -0.38 | -0.24 | -0.47 | -0.44 | -0.21 | -0.82 | -0.26 | -0.51 |


## 2. Noise-vs-gain split

For every obs present in both the DA and OL ens_avg ObsFcstAna files, matched on (time, species, tile, lon, lat) and using the DA run's
scaled obs O:

    inn = O − F_OL,  dF = F_DA − F_OL
    dMSE  = Σ(dF² − 2 dF·inn) / Σ inn²       (% of OL MSE; includes the bias, unlike the stdv table above)
    noise = Σ dF² / Σ inn²                   (cost of moving the forecast at all)
    gain  = −2 Σ dF·inn / Σ inn²             (benefit of moving it toward the obs)
    alpha_opt = Σ dF·inn / Σ dF²             (increment scale factor that would minimise MSE)
    corr  = corr(dF, inn)

- α_opt > 1: increments are too small.
- α_opt < 1: increments are too big.
- α_opt ≈ 0: the increments carry no information for that monitor.

Values are computed per month and for the pooled year.

### Key rows (SMAP Tb; SMOS is the same pattern)

| | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec | yr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| corr, benchmark | .29 | .14 | .12 | .14 | .09 | .08 | .12 | .10 | .19 | .22 | .16 | .15 | .15 |
| corr, filter | .38 | .20 | .20 | .19 | .03 | .05 | .11 | .15 | .22 | .34 | .19 | .17 | .19 |
| corr, filter+err39 | .36 | .19 | .18 | .17 | .02 | .04 | .10 | .14 | .21 | .35 | .18 | .15 | .17 |
| corr, L3-only | .45 | .28 | .16 | .32 | .16 | .15 | .09 | .14 | .27 | .37 | .27 | .16 | .19 |
| α_opt, benchmark | 1.12 | .82 | .53 | .82 | .28 | .36 | .36 | .26 | .57 | .31 | .31 | .39 | .48 |
| α_opt, filter+err39 | 2.26 | 1.80 | 1.57 | 1.58 | −.02 | .19 | .34 | .34 | .69 | .71 | .46 | .45 | .80 |
| noise, benchmark | 7.9 | 4.2 | 3.7 | 8.7 | 10.6 | 8.9 | 7.1 | 9.2 | 10.9 | 26.3 | 20.2 | 9.9 | 8.8 |
| noise, filter+err39 | 2.9 | 1.4 | 1.2 | 2.5 | 2.9 | 2.9 | 1.9 | 3.2 | 5.0 | 12.0 | 9.2 | 4.5 | 3.2 |

### Interpretation

1. **The coherency filter improves increment direction.** corr(dF, inn) rises for SMAP/SMOS/L3/L1 in nearly every month. The filter
   removes observations that pushed the state the wrong way, not just noisy ones.
2. **errstd 3.9 cuts noise about 2× but gain only about 1.5×, so the net effect is positive.** α_opt rises, and in Jan–Apr it overshoots
   to above 1.5 for Tb in the filtered arm: the winter increments are now too small.
3. **The late-spring failure is a signal failure.**
   - May–June Tb corr is 0.02–0.05 in every L1 arm, including the thinned one (≈0.00), and α_opt ≈ 0.
   - L1's own α_opt stays 0.6–0.7 in those months, so the updates fit the L1 obs but don't carry over to Tb.
   - L3-only keeps corr ≈ 0.15 in the same months.
   - Therefore scaling R cannot fix May–June. It needs seasonal gating or a physical explanation of why the L1 → soil-moisture → Tb
     chain breaks (vegetation green-up / VWC, roughness, operator, or scaling-clim season).
4. **October is the L1 arms' best month for Tb.** Filter / filter+err39 corr is .34/.35, matching L3-only's .37, and the SMAP O-F std
   change is −6.1/−6.8%. The noise is also highest in October–November (dry-down or monsoon-tail variability), so the filter matters most
   there.
5. **For L3 as a monitor, the L1 arms keep corr at 0.08–0.20 all year.** With filter+err39, α_opt is ≥ ~1 in most months, which is
   consistent with its L3 score being better than the OL in 11 of 12 months.

## 3. 2021–2022 extension: do the 2020 findings hold?

*Added 2026-09-28.* The four arms that decide the next experiment (coherency filter, filter + errstd 3.9, L3-only, SMAP-only) were
continued unchanged from their 2021-01-01 restarts to 2023-01-01 (CAP.rc END_DATE edit, same build, same exeinp; cap_restart 20230101
and 0 `LDAS ERROR`/`forrtl` in every run). The benchmark, benchmark + errstd 3.9, and thinned arms stop at 2020. Scoring is identical
to §1: same OL, the DA run's scaled obs, cross-masked, `% (DA−OL)/OL` of O-F stdv, negative = better.

### Per-year and pooled

| arm | period | SMOS | SMAP | ASCAT | L3 | L1 own |
|---|---|---:|---:|---:|---:|---:|
| full_xc015_coh040216_fixedop | 2020 | −1.29 | −1.51 | +1.99 | −3.01 | −0.81 |
| | 2021 | −0.42 | −0.68 | +1.81 | −1.35 | −0.64 |
| | 2022 | −1.11 | −1.47 | +2.22 | −0.75 | −0.79 |
| | 2021–22 | −0.77 | −1.02 | +2.10 | −1.15 | −0.71 |
| **full_xc015_coh040216_err39_fixedop** | 2020 | −1.35 | −1.52 | +0.23 | −2.99 | −0.71 |
| | 2021 | −0.55 | −0.70 | +0.35 | −1.62 | −0.62 |
| | 2022 | −0.88 | −1.21 | +0.77 | −1.23 | −0.65 |
| | 2021–22 | **−0.71** | **−0.93** | **+0.61** | **−1.51** | −0.64 |
| l3_fixedop | 2020 | −1.67 | −1.75 | −4.89 | −7.61 | −0.51 |
| | 2021 | −0.96 | −1.16 | −2.06 | −4.02 | −0.43 |
| | 2022 | −1.04 | −1.35 | −3.10 | −5.65 | −0.68 |
| | 2021–22 | −0.87 | −1.20 | −2.57 | −5.09 | −0.55 |
| smap_fixedop | 2020 | −14.41 | −11.81 | +8.29 | +1.80 | −0.67 |
| | 2021 | −14.68 | −9.09 | +13.05 | +3.04 | −0.53 |
| | 2022 | −11.58 | −11.79 | +10.54 | +4.73 | −0.82 |
| | 2021–22 | −13.02 | −10.26 | +12.70 | +3.88 | −0.65 |

Each year or two-year value is one pooled statistic, not the mean of the months (see §1). The whole-period 2020–2022 pooled values are
in `output/full_period_2020_2022/` (see Reproduction).

### Findings

1. **The best L1 arm still beats the OL on every monitor except ASCAT, but by about half as much as in 2020.** Filter+err39 goes from
   SMAP −1.52 / L3 −2.99 in 2020 to −0.93 / −1.51 pooled over 2021–22. Part of the 2020 gain looks like a one-off: Jan 2020 (SMOS/SMAP −5
   to −6%) was partly spin-up from the OL restart. 2021 is the weakest year; 2022 recovers on Tb (SMAP −1.21) but not on L3.
2. **Filter+err39 is still the L1 configuration to carry forward.** Against the filter alone it is equal on Tb (SMAP −0.93 vs −1.02),
   better on L3 (−1.51 vs −1.15), and cuts the ASCAT penalty from +2.1 to +0.6. It beats the OL on SMAP in 9 of 12 months in 2021 and 7 of
   12 in 2022, and on L3 in 11 of 12 and 8 of 12.
3. **The spring degradation is not a 2020 accident; it recurs every year, and its window moves.** Filter-arm SMAP months worse than
   the OL: May–Jun 2020, Apr–Jun 2021 (+3.1/+1.9/+0.4), Feb–May 2022 (+2.5/+0.7/+2.7/+3.4). L3 follows the same pattern (May 2022 +4.6).
   errstd 3.9 roughly halves each of these (Apr 2021 SMAP +3.1 → +1.0; May 2022 +3.4 → +1.3) but does not remove them. A second, weaker
   bad month appears in August (Aug 2022 SMAP +1.6 / +1.3). A fixed May–August gate would have missed Feb–Apr 2022.
4. **Oct–Jan is consistently the L1 arms' best period** (e.g. Jan 2022 SMAP −5.0 / −3.9, Oct and Dec 2022 −2.1 / −2.9 in the filter arm),
   as in 2020 (Oct 2020 −6.1 / −6.8).
5. **L3-only is still better than the best L1 arm on every independent monitor, in every year**, though its margin on Tb is small
   (2021–22 SMAP −1.20 vs −0.93) and its own spring is also weaker (May, Jun and Aug 2022 are its only SMAP months worse than the OL, by
   +0.1 to +0.6).
6. **SMAP-only's asymmetry got worse:** Tb −10 to −13%, but ASCAT +12.7% (2020: +8.3%) and L3 +3.9% (2020: +1.8%). The ASCAT damage
   peaks in autumn (Sep–Dec 2021 +16 to +26%).
7. **Assimilating L1 barely improves the fit to L1 itself more than assimilating another sensor does.** L1 own O-F is −0.64% for
   filter+err39, vs −0.55% for L3-only and −0.65% for SMAP-only (2021–22). Most of the L1 O-F reduction comes from a better soil-moisture
   state in general, not from L1-specific information.
8. **The spring failure has two different mechanisms, depending on the year** (noise-vs-gain split, SMAP Tb, same method as §2):

   | | 2020 May–Jun | 2021 Apr–May | 2022 Feb–May |
   |---|---|---|---|
   | corr(dF, inn), filter / filter+err39 | .03–.05 / .02–.04 | .13–.16 / .14–.16 | .08–.15 / .07–.15 |
   | α_opt, filter+err39 | −.02 / .19 | .62 / .74 | −.17 / .28 / .30 / .40 |
   | noise, filter+err39 (% of OL MSE) | 2.9 / 2.9 | 9.7 / 8.7 | 5.9 / 3.1 / 11.9 / 8.4 |

   - **2020 and 2022: signal collapse.** corr drops below ~0.1 and α_opt is near 0, as in §2.
   - **2021: overshoot.** corr stays at a normal 0.13–0.16, but the noise is 3–5× the winter value and α_opt ≈ 0.45–0.75, so the
     increments are about 1.5–2× too big. That is an R-amplitude problem, and a larger spring R would fix it.
   - Either way, the right move in spring is to strongly down-weight L1, and the window cannot be a fixed calendar block (it started in
     February in 2022).
   - **Winter increments are still too small.** α_opt for filter+err39 is 1.9–2.0 in Jan 2021, Jan 2022, Oct 2022 and Dec 2022 (and
     1.9–2.0 in July of both years, at monsoon onset). A smaller winter R would be worth more there.
   - Whole-period (36-month) pooled corr is 0.16–0.17 for both L1 arms and for L3-only, and 0.48 for SMAP-only.

### Whole period 2020–2022 (pooled over 36 months)

| arm | SMOS | SMAP | ASCAT | L3 | L1 own |
|---|---:|---:|---:|---:|---:|
| full_xc015_coh040216_fixedop | −1.02 | −1.29 | +2.18 | −1.86 | −0.76 |
| **full_xc015_coh040216_err39_fixedop** | **−0.97** | **−1.20** | **+0.55** | **−2.08** | −0.68 |
| l3_fixedop | −1.14 | −1.48 | −3.61 | −6.01 | −0.54 |
| smap_fixedop | −14.12 | −11.48 | +12.59 | +2.94 | −0.67 |

These are the species-group means of the per-species pooled % change, computed from the whole-period pkl files
(`spatial_stats_*_202001_202212.pkl`, N-weighted pooling of the 36 monthly rows, same formula as the bundle README).

### Monthly detail, 2021–2022

#### SMOS

| arm | year | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec | year |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_coh040216_fixedop | 2021 | -0.80 | -1.67 | +0.43 | +0.89 | +1.09 | +0.78 | -0.71 | +0.51 | -0.69 | -0.49 | -0.64 | -0.67 | -0.42 |
| full_xc015_coh040216_fixedop | 2022 | -3.65 | +0.73 | +1.30 | +2.04 | +4.87 | -0.10 | -1.49 | -0.42 | -0.40 | -2.02 | -0.26 | -3.46 | -1.11 |
| full_xc015_coh040216_err39_fixedop | 2021 | -0.56 | -1.86 | -0.06 | -0.28 | +0.42 | +0.11 | -0.61 | +0.20 | -0.72 | -0.62 | -1.07 | -0.62 | -0.55 |
| full_xc015_coh040216_err39_fixedop | 2022 | -2.76 | +0.19 | +0.28 | +0.34 | +2.27 | -0.05 | -1.17 | -0.25 | -0.38 | -1.69 | -0.43 | -2.52 | -0.88 |
| l3_fixedop | 2021 | -1.12 | -4.12 | -1.14 | -1.48 | -1.39 | -1.90 | -0.31 | -0.59 | -0.53 | +0.40 | -1.69 | -0.68 | -0.96 |
| l3_fixedop | 2022 | -4.20 | -0.87 | -0.49 | -1.39 | +0.58 | +0.48 | -0.58 | -0.76 | -1.14 | -0.24 | -5.02 | -3.80 | -1.04 |
| smap_fixedop | 2021 | -20.95 | -14.09 | -12.19 | -10.35 | -7.63 | -13.04 | -12.17 | -15.11 | -10.63 | -19.95 | -18.66 | -16.86 | -14.68 |
| smap_fixedop | 2022 | -19.35 | -6.55 | -8.65 | -8.96 | -14.40 | -8.34 | -16.12 | -1.23 | -0.23 | -22.43 | -12.98 | -25.48 | -11.58 |

#### SMAP

| arm | year | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec | year |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_coh040216_fixedop | 2021 | -1.68 | -1.61 | -0.33 | +3.09 | +1.94 | +0.42 | -0.63 | +0.27 | -0.59 | -1.55 | -1.84 | -0.35 | -0.68 |
| full_xc015_coh040216_fixedop | 2022 | -5.04 | +2.52 | +0.74 | +2.66 | +3.37 | +0.02 | -1.65 | +1.57 | -0.46 | -2.10 | -0.58 | -2.85 | -1.47 |
| full_xc015_coh040216_err39_fixedop | 2021 | -1.30 | -1.66 | -0.76 | +0.99 | +0.20 | -0.42 | -0.57 | +0.15 | -0.52 | -1.25 | -1.80 | -0.27 | -0.70 |
| full_xc015_coh040216_err39_fixedop | 2022 | -3.87 | +1.53 | +0.21 | +0.31 | +1.32 | -0.11 | -1.41 | +1.29 | -0.43 | -1.62 | -0.90 | -2.15 | -1.21 |
| l3_fixedop | 2021 | -1.56 | -5.13 | -2.00 | -1.56 | -1.70 | -0.80 | -0.32 | -0.69 | -0.92 | -0.57 | -1.93 | -0.45 | -1.16 |
| l3_fixedop | 2022 | -4.82 | -1.44 | -1.43 | -3.02 | +0.28 | +0.64 | -0.52 | +0.09 | -0.58 | -1.04 | -6.43 | -3.72 | -1.35 |
| smap_fixedop | 2021 | -17.72 | -16.95 | -11.28 | -10.45 | -8.19 | -4.86 | -4.19 | -5.04 | -5.32 | -14.80 | -15.67 | -11.29 | -9.09 |
| smap_fixedop | 2022 | -18.21 | -5.67 | -8.36 | -6.55 | -13.92 | -6.21 | -5.77 | -5.05 | -10.31 | -13.20 | -13.26 | -18.19 | -11.79 |

#### ASCAT

| arm | year | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec | year |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_coh040216_fixedop | 2021 | -0.23 | +3.63 | +1.42 | +5.77 | +2.17 | +1.12 | +1.01 | +1.59 | +3.19 | +6.37 | +4.11 | +5.12 | +1.81 |
| full_xc015_coh040216_fixedop | 2022 | +0.46 | +7.08 | +1.69 | +5.12 | +2.00 | +2.56 | +3.01 | -0.14 | +2.76 | +2.67 | +2.38 | +1.93 | +2.23 |
| full_xc015_coh040216_err39_fixedop | 2021 | -1.26 | +1.07 | +0.35 | +2.93 | +0.89 | +0.38 | +0.33 | +0.37 | +1.28 | +2.98 | +1.65 | +2.92 | +0.35 |
| full_xc015_coh040216_err39_fixedop | 2022 | -0.35 | +3.77 | +0.76 | +2.35 | +0.64 | +1.25 | +1.33 | -0.51 | +0.72 | +0.93 | +0.90 | +0.48 | +0.77 |
| l3_fixedop | 2021 | -3.44 | -2.48 | -1.05 | +0.37 | -0.24 | -1.18 | -0.98 | -2.06 | -0.56 | +1.09 | -1.06 | -0.47 | -2.06 |
| l3_fixedop | 2022 | -2.92 | -0.29 | -1.02 | -0.98 | +0.41 | -1.17 | -1.01 | -1.10 | -1.93 | -1.93 | -3.57 | -1.63 | -3.10 |
| smap_fixedop | 2021 | +10.79 | +4.81 | +4.72 | +4.13 | +4.40 | +2.70 | +6.55 | +11.29 | +20.59 | +25.57 | +24.55 | +15.96 | +13.05 |
| smap_fixedop | 2022 | +6.11 | +11.07 | +8.14 | +7.43 | -2.10 | +9.80 | +12.11 | +0.20 | +2.62 | +13.38 | +15.07 | +16.23 | +10.54 |

#### L3

| arm | year | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec | year |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_coh040216_fixedop | 2021 | -0.94 | -2.55 | -3.01 | +0.81 | +1.20 | -0.35 | -0.87 | +0.08 | -0.95 | +0.75 | -1.87 | +0.81 | -1.35 |
| full_xc015_coh040216_fixedop | 2022 | -3.55 | +1.79 | -0.05 | +2.38 | +4.60 | -0.45 | -0.01 | +0.48 | -0.61 | -1.44 | +0.68 | -1.72 | -0.75 |
| full_xc015_coh040216_err39_fixedop | 2021 | -1.54 | -2.41 | -2.91 | -0.01 | -0.08 | -1.56 | -0.73 | -0.22 | -1.21 | -0.17 | -2.43 | +0.06 | -1.62 |
| full_xc015_coh040216_err39_fixedop | 2022 | -3.22 | +0.81 | -0.29 | +0.89 | +2.50 | -1.12 | -0.73 | +0.16 | -1.25 | -1.80 | -0.31 | -1.94 | -1.23 |
| l3_fixedop | 2021 | -4.80 | -4.22 | -3.78 | -1.27 | -1.82 | -5.35 | -1.01 | -4.22 | -3.75 | -1.57 | -4.44 | -1.99 | -4.02 |
| l3_fixedop | 2022 | -5.38 | +0.80 | -1.50 | -1.37 | -0.31 | -1.74 | -2.43 | -3.28 | -6.47 | -5.58 | -7.53 | -4.71 | -5.65 |
| smap_fixedop | 2021 | +5.77 | -1.91 | -2.16 | -0.83 | -1.32 | +0.15 | +2.33 | +5.52 | +7.72 | +10.96 | -0.89 | -0.71 | +3.04 |
| smap_fixedop | 2022 | -6.51 | +0.76 | +0.52 | -0.11 | +0.61 | +7.90 | +6.61 | +2.47 | +1.28 | +10.52 | -2.54 | +6.27 | +4.73 |

#### L1

| arm | year | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec | year |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_coh040216_fixedop | 2021 | -0.20 | -0.62 | -1.21 | -0.74 | -0.25 | +0.18 | -0.72 | -0.19 | -0.64 | -0.47 | -0.74 | -0.02 | -0.64 |
| full_xc015_coh040216_fixedop | 2022 | -1.10 | +0.01 | -0.08 | -0.04 | +0.11 | -0.30 | -0.61 | -0.39 | -0.57 | -1.84 | -0.11 | -1.70 | -0.79 |
| full_xc015_coh040216_err39_fixedop | 2021 | -0.29 | -0.58 | -1.10 | -0.75 | -0.39 | -0.07 | -0.56 | -0.17 | -0.51 | -0.44 | -0.70 | -0.10 | -0.62 |
| full_xc015_coh040216_err39_fixedop | 2022 | -0.91 | -0.06 | -0.17 | -0.19 | -0.00 | -0.31 | -0.52 | -0.29 | -0.48 | -1.35 | -0.20 | -1.32 | -0.65 |
| l3_fixedop | 2021 | -0.67 | -0.76 | -1.05 | -0.24 | -0.46 | -0.57 | -0.15 | -0.37 | -0.57 | -0.15 | -0.81 | -0.23 | -0.43 |
| l3_fixedop | 2022 | -1.26 | -0.48 | -0.30 | -0.16 | -0.27 | -0.55 | -0.53 | -0.43 | -0.82 | -0.65 | -0.95 | -1.45 | -0.68 |
| smap_fixedop | 2021 | -1.36 | -0.54 | -1.08 | +0.04 | -0.52 | -0.47 | -0.40 | -0.39 | -1.30 | -0.11 | -0.73 | -1.20 | -0.53 |
| smap_fixedop | 2022 | -2.73 | -0.13 | +0.12 | -0.14 | +0.17 | -0.07 | -1.46 | -0.55 | -0.46 | -2.03 | -0.82 | -3.31 | -0.82 |

## 4. Recommended next steps

*Revised 2026-09-28 after the 2021–2022 extension.*

1. **Seasonal / adaptive L1 weighting experiment**, as one intervention against filter+err39. The 2021–22 results rule out a fixed
   May–August gate (2022's bad window was Feb–May, and 2021's was an overshoot, not a signal loss). Options, in order of simplicity:
   - a. Larger R (errstd ×2, i.e. 7.8) from February through June, errstd 3.9 otherwise. It would cover the bad windows of all three
     years, at the cost of some Feb–Mar gain in good years.
   - b. The same plus a smaller winter R (errstd ≈ 2.75) in Oct–Jan, where α_opt is 1.5–2. That is a second knob, so run it as its own arm
     after (a).
   - c. Adaptive: scale R per month (or per pentad) from an innovation-based statistic, e.g. a rolling corr(L1 innovation, Tb innovation)
     or Desroziers ratio. More work; only if (a)/(b) show the season is the right axis.
2. **Diagnose the L1 → Tb decoupling in 2020/2022 spring:** corr(dF, inn) stratified by vegetation / NDVI / VWC and by land cover, plus L1
   innovation vs Tb innovation correlation by month (from the obs-quality pipeline). 2021 is the useful control: a spring where the
   signal did not collapse.
3. Extend the obs-quality report (`cygl1_obs_quality_vs_qc_report.md`) to 2020–2022 using the filtered arm.
4. Still open: along-track superobbing / nugget options (`cygl1_obs_error_correlation_report.md`), why the realized L1-space HPH is
   only 30–50% of the ensemble's, and why SMAP-only degrades ASCAT so much (worse in 2021–22, peaking in autumn).

## Reproduction

- Month × arm table: `scripts/postproc_drivers/build_month_arm_table.py --start 202001 --end 202012 --arms <EXPID:tag ...>
  --log-dir output/month_arm_table_2020 --reuse-dirs output/overnight_20260927 --out-prefix output/month_arm_table_2020/table`.
  It wraps `score_cygl1_arm.py` and reuses existing score logs and cached monthly sums; a full rebuild takes under a minute.
  The output is `output/month_arm_table_2020/table.{md,csv}`.
- Noise-vs-gain split: `scripts/cygl1_noise_gain_split.py --da-expid <EXP> --start 202001 --end 202012 --csv ...`. It was run as SLURM
  array `output/noise_gain_split/run_ngs_2020.sh` (job 58598164, about 3 min per arm). The output is
  `output/noise_gain_split/<tag>_202001_202012.csv`.
- Overnight scoring log (individual per-month score logs): `output/overnight_20260927/`.
- 2021–2022 month × arm tables (§3): `output/month_arm_table_2021_2022/score_arm.sh` (one SLURM job per arm, jobs 58620093–96, 5–7 min
  each; `--time=1:00:00` is plenty). The outputs are `table_<SHORT>_{months,2021,2022}.{md,csv}` and the combined
  `month_arm_omf_stdv_pct_vs_OL_2021_2022.csv`.
- Whole-period 2020–2022 stats and noise-vs-gain split (§3): SLURM array `output/full_period_2020_2022/run_full_period.sh` (job 58621001,
  about 10 min per arm). It runs `score_cygl1_arm.py --start 20200101 --cap-end 20230101`, which writes
  `spatial_stats_*_202001_202212.pkl` and `temporal_stats_*_20200101_20221231.nc4` to `output/postproc_paired_density/stats_output/`,
  and writes `noise_gain_<tag>_202001_202212.csv`.

## Appendix: full noise-vs-gain tables (all arms, all monitors)

#### SMAP — corr

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 0.29 | 0.14 | 0.12 | 0.14 | 0.09 | 0.08 | 0.12 | 0.10 | 0.19 | 0.22 | 0.16 | 0.15 | 0.15 |
| full_xc015_err39_fixedop | 0.28 | 0.11 | 0.10 | 0.13 | 0.08 | 0.08 | 0.10 | 0.10 | 0.18 | 0.23 | 0.15 | 0.13 | 0.14 |
| full_xc015_coh040216_fixedop | 0.38 | 0.20 | 0.20 | 0.19 | 0.03 | 0.05 | 0.11 | 0.15 | 0.22 | 0.34 | 0.19 | 0.17 | 0.19 |
| full_xc015_coh040216_err39_fixedop | 0.36 | 0.19 | 0.18 | 0.17 | 0.02 | 0.04 | 0.10 | 0.14 | 0.21 | 0.35 | 0.18 | 0.15 | 0.17 |
| dense075_coh05_fixedop | 0.11 | 0.09 | 0.06 | 0.00 | 0.00 | 0.01 | 0.10 | 0.08 | 0.14 | 0.16 | 0.09 | 0.07 | 0.08 |
| smap_fixedop | 0.67 | 0.52 | 0.48 | 0.46 | 0.21 | 0.33 | 0.36 | 0.33 | 0.56 | 0.64 | 0.56 | 0.58 | 0.48 |
| l3_fixedop | 0.45 | 0.28 | 0.16 | 0.32 | 0.16 | 0.15 | 0.09 | 0.14 | 0.27 | 0.37 | 0.27 | 0.16 | 0.19 |

#### SMAP — alpha_opt

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 1.12 | 0.82 | 0.53 | 0.82 | 0.28 | 0.36 | 0.36 | 0.26 | 0.57 | 0.31 | 0.31 | 0.39 | 0.48 |
| full_xc015_err39_fixedop | 1.42 | 0.93 | 0.58 | 0.96 | 0.25 | 0.37 | 0.35 | 0.27 | 0.69 | 0.40 | 0.36 | 0.40 | 0.55 |
| full_xc015_coh040216_fixedop | 1.69 | 1.38 | 1.21 | 1.22 | 0.14 | 0.27 | 0.38 | 0.38 | 0.61 | 0.56 | 0.40 | 0.46 | 0.68 |
| full_xc015_coh040216_err39_fixedop | 2.26 | 1.80 | 1.57 | 1.58 | -0.02 | 0.19 | 0.34 | 0.34 | 0.69 | 0.71 | 0.46 | 0.45 | 0.80 |
| dense075_coh05_fixedop | 0.77 | 0.72 | 0.62 | -0.00 | -0.25 | -0.33 | 0.44 | 0.17 | 0.64 | 0.48 | 0.25 | 0.24 | 0.32 |
| smap_fixedop | 1.87 | 1.72 | 1.70 | 1.42 | 0.91 | 1.39 | 1.27 | 1.28 | 1.39 | 1.31 | 1.29 | 1.59 | 1.44 |
| l3_fixedop | 2.38 | 2.22 | 0.05 | 2.63 | 1.97 | 1.67 | 0.66 | -0.13 | 0.89 | 0.43 | 0.57 | 0.34 | 1.03 |

#### SMAP — noise

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 7.9 | 4.2 | 3.7 | 8.7 | 10.6 | 8.9 | 7.1 | 9.2 | 10.9 | 26.3 | 20.2 | 9.9 | 8.8 |
| full_xc015_err39_fixedop | 4.6 | 2.3 | 2.0 | 5.0 | 6.0 | 4.7 | 4.0 | 4.9 | 6.3 | 15.4 | 11.9 | 6.1 | 5.0 |
| full_xc015_coh040216_fixedop | 5.4 | 2.8 | 2.4 | 5.0 | 5.8 | 5.7 | 3.6 | 6.0 | 8.6 | 19.7 | 15.6 | 7.4 | 5.9 |
| full_xc015_coh040216_err39_fixedop | 2.9 | 1.4 | 1.2 | 2.5 | 2.9 | 2.9 | 1.9 | 3.2 | 5.0 | 12.0 | 9.2 | 4.5 | 3.2 |
| dense075_coh05_fixedop | 2.3 | 1.6 | 1.2 | 1.8 | 3.6 | 3.1 | 1.8 | 2.5 | 3.2 | 6.2 | 7.0 | 3.6 | 2.6 |
| smap_fixedop | 13.1 | 9.7 | 8.7 | 10.9 | 12.5 | 11.8 | 11.2 | 12.7 | 18.4 | 25.7 | 19.3 | 13.7 | 12.7 |
| l3_fixedop | 3.6 | 1.9 | 1.7 | 4.1 | 3.1 | 2.5 | 1.1 | 3.6 | 2.4 | 13.7 | 9.3 | 3.5 | 3.3 |

#### SMAP — gain

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -17.8 | -6.9 | -3.9 | -14.2 | -5.8 | -6.4 | -5.1 | -4.8 | -12.4 | -16.3 | -12.4 | -7.6 | -8.5 |
| full_xc015_err39_fixedop | -13.0 | -4.2 | -2.3 | -9.7 | -3.0 | -3.5 | -2.8 | -2.7 | -8.7 | -12.5 | -8.6 | -4.8 | -5.5 |
| full_xc015_coh040216_fixedop | -18.4 | -7.7 | -5.7 | -12.2 | -1.6 | -3.1 | -2.8 | -4.5 | -10.4 | -22.0 | -12.4 | -6.8 | -8.0 |
| full_xc015_coh040216_err39_fixedop | -12.9 | -5.0 | -3.7 | -7.8 | 0.1 | -1.1 | -1.3 | -2.1 | -6.9 | -17.1 | -8.5 | -4.1 | -5.2 |
| dense075_coh05_fixedop | -3.5 | -2.3 | -1.5 | 0.0 | 1.8 | 2.0 | -1.6 | -0.9 | -4.1 | -6.0 | -3.5 | -1.7 | -1.6 |
| smap_fixedop | -49.2 | -33.3 | -29.7 | -31.0 | -22.6 | -32.9 | -28.4 | -32.5 | -50.9 | -67.2 | -50.0 | -43.5 | -36.4 |
| l3_fixedop | -17.0 | -8.7 | -0.2 | -21.3 | -12.4 | -8.3 | -1.4 | 0.9 | -4.3 | -11.7 | -10.6 | -2.4 | -6.9 |

#### SMAP — dMSE

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -9.9 | -2.7 | -0.2 | -5.5 | +4.8 | +2.5 | +2.0 | +4.4 | -1.5 | +10.0 | +7.8 | +2.3 | +0.3 |
| full_xc015_err39_fixedop | -8.4 | -2.0 | -0.3 | -4.6 | +2.9 | +1.2 | +1.2 | +2.2 | -2.4 | +2.9 | +3.4 | +1.3 | -0.5 |
| full_xc015_coh040216_fixedop | -13.0 | -4.9 | -3.3 | -7.2 | +4.1 | +2.6 | +0.8 | +1.4 | -1.8 | -2.3 | +3.2 | +0.6 | -2.2 |
| full_xc015_coh040216_err39_fixedop | -10.1 | -3.6 | -2.5 | -5.3 | +3.0 | +1.8 | +0.6 | +1.0 | -1.9 | -5.0 | +0.7 | +0.5 | -1.9 |
| dense075_coh05_fixedop | -1.2 | -0.7 | -0.3 | +1.8 | +5.4 | +5.1 | +0.2 | +1.6 | -0.9 | +0.2 | +3.4 | +1.9 | +0.9 |
| smap_fixedop | -36.1 | -23.6 | -21.0 | -20.1 | -10.1 | -21.1 | -17.3 | -19.9 | -32.6 | -41.4 | -30.7 | -29.8 | -23.8 |
| l3_fixedop | -13.4 | -6.7 | +1.6 | -17.3 | -9.3 | -5.8 | -0.3 | +4.5 | -1.9 | +1.9 | -1.3 | +1.1 | -3.6 |

#### SMOS — corr

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 0.33 | 0.20 | 0.09 | 0.12 | 0.06 | 0.10 | 0.10 | 0.08 | 0.18 | 0.22 | 0.15 | 0.14 | 0.14 |
| full_xc015_err39_fixedop | 0.32 | 0.17 | 0.07 | 0.11 | 0.05 | 0.10 | 0.09 | 0.07 | 0.18 | 0.23 | 0.14 | 0.13 | 0.13 |
| full_xc015_coh040216_fixedop | 0.41 | 0.23 | 0.12 | 0.18 | 0.04 | 0.08 | 0.10 | 0.11 | 0.21 | 0.32 | 0.17 | 0.18 | 0.18 |
| full_xc015_coh040216_err39_fixedop | 0.39 | 0.22 | 0.10 | 0.16 | 0.03 | 0.07 | 0.08 | 0.10 | 0.20 | 0.33 | 0.16 | 0.16 | 0.17 |
| dense075_coh05_fixedop | 0.10 | 0.12 | 0.06 | 0.01 | 0.01 | 0.04 | 0.05 | 0.09 | 0.14 | 0.17 | 0.06 | 0.04 | 0.08 |
| smap_fixedop | 0.70 | 0.60 | 0.48 | 0.54 | 0.28 | 0.39 | 0.47 | 0.52 | 0.67 | 0.66 | 0.60 | 0.53 | 0.55 |
| l3_fixedop | 0.45 | 0.24 | 0.09 | 0.32 | 0.14 | 0.11 | 0.11 | 0.17 | 0.30 | 0.35 | 0.26 | 0.11 | 0.18 |

#### SMOS — alpha_opt

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 1.32 | 1.17 | 0.40 | 0.71 | 0.24 | 0.50 | 0.27 | 0.18 | 0.56 | 0.31 | 0.34 | 0.44 | 0.49 |
| full_xc015_err39_fixedop | 1.70 | 1.43 | 0.40 | 0.79 | 0.18 | 0.54 | 0.24 | 0.08 | 0.69 | 0.39 | 0.41 | 0.48 | 0.55 |
| full_xc015_coh040216_fixedop | 1.85 | 1.59 | 0.72 | 1.24 | 0.25 | 0.47 | 0.26 | 0.17 | 0.59 | 0.60 | 0.40 | 0.55 | 0.68 |
| full_xc015_coh040216_err39_fixedop | 2.42 | 2.04 | 0.85 | 1.52 | 0.07 | 0.47 | 0.15 | 0.01 | 0.69 | 0.76 | 0.49 | 0.59 | 0.78 |
| dense075_coh05_fixedop | 0.73 | 1.05 | 0.81 | 0.16 | -0.22 | -0.08 | 0.08 | 0.16 | 0.69 | 0.65 | 0.17 | 0.11 | 0.34 |
| smap_fixedop | 1.71 | 1.82 | 1.66 | 1.65 | 1.00 | 1.48 | 1.68 | 1.51 | 1.61 | 1.60 | 1.36 | 1.39 | 1.57 |
| l3_fixedop | 2.63 | 2.30 | -0.38 | 2.77 | 1.80 | 1.45 | 0.77 | -0.14 | 1.07 | 0.41 | 0.74 | 0.25 | 1.08 |

#### SMOS — noise

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 7.9 | 4.2 | 2.9 | 8.4 | 7.9 | 7.1 | 6.0 | 8.2 | 9.2 | 21.3 | 16.5 | 8.2 | 7.7 |
| full_xc015_err39_fixedop | 4.5 | 2.3 | 1.6 | 4.9 | 4.4 | 3.7 | 3.5 | 4.8 | 5.3 | 12.6 | 9.8 | 5.1 | 4.4 |
| full_xc015_coh040216_fixedop | 5.6 | 3.0 | 2.0 | 4.4 | 4.2 | 4.9 | 3.2 | 5.5 | 7.4 | 12.4 | 13.0 | 6.5 | 5.1 |
| full_xc015_coh040216_err39_fixedop | 3.0 | 1.6 | 1.0 | 2.3 | 2.2 | 2.5 | 1.8 | 3.0 | 4.2 | 7.6 | 7.6 | 4.0 | 2.9 |
| dense075_coh05_fixedop | 2.0 | 1.6 | 0.6 | 1.5 | 2.7 | 3.1 | 1.8 | 2.6 | 2.5 | 3.9 | 5.5 | 2.4 | 2.1 |
| smap_fixedop | 17.4 | 11.7 | 9.2 | 14.1 | 12.6 | 13.6 | 10.3 | 15.8 | 19.1 | 18.2 | 19.9 | 14.9 | 13.7 |
| l3_fixedop | 3.5 | 2.0 | 1.2 | 3.8 | 2.9 | 2.0 | 1.0 | 3.3 | 1.8 | 9.1 | 6.6 | 2.9 | 2.8 |

#### SMOS — gain

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -20.7 | -9.9 | -2.3 | -11.9 | -3.8 | -7.0 | -3.3 | -2.9 | -10.3 | -13.2 | -11.3 | -7.3 | -7.5 |
| full_xc015_err39_fixedop | -15.4 | -6.6 | -1.2 | -7.7 | -1.6 | -4.0 | -1.7 | -0.7 | -7.3 | -10.0 | -8.0 | -4.9 | -4.9 |
| full_xc015_coh040216_fixedop | -20.9 | -9.6 | -2.9 | -11.1 | -2.1 | -4.6 | -1.7 | -1.9 | -8.7 | -14.9 | -10.5 | -7.2 | -6.9 |
| full_xc015_coh040216_err39_fixedop | -14.5 | -6.5 | -1.8 | -7.0 | -0.3 | -2.4 | -0.5 | -0.1 | -5.8 | -11.6 | -7.5 | -4.7 | -4.5 |
| dense075_coh05_fixedop | -2.9 | -3.3 | -1.0 | -0.5 | 1.2 | 0.5 | -0.3 | -0.8 | -3.5 | -5.0 | -1.9 | -0.5 | -1.4 |
| smap_fixedop | -59.6 | -42.5 | -30.4 | -46.3 | -25.4 | -40.3 | -34.6 | -47.7 | -61.3 | -58.1 | -54.2 | -41.5 | -43.1 |
| l3_fixedop | -18.5 | -9.2 | 0.9 | -21.0 | -10.6 | -5.8 | -1.6 | 0.9 | -3.9 | -7.4 | -9.8 | -1.5 | -6.0 |

#### SMOS — dMSE

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -12.8 | -5.7 | +0.6 | -3.5 | +4.1 | +0.0 | +2.7 | +5.3 | -1.1 | +8.0 | +5.2 | +1.0 | +0.2 |
| full_xc015_err39_fixedop | -10.9 | -4.3 | +0.3 | -2.8 | +2.7 | -0.3 | +1.8 | +4.0 | -2.0 | +2.7 | +1.8 | +0.2 | -0.5 |
| full_xc015_coh040216_fixedop | -15.2 | -6.6 | -0.9 | -6.6 | +2.1 | +0.3 | +1.6 | +3.7 | -1.3 | -2.5 | +2.6 | -0.7 | -1.8 |
| full_xc015_coh040216_err39_fixedop | -11.5 | -4.9 | -0.7 | -4.7 | +1.9 | +0.1 | +1.2 | +2.9 | -1.6 | -4.0 | +0.1 | -0.7 | -1.6 |
| dense075_coh05_fixedop | -0.9 | -1.7 | -0.4 | +1.0 | +3.9 | +3.6 | +1.5 | +1.7 | -0.9 | -1.2 | +3.6 | +1.8 | +0.7 |
| smap_fixedop | -42.2 | -30.8 | -21.3 | -32.3 | -12.7 | -26.7 | -24.3 | -31.9 | -42.2 | -39.9 | -34.3 | -26.7 | -29.4 |
| l3_fixedop | -15.0 | -7.2 | +2.2 | -17.2 | -7.7 | -3.8 | -0.6 | +4.1 | -2.1 | +1.7 | -3.2 | +1.5 | -3.2 |

#### L3 — corr

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 0.20 | 0.11 | 0.13 | 0.13 | 0.12 | 0.17 | 0.19 | 0.17 | 0.19 | 0.19 | 0.24 | 0.17 | 0.22 |
| full_xc015_err39_fixedop | 0.18 | 0.09 | 0.11 | 0.12 | 0.11 | 0.17 | 0.19 | 0.16 | 0.18 | 0.18 | 0.24 | 0.16 | 0.22 |
| full_xc015_coh040216_fixedop | 0.25 | 0.17 | 0.16 | 0.12 | 0.09 | 0.15 | 0.18 | 0.20 | 0.22 | 0.16 | 0.30 | 0.15 | 0.25 |
| full_xc015_coh040216_err39_fixedop | 0.24 | 0.15 | 0.15 | 0.11 | 0.08 | 0.15 | 0.18 | 0.18 | 0.21 | 0.16 | 0.29 | 0.15 | 0.25 |
| dense075_coh05_fixedop | 0.08 | 0.11 | 0.10 | 0.05 | 0.03 | 0.07 | 0.09 | 0.08 | 0.12 | 0.07 | 0.15 | 0.10 | 0.11 |
| smap_fixedop | 0.39 | 0.26 | 0.20 | 0.16 | 0.00 | 0.03 | 0.20 | 0.11 | 0.29 | 0.22 | 0.33 | 0.21 | 0.16 |
| l3_fixedop | 0.38 | 0.29 | 0.18 | 0.26 | 0.19 | 0.33 | 0.30 | 0.27 | 0.31 | 0.24 | 0.37 | 0.29 | 0.46 |

#### L3 — alpha_opt

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 0.96 | 0.93 | 0.63 | 0.94 | 0.32 | 0.46 | 0.51 | 0.54 | 0.44 | 0.80 | 0.65 | 0.50 | 0.64 |
| full_xc015_err39_fixedop | 1.24 | 1.16 | 0.82 | 1.17 | 0.38 | 0.61 | 0.68 | 0.79 | 0.59 | 1.08 | 0.86 | 0.65 | 0.85 |
| full_xc015_coh040216_fixedop | 1.29 | 1.35 | 0.78 | 1.11 | 0.30 | 0.50 | 0.65 | 0.91 | 0.75 | 1.01 | 0.92 | 0.59 | 0.86 |
| full_xc015_coh040216_err39_fixedop | 1.73 | 1.76 | 1.03 | 1.40 | 0.33 | 0.67 | 0.86 | 1.29 | 1.03 | 1.38 | 1.22 | 0.79 | 1.15 |
| dense075_coh05_fixedop | 0.59 | 0.66 | 0.57 | 0.30 | -0.00 | 0.24 | 0.49 | 0.77 | 0.64 | 0.81 | 0.81 | 0.57 | 0.56 |
| smap_fixedop | 1.02 | 0.83 | 0.00 | 0.88 | 0.22 | 0.17 | 0.36 | -0.04 | 0.41 | 0.55 | 0.66 | 0.30 | 0.42 |
| l3_fixedop | 3.10 | 2.96 | 2.74 | 3.17 | 1.61 | 1.78 | 2.10 | 2.12 | 2.29 | 1.92 | 1.49 | 1.63 | 2.28 |

#### L3 — noise

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 7.9 | 5.5 | 6.8 | 7.6 | 14.9 | 14.9 | 14.6 | 9.4 | 16.5 | 15.4 | 19.6 | 16.1 | 11.2 |
| full_xc015_err39_fixedop | 4.5 | 3.0 | 3.7 | 4.1 | 8.0 | 7.8 | 8.0 | 5.0 | 9.6 | 8.9 | 11.3 | 9.8 | 6.2 |
| full_xc015_coh040216_fixedop | 5.6 | 4.0 | 4.8 | 4.5 | 8.5 | 9.7 | 8.1 | 6.3 | 12.9 | 10.7 | 15.0 | 12.0 | 7.6 |
| full_xc015_coh040216_err39_fixedop | 2.9 | 2.0 | 2.4 | 2.2 | 4.3 | 4.9 | 4.3 | 3.3 | 7.1 | 6.5 | 8.8 | 7.0 | 4.1 |
| dense075_coh05_fixedop | 2.2 | 2.0 | 2.0 | 1.7 | 4.9 | 5.2 | 3.7 | 2.3 | 4.6 | 3.6 | 6.4 | 4.8 | 3.1 |
| smap_fixedop | 17.3 | 15.0 | 20.3 | 10.4 | 12.5 | 21.9 | 22.9 | 15.1 | 23.6 | 10.2 | 23.6 | 21.1 | 16.9 |
| l3_fixedop | 3.4 | 2.9 | 3.0 | 3.8 | 3.9 | 3.7 | 2.1 | 3.2 | 3.4 | 7.4 | 9.6 | 5.2 | 4.1 |

#### L3 — gain

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -15.2 | -10.3 | -8.6 | -14.2 | -9.6 | -13.7 | -14.9 | -10.2 | -14.6 | -24.6 | -25.5 | -16.0 | -14.4 |
| full_xc015_err39_fixedop | -11.1 | -7.0 | -6.0 | -9.5 | -6.0 | -9.5 | -10.8 | -7.9 | -11.4 | -19.3 | -19.5 | -12.6 | -10.6 |
| full_xc015_coh040216_fixedop | -14.5 | -10.7 | -7.4 | -10.0 | -5.1 | -9.7 | -10.6 | -11.5 | -19.4 | -21.6 | -27.6 | -14.1 | -13.0 |
| full_xc015_coh040216_err39_fixedop | -10.1 | -7.0 | -4.9 | -6.2 | -2.9 | -6.5 | -7.4 | -8.4 | -14.7 | -17.8 | -21.4 | -11.2 | -9.5 |
| dense075_coh05_fixedop | -2.6 | -2.6 | -2.3 | -1.0 | 0.0 | -2.4 | -3.6 | -3.5 | -5.9 | -5.9 | -10.4 | -5.5 | -3.5 |
| smap_fixedop | -35.2 | -24.8 | -0.2 | -18.3 | -5.5 | -7.6 | -16.3 | 1.1 | -19.3 | -11.3 | -31.1 | -12.6 | -14.0 |
| l3_fixedop | -21.1 | -17.0 | -16.3 | -24.2 | -12.4 | -13.2 | -8.7 | -13.4 | -15.6 | -28.4 | -28.6 | -16.8 | -18.7 |

#### L3 — dMSE

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -7.3 | -4.8 | -1.8 | -6.6 | +5.2 | +1.2 | -0.3 | -0.8 | +1.9 | -9.2 | -5.8 | +0.1 | -3.2 |
| full_xc015_err39_fixedop | -6.6 | -4.0 | -2.3 | -5.5 | +2.0 | -1.7 | -2.8 | -2.9 | -1.8 | -10.4 | -8.2 | -2.9 | -4.3 |
| full_xc015_coh040216_fixedop | -8.9 | -6.8 | -2.6 | -5.5 | +3.4 | -0.0 | -2.5 | -5.2 | -6.5 | -10.9 | -12.6 | -2.1 | -5.5 |
| full_xc015_coh040216_err39_fixedop | -7.2 | -5.0 | -2.5 | -4.0 | +1.4 | -1.6 | -3.1 | -5.2 | -7.6 | -11.3 | -12.6 | -4.2 | -5.4 |
| dense075_coh05_fixedop | -0.4 | -0.6 | -0.3 | +0.7 | +4.9 | +2.7 | +0.1 | -1.2 | -1.3 | -2.2 | -4.0 | -0.7 | -0.3 |
| smap_fixedop | -17.9 | -9.8 | +20.1 | -7.9 | +7.0 | +14.3 | +6.6 | +16.2 | +4.3 | -1.1 | -7.5 | +8.5 | +2.8 |
| l3_fixedop | -17.7 | -14.1 | -13.3 | -20.4 | -8.6 | -9.5 | -6.7 | -10.2 | -12.2 | -21.0 | -19.0 | -11.6 | -14.6 |

#### L1 — corr

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 0.12 | 0.07 | 0.11 | 0.06 | 0.07 | 0.07 | 0.06 | 0.07 | 0.12 | 0.05 | 0.10 | 0.10 | 0.09 |
| full_xc015_err39_fixedop | 0.11 | 0.06 | 0.10 | 0.05 | 0.06 | 0.06 | 0.05 | 0.07 | 0.12 | 0.05 | 0.09 | 0.10 | 0.08 |
| full_xc015_coh040216_fixedop | 0.14 | 0.10 | 0.17 | 0.09 | 0.05 | 0.07 | 0.09 | 0.12 | 0.13 | 0.13 | 0.14 | 0.12 | 0.13 |
| full_xc015_coh040216_err39_fixedop | 0.12 | 0.10 | 0.15 | 0.08 | 0.05 | 0.06 | 0.08 | 0.11 | 0.13 | 0.13 | 0.14 | 0.11 | 0.12 |
| dense075_coh05_fixedop | 0.03 | 0.07 | 0.09 | 0.04 | 0.05 | 0.04 | 0.06 | 0.07 | 0.08 | 0.08 | 0.10 | 0.12 | 0.07 |
| smap_fixedop | 0.18 | 0.11 | 0.22 | 0.10 | 0.02 | 0.05 | 0.08 | 0.10 | 0.16 | 0.09 | 0.10 | 0.11 | 0.11 |
| l3_fixedop | 0.13 | 0.14 | 0.15 | 0.11 | 0.10 | 0.10 | 0.08 | 0.11 | 0.11 | 0.06 | 0.14 | 0.12 | 0.11 |

#### L1 — alpha_opt

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 1.08 | 1.13 | 0.95 | 1.10 | 0.67 | 0.66 | 0.69 | 0.61 | 0.85 | 0.55 | 0.60 | 0.55 | 0.73 |
| full_xc015_err39_fixedop | 1.30 | 1.43 | 1.10 | 1.28 | 0.76 | 0.74 | 0.78 | 0.74 | 1.00 | 0.71 | 0.72 | 0.65 | 0.86 |
| full_xc015_coh040216_fixedop | 1.21 | 1.43 | 1.48 | 1.66 | 0.57 | 0.59 | 0.63 | 0.76 | 0.87 | 0.83 | 0.81 | 0.57 | 0.88 |
| full_xc015_coh040216_err39_fixedop | 1.49 | 1.87 | 1.87 | 2.10 | 0.62 | 0.66 | 0.77 | 0.93 | 1.09 | 1.05 | 1.03 | 0.61 | 1.07 |
| dense075_coh05_fixedop | 0.54 | 1.25 | 1.34 | 0.62 | 0.69 | 0.42 | 0.75 | 0.91 | 1.05 | 0.90 | 1.00 | 1.15 | 0.88 |
| smap_fixedop | 1.37 | 1.15 | 1.46 | 1.24 | 0.31 | 0.99 | 0.69 | 1.08 | 1.22 | 0.85 | 0.69 | 0.71 | 1.01 |
| l3_fixedop | 2.42 | 4.09 | 1.32 | 3.93 | 2.57 | 2.29 | 2.21 | 1.09 | 0.70 | 1.03 | 0.99 | 0.80 | 1.51 |

#### L1 — noise

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 1.3 | 0.9 | 1.3 | 1.2 | 1.1 | 1.3 | 0.8 | 1.2 | 1.9 | 1.8 | 2.4 | 2.9 | 1.4 |
| full_xc015_err39_fixedop | 0.7 | 0.5 | 0.7 | 0.7 | 0.6 | 0.7 | 0.5 | 0.7 | 1.2 | 1.1 | 1.4 | 1.8 | 0.8 |
| full_xc015_coh040216_fixedop | 1.3 | 1.0 | 1.2 | 1.0 | 1.0 | 1.4 | 1.8 | 2.1 | 2.3 | 4.3 | 3.8 | 3.9 | 2.0 |
| full_xc015_coh040216_err39_fixedop | 0.7 | 0.5 | 0.6 | 0.5 | 0.5 | 0.7 | 1.0 | 1.1 | 1.3 | 2.6 | 2.3 | 2.3 | 1.1 |
| dense075_coh05_fixedop | 0.4 | 0.3 | 0.4 | 0.4 | 0.5 | 0.5 | 0.8 | 0.5 | 0.5 | 1.0 | 1.0 | 1.1 | 0.6 |
| smap_fixedop | 1.7 | 1.3 | 2.4 | 1.2 | 0.7 | 1.1 | 0.8 | 1.1 | 1.8 | 1.2 | 1.9 | 2.5 | 1.3 |
| l3_fixedop | 0.3 | 0.2 | 0.3 | 0.4 | 0.2 | 0.3 | 0.1 | 0.4 | 0.4 | 1.0 | 1.1 | 0.8 | 0.5 |

#### L1 — gain

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -2.8 | -2.0 | -2.4 | -2.6 | -1.5 | -1.7 | -1.2 | -1.5 | -3.3 | -2.0 | -2.9 | -3.2 | -2.0 |
| full_xc015_err39_fixedop | -1.9 | -1.4 | -1.5 | -1.7 | -0.9 | -1.0 | -0.7 | -1.0 | -2.3 | -1.5 | -2.1 | -2.3 | -1.4 |
| full_xc015_coh040216_fixedop | -3.2 | -2.9 | -3.5 | -3.3 | -1.2 | -1.7 | -2.3 | -3.2 | -4.0 | -7.0 | -6.1 | -4.4 | -3.4 |
| full_xc015_coh040216_err39_fixedop | -2.1 | -1.9 | -2.2 | -2.1 | -0.6 | -1.0 | -1.5 | -2.1 | -2.9 | -5.4 | -4.7 | -2.8 | -2.4 |
| dense075_coh05_fixedop | -0.4 | -0.8 | -1.2 | -0.5 | -0.6 | -0.4 | -1.2 | -0.9 | -1.2 | -1.8 | -2.0 | -2.5 | -1.0 |
| smap_fixedop | -4.8 | -3.0 | -6.9 | -3.0 | -0.4 | -2.2 | -1.2 | -2.4 | -4.4 | -2.1 | -2.7 | -3.5 | -2.7 |
| l3_fixedop | -1.3 | -2.0 | -0.9 | -3.3 | -1.1 | -1.2 | -0.7 | -0.9 | -0.6 | -2.1 | -2.2 | -1.3 | -1.4 |

#### L1 — dMSE

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -1.5 | -1.1 | -1.1 | -1.4 | -0.4 | -0.4 | -0.3 | -0.3 | -1.4 | -0.2 | -0.5 | -0.3 | -0.6 |
| full_xc015_err39_fixedop | -1.2 | -0.9 | -0.8 | -1.0 | -0.3 | -0.3 | -0.3 | -0.3 | -1.2 | -0.5 | -0.6 | -0.5 | -0.6 |
| full_xc015_coh040216_fixedop | -1.9 | -1.9 | -2.3 | -2.3 | -0.1 | -0.3 | -0.5 | -1.1 | -1.7 | -2.8 | -2.3 | -0.5 | -1.5 |
| full_xc015_coh040216_err39_fixedop | -1.4 | -1.4 | -1.6 | -1.6 | -0.1 | -0.2 | -0.5 | -1.0 | -1.6 | -2.8 | -2.4 | -0.5 | -1.3 |
| dense075_coh05_fixedop | -0.0 | -0.5 | -0.7 | -0.1 | -0.2 | +0.1 | -0.4 | -0.4 | -0.6 | -0.8 | -1.0 | -1.4 | -0.4 |
| smap_fixedop | -3.0 | -1.7 | -4.6 | -1.8 | +0.3 | -1.1 | -0.3 | -1.3 | -2.6 | -0.9 | -0.7 | -1.0 | -1.4 |
| l3_fixedop | -1.0 | -1.8 | -0.6 | -2.9 | -0.9 | -0.9 | -0.5 | -0.5 | -0.2 | -1.1 | -1.1 | -0.5 | -0.9 |

#### ASCAT — corr

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 0.21 | 0.12 | 0.06 | -0.02 | 0.04 | 0.03 | 0.10 | 0.10 | 0.07 | 0.02 | 0.11 | 0.11 | 0.09 |
| full_xc015_err39_fixedop | 0.20 | 0.12 | 0.06 | -0.03 | 0.03 | 0.03 | 0.12 | 0.10 | 0.07 | 0.03 | 0.11 | 0.12 | 0.10 |
| full_xc015_coh040216_fixedop | 0.22 | 0.13 | 0.05 | -0.07 | 0.03 | 0.14 | 0.10 | 0.08 | 0.05 | 0.06 | 0.11 | 0.10 | 0.10 |
| full_xc015_coh040216_err39_fixedop | 0.22 | 0.13 | 0.05 | -0.07 | 0.02 | 0.13 | 0.10 | 0.08 | 0.06 | 0.07 | 0.12 | 0.10 | 0.10 |
| dense075_coh05_fixedop | 0.05 | 0.06 | 0.04 | 0.01 | 0.00 | 0.10 | 0.05 | 0.03 | 0.04 | 0.05 | 0.07 | 0.10 | 0.05 |
| smap_fixedop | 0.39 | 0.30 | 0.20 | 0.11 | 0.07 | 0.00 | 0.11 | 0.12 | 0.21 | 0.12 | 0.16 | -0.04 | 0.09 |
| l3_fixedop | 0.30 | 0.24 | 0.09 | 0.07 | 0.17 | 0.25 | 0.29 | 0.20 | 0.18 | 0.17 | 0.22 | 0.19 | 0.33 |

#### ASCAT — alpha_opt

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 0.92 | 0.61 | 0.14 | 0.01 | 0.17 | 0.15 | 0.30 | 0.24 | 0.11 | 0.21 | 0.18 | 0.46 | 0.25 |
| full_xc015_err39_fixedop | 1.23 | 0.86 | 0.20 | -0.04 | 0.23 | 0.22 | 0.48 | 0.42 | 0.19 | 0.44 | 0.29 | 0.74 | 0.40 |
| full_xc015_coh040216_fixedop | 1.05 | 0.68 | 0.08 | -0.13 | 0.17 | 0.72 | 0.40 | 0.31 | 0.15 | 0.37 | 0.32 | 0.60 | 0.31 |
| full_xc015_coh040216_err39_fixedop | 1.45 | 0.95 | 0.09 | -0.26 | 0.22 | 0.98 | 0.59 | 0.52 | 0.25 | 0.63 | 0.48 | 0.94 | 0.49 |
| dense075_coh05_fixedop | 0.45 | 0.44 | 0.06 | -0.02 | 0.19 | 0.93 | 0.38 | 0.26 | 0.18 | 0.35 | 0.43 | 0.69 | 0.29 |
| smap_fixedop | 1.06 | 0.61 | 0.16 | 0.32 | -0.04 | -0.03 | 0.01 | -0.15 | 0.16 | -0.44 | 0.10 | -0.24 | 0.11 |
| l3_fixedop | 2.55 | 1.65 | 1.55 | 0.89 | 0.69 | 1.22 | 1.81 | 1.41 | 1.10 | 1.76 | 1.18 | 1.71 | 1.52 |

#### ASCAT — noise

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | 8.5 | 6.5 | 11.4 | 10.7 | 12.9 | 5.7 | 12.4 | 12.7 | 19.3 | 15.7 | 24.4 | 12.3 | 13.1 |
| full_xc015_err39_fixedop | 4.5 | 3.2 | 5.9 | 5.3 | 6.5 | 3.2 | 6.9 | 6.4 | 11.0 | 8.5 | 12.7 | 6.8 | 7.0 |
| full_xc015_coh040216_fixedop | 6.3 | 5.1 | 8.6 | 7.3 | 6.6 | 3.8 | 7.9 | 9.5 | 16.2 | 12.3 | 19.0 | 9.7 | 9.9 |
| full_xc015_coh040216_err39_fixedop | 3.2 | 2.4 | 4.3 | 3.3 | 3.0 | 1.8 | 4.4 | 4.6 | 8.8 | 7.3 | 10.1 | 5.6 | 5.2 |
| dense075_coh05_fixedop | 2.5 | 2.4 | 2.6 | 1.6 | 2.6 | 1.3 | 3.1 | 3.1 | 5.6 | 4.2 | 5.5 | 3.3 | 3.3 |
| smap_fixedop | 17.3 | 26.6 | 27.9 | 24.3 | 19.8 | 23.3 | 31.8 | 26.8 | 42.8 | 21.5 | 57.0 | 21.9 | 28.2 |
| l3_fixedop | 4.2 | 4.5 | 2.9 | 5.2 | 4.9 | 4.0 | 2.8 | 4.9 | 4.4 | 7.8 | 7.9 | 4.5 | 4.8 |

#### ASCAT — gain

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -15.6 | -8.0 | -3.2 | -0.3 | -4.5 | -1.7 | -7.5 | -6.2 | -4.4 | -6.7 | -9.0 | -11.3 | -6.6 |
| full_xc015_err39_fixedop | -11.1 | -5.6 | -2.4 | 0.4 | -3.0 | -1.4 | -6.7 | -5.5 | -4.3 | -7.4 | -7.4 | -10.1 | -5.6 |
| full_xc015_coh040216_fixedop | -13.4 | -7.0 | -1.4 | 1.9 | -2.2 | -5.5 | -6.3 | -5.9 | -4.8 | -9.0 | -12.0 | -11.6 | -6.2 |
| full_xc015_coh040216_err39_fixedop | -9.4 | -4.5 | -0.7 | 1.7 | -1.3 | -3.6 | -5.2 | -4.8 | -4.4 | -9.2 | -9.7 | -10.6 | -5.1 |
| dense075_coh05_fixedop | -2.3 | -2.2 | -0.3 | 0.0 | -1.0 | -2.4 | -2.4 | -1.6 | -2.0 | -2.9 | -4.8 | -4.6 | -1.9 |
| smap_fixedop | -36.4 | -32.5 | -8.9 | -15.4 | 1.6 | 1.3 | -0.5 | 7.8 | -14.0 | 18.9 | -11.5 | 10.4 | -6.2 |
| l3_fixedop | -21.2 | -14.9 | -9.0 | -9.2 | -6.7 | -9.7 | -10.0 | -13.9 | -9.6 | -27.4 | -18.7 | -15.3 | -14.6 |

#### ASCAT — dMSE

| arm | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 | 2020 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full_xc015_fixedop | -7.1 | -1.5 | +8.2 | +10.5 | +8.4 | +4.0 | +4.9 | +6.5 | +14.9 | +9.0 | +15.4 | +1.0 | +6.5 |
| full_xc015_err39_fixedop | -6.6 | -2.3 | +3.6 | +5.7 | +3.5 | +1.8 | +0.2 | +1.0 | +6.7 | +1.1 | +5.2 | -3.3 | +1.4 |
| full_xc015_coh040216_fixedop | -7.0 | -1.8 | +7.2 | +9.2 | +4.4 | -1.7 | +1.6 | +3.6 | +11.4 | +3.3 | +7.0 | -2.0 | +3.7 |
| full_xc015_coh040216_err39_fixedop | -6.1 | -2.1 | +3.6 | +5.0 | +1.7 | -1.8 | -0.8 | -0.2 | +4.4 | -1.9 | +0.4 | -4.9 | +0.1 |
| dense075_coh05_fixedop | +0.2 | +0.3 | +2.3 | +1.7 | +1.6 | -1.1 | +0.7 | +1.5 | +3.6 | +1.3 | +0.8 | -1.3 | +1.4 |
| smap_fixedop | -19.2 | -5.9 | +19.0 | +8.9 | +21.3 | +24.5 | +31.2 | +34.6 | +28.8 | +40.4 | +45.5 | +32.3 | +22.1 |
| l3_fixedop | -17.1 | -10.4 | -6.1 | -4.0 | -1.8 | -5.7 | -7.2 | -9.0 | -5.2 | -19.6 | -10.7 | -10.8 | -9.8 |

