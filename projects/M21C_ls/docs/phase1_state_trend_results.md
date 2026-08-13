# M21C Phase 1 OL and DA state-trend results

## Production contract

Twelve matched tile-level control fields were generated from the audited OL v2
and DA v3 monthly products for June 2000-May 2024. OL and DA use identical
paired-finite samples and the same masks as the corresponding primary
`DA - OL` result. All 12 fields pass the production audit. Spatial significance
uses modified Mann-Kendall followed by BH-FDR at 0.05 separately within each
complete tile field.

Area-weighted domain means use the same monthly loader and trend engine. Their
18 tests (six variables times OL, DA, and `DA - OL`) form one additional
BH-FDR family. The domain summaries describe a net spatial mean; tile summaries
describe geographically varying trends that may cancel in that mean.

## Background trends and DA effects

| Variable | OL significant tiles | DA significant tiles | DA - OL significant tiles | OL/DA slope correlation |
|---|---:|---:|---:|---:|
| Precipitation | 3,719 (3.30%) | 3,726 (3.31%) | 0 | 0.9998 |
| Surface soil moisture | 6,992 (6.21%) | 8,966 (7.96%) | 1,412 (1.25%) | 0.9304 |
| Root-zone soil moisture | 4,909 (4.36%) | 10,329 (9.18%) | 7,892 (7.01%) | 0.9137 |
| Snow mass | 12 (0.025%) | 5 (0.010%) | 0 | 0.9861 |
| Snow depth | 7 (0.015%) | 6 (0.012%) | 0 | 0.9633 |
| Snow-covered fraction | 1 (0.002%) | 0 | 0 | 0.6777 |

Precipitation contains geographically mixed background trends, split almost
equally between positive and negative tiles. OL and DA are nearly identical:
3,603 significant tiles overlap and all have the same direction. The
area-weighted domain trend is not significant, and no `DA - OL` tile survives
FDR. DA therefore does not measurably alter the precipitation trend field.

Surface soil moisture has a modest background trend pattern, predominantly
positive where significant. DA preserves its broad structure but increases the
number of significant tiles. The assimilation-specific `DA - OL` trend remains
limited to 1.25% of tiles and has no significant area-weighted domain trend.

Root-zone soil moisture is the qualified exception to a complete non-result.
DA increases significant spatial coverage from 4.36% in OL to 9.18%, while the
`DA - OL` trend is significant over 7.01% and predominantly positive. The
area-weighted OL, DA, and `DA - OL` trends are nevertheless all nonsignificant.
Thus DA modifies regional root-zone trends without producing a resolved global
land-mean secular trend.

## Snow

Snow mass and snow depth have effectively no tile-level trend in either OL or
DA, and no `DA - OL` tile survives FDR. Their area-weighted domain trends are
also nonsignificant.

Area-weighted snow-covered fraction declines in both OL and DA:

| Series | Sen slope (yr-1) | Autocorrelation-adjusted 95% interval | BH q |
|---|---:|---:|---:|
| OL | -0.000554 | -0.000911 to -0.000198 | 0.0061 |
| DA | -0.000549 | -0.000856 to -0.000233 | 0.0016 |
| DA - OL | -0.000006 | -0.000073 to +0.000065 | 0.9934 |

The OL and DA slopes correspond to a reduction of approximately 0.013 in
domain-mean snow-covered fraction over the 24-year record. The nearly identical
slopes and null `DA - OL` trend indicate that DA retains rather than creates
this background decline. The apparent contrast between a significant domain
mean and almost no individually FDR-significant tiles is possible because many
small, spatially coherent tile effects can be detectable after aggregation but
too weak to survive a 48,067-tile multiple-testing correction individually.

## Main reading

1. OL does contain trend structure, but not a significant net domain trend in
   precipitation, surface soil moisture, root-zone soil moisture, snow mass,
   or snow depth.
2. DA does not introduce a broad snow-state trend and does not alter the
   background precipitation trend.
3. DA has a limited surface-soil-moisture trend effect and a clearer regional
   root-zone effect, without a significant global land-mean trend.
4. Snow-covered fraction declines coherently in both OL and DA, while their
   difference has no trend. This is the cleanest example of DA preserving an OL
   background trend.

Derived NetCDFs and CSV summaries are under `output/trends_breakpoints/` and
are intentionally not versioned.
