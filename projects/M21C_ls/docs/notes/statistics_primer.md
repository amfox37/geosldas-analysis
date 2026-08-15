# Statistics primer (plain-language)

Plain-language description of the statistical methods used in the trend and
transition analyses. Intended as orientation material that can be adapted into
a short introductory paragraph in the Methods, or into the Discussion where the
conservative posture is defended. It is not a substitute for the formal Methods
text in `../manuscript/manuscript.md` (§2.7); the numbers, citations, and exact
procedures live there.

---

**Two questions, two tools.** The analysis keeps two different questions apart.
A *long-term trend* asks whether a quantity drifted steadily over 24 years. A
*transition* asks whether it stepped or changed pace at a specific moment — when
a new satellite entered. These are different estimands: a record can have a
sharp transition and no net trend, or a trend and no transition, so we never use
one method's answer for the other's question.

**Robust trends (Theil–Sen).** Rather than least-squares regression, which a few
extreme months can tilt, we use the Theil–Sen slope: the median of the slopes
between all pairs of points. It is insensitive to outliers and makes no
assumption that the noise is Gaussian — appropriate for monthly land fields.

**Honest significance under memory (Mann–Kendall + Hamed–Rao).** Monthly soil
moisture and snow are strongly autocorrelated — each month resembles the last —
so a naive test sees far more "significant" trends than really exist, because it
counts correlated months as independent evidence. We use the Mann–Kendall rank
test with a conservative Hamed–Rao correction that inflates the variance in
proportion to the series' own persistence, and by construction this correction
can only make a result *less* significant, never more.

**Controlling false discoveries (Benjamini–Hochberg).** Testing ~112,000 tiles,
5% would look "significant" by chance alone — thousands of false positives.
Benjamini–Hochberg false-discovery-rate control caps the expected fraction of
false positives among the tiles we flag, so a map of significant trends means
what it appears to mean.

**Transitions as interrupted time series.** At each known sensor-onset date we
fit one regression to the whole record containing a baseline trend, a seasonal
cycle, a *level jump* at the date, and — where a period is long enough — a
*change of slope*. This separates an abrupt offset from a change in rate.
Because the residuals are autocorrelated, we whiten them with an iterative
Prais–Winsten AR(1) fit and get uncertainty from a bootstrap that resamples the
model's own innovations under a deliberately conservative persistence — wide,
cautious intervals.

**An independent second opinion (PELT changepoints).** As a check that we aren't
just finding what we went looking for, a changepoint detector scans the same
series *without being told the sensor dates* and reports where the data
themselves break. Where an independent break lands on a known sensor date, the
evidence for a real regime change is much stronger. This detector sees abrupt
steps well but gradual hinges poorly, so it corroborates the level changes and
the known-date model remains authoritative for slopes.

**Spatial-block bootstrap for the water budget.** For the snow-water partition,
neighboring tiles are not independent — nearby locations share weather.
Resampling individual tiles would understate uncertainty. We instead resample
~5° spatial blocks, preserving local spatial dependence, so the confidence
intervals on the runoff/ET/storage fractions are honest.

**A recurring theme: conservatism.** Every choice — the variance inflation
floored at one, FDR over pointwise stippling, the conservative bootstrap
persistence, block resampling — is deliberately tilted toward under-claiming.
Consequently, where the analysis finds nothing, that means "not resolved," not
"no effect."
