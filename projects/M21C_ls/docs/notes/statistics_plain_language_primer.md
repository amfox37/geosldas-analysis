# Plain-Language Statistics Primer for the M21C Land-Reanalysis Paper

This document is an internal learning and reference guide. Its purpose is to explain what each statistical method in Sections 2.6 and 2.7 is doing, why it is needed, what assumptions it makes, and what conclusions it does and does not support. It is intentionally more repetitive and pedagogical than the manuscript or Supporting Information.

It is not manuscript prose and should never be pasted into the paper. The formal, publication-facing versions live in `../manuscript/manuscript.md` (§2.6, §2.7) and `../manuscript/supporting_information.md`. This file supersedes the earlier short orientation note `statistics_primer.md`, which was retired once §2.7 had been drafted and this document written.

---

## 1. The five ideas to understand first

Almost every statistical argument in this paper is an application of five ideas. If these are solid, the rest is detail.

### 1.1 Estimand versus estimator

**What question are we asking?** What exactly are we trying to learn, and what calculation are we using to learn it?

The **estimand** is the scientific quantity you actually want to know. The **estimator** is the calculation you perform to approximate it. They are not the same thing, and confusing them is the most common way that a technically correct calculation ends up answering the wrong question.

Two examples from this paper:

| Estimand (what we want) | Estimator (what we compute) |
|---|---|
| The fraction of snow-DA-added water that leaves as runoff | The area-weighted ratio of aggregate runoff response to aggregate snow input |
| The long-term DA-induced change in root-zone soil moisture | The Theil–Sen slope of the paired DA−OL monthly series |

**Why would a simpler approach be risky?** Because there is always more than one estimator for a given estimand, and they can disagree. For the snow partition, we could have averaged the individual tile-year ratios instead of taking the ratio of the sums. That is also a perfectly valid calculation — it just answers a different question, and in this case a badly behaved one (see §2.4).

**How would I explain this to a coauthor or reviewer?** "Before we discuss whether the number is right, let's agree on what number we were trying to produce." Most methodological arguments in reanalysis evaluation are estimand disagreements dressed up as estimator disagreements.

The practical lesson: choosing the right estimand matters more than choosing a sophisticated method. A Theil–Sen slope of the wrong series is worse than an ordinary least-squares slope of the right one.

### 1.2 Effect size versus uncertainty

**What question are we asking?** How big is the effect, and how well do we know it?

These are separate questions and they come apart in both directions. An estimate can be large but so uncertain that we cannot rule out zero. An estimate can be small but so precisely determined that zero is comfortably excluded.

The paper contains a clean pair of contrasting examples at the April 2015 (P6) boundary:

- **P6 total storage:** +2.144 kg m⁻², interval 0.932 to 3.266, q = 0.009. Large *and* well-resolved. This survives family correction.
- **P6 total runoff:** +0.230 kg m⁻² month⁻¹, interval −0.170 to +0.623, q = 0.478. The point estimate is positive, but the interval comfortably contains zero. We cannot distinguish this from no change.

**What does it NOT mean?** The runoff result does not mean runoff did not change at P6. It means our analysis cannot resolve whether it did. Those are very different claims, and the paper must say the second one.

### 1.3 Dependence

**What question are we asking?** How much genuinely independent information do we actually have?

There are two kinds of dependence in this analysis, and both inflate the *apparent* sample size well beyond the real one.

**Spatial dependence.** Neighbouring M36 tiles share weather. A Colorado mountain tile and the tile next to it experience much the same storms, much the same snowpack, and much the same melt timing. We have 48,067 tiles in the snow domain, but we emphatically do not have 48,067 independent experiments.

**Temporal dependence.** This month's soil moisture strongly resembles last month's. Land states have memory: they integrate past forcing rather than resetting each month. We have 288 monthly values, but they carry far less independent information than 288 independent draws would.

**Why would a simpler approach be risky?** Because standard formulas for uncertainty almost all assume independence, and the sample size enters those formulas directly. If you feed in 288 when the effective independent count is closer to 40, your standard errors come out roughly $\sqrt{288/40} \approx 2.7$ times too small, and things look significant that are not.

**How would I explain this to a coauthor?** "Asking the same person the same question 288 times does not give you 288 opinions." Autocorrelated data is closer to that than to 288 separate surveys.

Every uncertainty method in the current paper exists to handle one of these two dependences: the spatial-block bootstrap handles the first, while the Hamed–Rao correction and regional AR(1) effective sample size handle the second.

### 1.4 Pointwise significance versus multiple testing

**What question are we asking?** Is *this one* result surprising, or is it merely the most surprising of very many results?

Suppose you test one tile for a trend and get p = 0.03. That is reasonably surprising under the null hypothesis. Now suppose you test 112,573 tiles and flag every one with p < 0.05. Even if not a single tile has a real trend, you should expect roughly 5,600 tiles to clear that bar by chance alone. A map of them would look strikingly organised, because spatial dependence makes the false positives cluster.

So there are two different notions of "significant":

- a **p-value** answers "how surprising is this individual result?";
- a **q-value** under Benjamini–Hochberg answers "if I flag everything at least this extreme, what fraction of my flagged set is expected to be false?"

We use q-values for anything mapped. §3.6 goes into what q actually means, because it is very commonly misdescribed.

### 1.5 Period differences versus boundary jumps

**What question are we asking?** Did average assimilation influence differ between two predefined observing-system periods?

Figure 17 compares the mean of the adjusted regional series in adjacent periods, such as P6 minus P5. The observing-system dates define which months belong to each group, but the estimator is a difference between two averages.

That is different from estimating an instantaneous jump at the boundary. If the series changes gradually through P6, its P6 mean can differ from P5 even though nothing discontinuous happened on 1 April 2015.

**What does it NOT mean?** A resolved P6−P5 difference is not causal proof about SMAP. The correct phrasing is that the difference coincides with SMAP entry and the broader multi-sensor reorganization, not that SMAP uniquely caused it.

---

## 2. Section 2.6 — following snow-DA water through the system

### 2.1 Why the differential budget is unusually powerful

**What question are we asking?** When snow assimilation adds water to the snowpack, where does that water end up?

**Why is this a better question than it first appears?** Because snow data assimilation does something unusual: it literally adds and removes water mass from the model. This is not a subtle statistical influence on a parameter. It is a mass injection, and mass is conserved, so the added water has to go somewhere and we can go looking for it.

The OL and DA experiments share the same precipitation forcing. So when we subtract OL from DA, the common meteorological water input cancels to first order. What survives the subtraction is the assimilation effect and its consequences. The scientific question then becomes concrete and almost accounting-like:

> If snow DA added $X$ kg m⁻² of water to this tile this water year, how much larger was DA's runoff than OL's by the end of that water year? How much larger was its ET? How much more water was still sitting in storage?

**Why would a simpler approach be risky?** The obvious alternative is a correlation study: do tiles with more snow assimilation have more runoff? That is much weaker. Correlations do not have to sum to anything, they carry no units that let you check closure, and they cannot distinguish "assimilation added water which became runoff" from "snowy places have more runoff and also more assimilation." The budget can be checked for closure, which is a genuine constraint that a correlation study simply does not have.

**What assumptions remain?** Chiefly that OL and DA really do share forcing (we check this numerically rather than assuming it — the six-year integrated domain-mean precipitation difference is −5.5 × 10⁻⁵ kg m⁻², i.e. nothing), and that the storage endpoints from monthly means adequately approximate true water-year endpoints (they do not exactly, and that approximation is part of why there is a residual).

### 2.2 Why snowmelt is not a partition

**What question are we asking?** Which terms are *sinks* — places water ends up — and which are merely *waypoints*?

Suppose snow assimilation adds 100 units of water to the snowpack in January. Over the following spring, essentially all 100 units melt. Later, 64 of those units run off and 36 evaporate.

Now: how much water did assimilation add to the system? 100 units. Where did it go? 64 to runoff, 36 to ET.

What would happen if we included snowmelt as a partition term? We would report 100 units to melt, 64 to runoff, 36 to ET — a total of 200 units out of a 100-unit input. We would have counted the same water twice, once as it left the snowpack and again as it left the land surface.

The same logic applies to infiltration. Water that infiltrates has moved from the surface into the soil column. It has not left the land surface. It will subsequently either evaporate, drain as baseflow, or still be sitting in storage at the end of the year — and *those* are the terms we count.

So the terminal accounting terms are:

- **runoff** (surface runoff + baseflow) — water leaves laterally;
- **ET** — water leaves vertically to the atmosphere;
- **end-of-year storage change** — water has not left, it is still here;
- **residual** — the part that does not close, which we report honestly rather than hide.

And the pathway diagnostic terms are **snowmelt** and **infiltration**. They tell us *how* and *when* the water moved, which is scientifically interesting and appears in Figure 15. They are not sinks and never enter the closing equation.

**How would I explain this to a coauthor?** "Melt is a door, not a room. We're counting rooms."

### 2.3 Why RZMC and SFMC are not added to storage

**What question are we asking?** Is soil moisture a separate place water can be, or is it already counted?

`TWLAND` is total land-water storage. It already contains soil water. RZMC (root-zone soil moisture) and SFMC (surface soil moisture) are diagnostic states describing part of that same water. Adding a RZMC term to a budget that already has a `TWLAND` term would count the soil water twice.

This is worth stating explicitly because RZMC is scientifically the most interesting response variable in the whole snow analysis — it is where you actually see the snow-DA signal arrive in the soil column, with a clear seasonal timing and a clear peak in May. It is central to Figure 15. It is just not a mass term.

**What does a positive RZMC response mean?** That assimilation-added snow water reached the root zone and stayed long enough to be visible in monthly means. **What does it NOT mean?** It is not an additional quantity of water beyond what `TWLAND` already accounts for.

There is also a practical constraint: the compressed monthly archive contains SFMC, RZMC, and total `TWLAND`, but no separate native prognostic soil-water storage state. Even if we wanted a distinct soil-water mass term, the input product does not contain one.

### 2.4 What the direct partition actually calculates

**What question are we asking?** Of all the water that snow assimilation added across the domain over six years, what fraction ended up in each sink?

**What do we actually calculate?** The ratio of area-weighted sums:

$$
f_k = \frac{\sum A_i Y_{k,i,t}}{\sum A_i I_{\mathrm{snow},i,t}}
$$

Total response across the sample, divided by total input across the sample.

**Why would a simpler approach be risky?** The tempting alternative is to compute a ratio for each tile-year and average those ratios. Here is a two-tile toy example showing why that fails.

| Tile-year | Input | Runoff response | Individual ratio |
|---|---:|---:|---:|
| A | 100 kg m⁻² | 64 kg m⁻² | 0.64 |
| B | 0.1 kg m⁻² | 0.5 kg m⁻² | 5.0 |

Average of individual ratios: (0.64 + 5.0) / 2 = **2.82**. Ratio of sums: (64 + 0.5) / (100 + 0.1) = **0.64**.

The first answer is nonsense — it claims 282% of the added water became runoff. It happened because tile-year B has almost no water in it, so its ratio is numerically unstable and dominated by noise, yet the averaging step gives it exactly the same weight as tile-year A, which carries a thousand times more water. The ratio of sums weights each tile-year by the water it actually contributes, which is what "where did the water go" means.

**What are the real numbers?** For the 247,545 snow-addition tile-years:

- runoff **64.3%** (surface runoff 43.1%, baseflow 21.2%)
- ET **35.9%**
- storage **4.2%**
- peatland free-standing water **−2.7%**
- residual **−1.7%**

**Why do these sum to 100% only after the residual is included?** Because the residual is defined as whatever is left over: $\epsilon = I_{\mathrm{snow}} - \Delta R - \Delta ET - \Delta S$. Once you include it, closure is arithmetic, not evidence. The informative part is that the residual is *small* — about 4% of the input — which tells you the accounting is nearly complete.

**What does the negative residual mean?** That the accounted sinks and storage change slightly *exceed* the estimated input: we found about 4% more water leaving than we think went in. Plausible contributors include the monthly-mean storage endpoint approximation (we do not have true October 1 and September 30 states) and the September boundary issue — mean September signed input is 9.01 kg m⁻², 15.5% of the annual total, and much of it appears to melt and redistribute within September itself. The boundary sensitivity test (§S1.5) shows moving to a September–August year changes the residual by only −0.6 percentage points, so the boundary is not the main explanation.

### 2.5 Why use a spatial-block bootstrap?

**What question are we asking?** How much would our partition estimate move if we had happened to sample a different collection of places?

**Why would a simpler approach be risky?** The naive bootstrap resamples individual tile-years. With 247,545 of them, that produces gloriously narrow confidence intervals — and they are fiction. Those 247,545 tile-years are not 247,545 independent experiments. A Colorado mountain tile behaves almost identically to the Colorado mountain tile beside it. Resampling individually treats each as fresh evidence and drives the interval toward zero width.

**What do we actually calculate?** We bin tiles into approximately 5° × 5° geographic blocks. Then:

1. every tile-year in a block stays with its block — they travel together;
2. each bootstrap realization draws the original number of *blocks*, with replacement;
3. both numerator and denominator are recomputed from the resampled block sums, so each replicate is a complete re-estimate of the fraction;
4. 1,000 replicates, fixed seed;
5. the 95% interval is the 2.5th to 97.5th percentile of the replicate distribution.

The mental model: *would our inferred partition look similar if the geographic collection of regions in our sample had come out differently?* Some replicates will contain three copies of the Rockies and no Scandinavia; others the reverse. The spread across those alternative worlds is the uncertainty.

**What does this protect against?** Over-confidence arising from treating spatially correlated tiles as independent, and sensitivity of the result to any single region.

**What does it NOT protect against?** Quite a lot, and this is worth being honest about. It does not protect against a systematic model bias present everywhere — if the Catchment model partitions snowmelt wrongly in all regions equally, every bootstrap replicate inherits that error identically. It does not protect against errors in the increment product itself. It does not protect against the six-year sample being unrepresentative of other decades. It is an uncertainty estimate for *spatial sampling*, nothing more.

### 2.6 What "within-tile anomaly" means

This is the single most important conceptual idea in the controlled regression, so it gets a concrete example.

**What question are we asking?** Does more snow-DA input predict more hydrologic response *at the same place*?

Suppose two tiles:

- **Tile A**, high in the Rockies. Typically receives large snow-DA corrections. Typically generates a lot of runoff.
- **Tile B**, in a dry lowland. Typically receives small snow-DA corrections. Typically generates little runoff.

Now run a naive pooled regression of runoff response on snow input across all tile-years. You will find a strong positive relationship — and it may be entirely an artifact. Tile A is bigger than Tile B on *both* axes, always, in every year. The regression is really just measuring "the Rockies are not the lowlands." It would produce the same positive slope even if, within each tile, snow input and runoff were completely unrelated year to year.

**What do we actually calculate?** For each tile, subtract that tile's own six-year mean from every variable. Tile A's snow input is now expressed as "how much more or less than Tile A usually gets," and Tile A's runoff as "how much more or less than Tile A usually produces."

The regression then asks a much sharper question:

> In years when Tile A receives *more* snow-DA input than is normal for Tile A, does Tile A also produce *more* runoff response than is normal for Tile A?

And it asks the same of Tile B, and of all 48,067 tiles, and pools those within-location relationships.

**What does a positive coefficient mean?** That the relationship holds *within* locations, not merely across them. That is a substantially stronger claim than the pooled version.

**What does it NOT mean?** It still is not causal. Within-tile demeaning removes anything permanently characteristic of a location; it does nothing about a time-varying confounder that moves both snow input and runoff together within a tile. §2.7 and §2.8 address two specific such confounders.

### 2.7 What year fixed effects do

**What question are we asking?** Is the relationship driven by whole-domain year-to-year swings rather than by local variation?

Suppose WY2005 was unusually snowy across the entire Northern Hemisphere. Then in WY2005 nearly every tile has large snow-DA corrections *and* large runoff — not because one caused the other at any location, but because a single hemispheric anomaly drove both everywhere at once. Even after within-tile demeaning, that would still show up as a positive relationship.

**What do we actually calculate?** We include a dummy variable for each water year (five dummies for six years, one being the reference). These absorb whatever is common to a given year across the whole domain.

**What is the coefficient identified from, after this?** Only from *spatially varying departures within each year*, after each tile's climatology has already been removed. In plain terms: in WY2005, some tiles got more snow-DA input than the WY2005 domain-wide average, relative to their own norms. Did *those* tiles show a correspondingly larger runoff response? That is what $\beta$ now measures.

**How would I explain this to a coauthor?** "Tile demeaning removes 'where you are.' Year effects remove 'what year it was.' What's left is 'this place, this year, relative to both.'"

### 2.8 What the OL MAM snow covariate does

**What question are we asking?** Is snow-DA input just a proxy for how much snow the model had anyway?

Here is the concern. Assimilation corrections tend to scale with the thing being corrected. Snow-DA input might simply be larger in years when the underlying modeled snowpack is larger. And larger modeled snowpack independently produces more melt and more runoff. So we might be measuring "big snow years produce big runoff," which is true, uninteresting, and has nothing to do with assimilation.

**What do we actually calculate?** We add the day-weighted March–May mean OL `SNOMASLAND`, expressed as a within-tile anomaly, as a covariate. "Day-weighted" means the March, April, and May monthly means are averaged with weights proportional to days in each month rather than treated as equal — a minor detail, but it is what the production code does. Critically, this is the **OL** snowpack: the model with no assimilation, i.e. the background snow climate at that tile in that year, uncontaminated by the correction whose effect we are trying to measure.

The regression now asks:

> Given two years at the same tile that had the *same* background OL spring snowpack, did the year with more snow-DA input also show more runoff response?

**What does a positive coefficient mean now?** That snow-DA input predicts hydrologic response above and beyond the effect of background snowiness. That is a meaningfully controlled comparison.

**What does it NOT mean?** It is still not causal proof. Three specific caveats:

1. Snow-DA input is not randomly assigned. It is larger where and when the model disagrees with MODIS, which is itself informative about conditions.
2. MAM mean OL snow is one summary of background snowiness. It is not the only one, and the choice is somewhat arbitrary. In practice this turns out not to matter much — substituting the water-year maximum or the March value changes the runoff slope from 0.749 to 0.751 and 0.750 respectively, which is negligible. That insensitivity is reassuring but is not itself proof of correct specification.
3. The regression corroborates the direct accounting. It does not independently establish causality, and the paper says so.

**Note on the manuscript text:** an earlier draft of §2.6 described this covariate as the "OL water-year peak snow mass." That is the sensitivity variant, not the production specification. The production code uses the day-weighted MAM mean. This has been corrected in the manuscript.

### 2.9 Putting the M3 equation into words

$$
\tilde{Y}_{i,t} = \alpha + \beta\, \tilde{I}_{i,t} + \delta\, \tilde{S}^{\mathrm{OL}}_{i,t} + \gamma_t + \varepsilon_{i,t}
$$

Term by term:

| Symbol | In English |
|---|---|
| $\tilde{Y}_{i,t}$ | How much more (or less) runoff response tile $i$ showed in year $t$ than tile $i$ normally shows |
| $\tilde{I}_{i,t}$ | How much more (or less) snow-DA water tile $i$ received in year $t$ than tile $i$ normally receives |
| $\tilde{S}^{\mathrm{OL}}_{i,t}$ | How much more (or less) spring snow the *unassimilated* model had at tile $i$ in year $t$ than normal |
| $\gamma_t$ | Whatever was distinctive about year $t$ across the entire domain |
| $\beta$ | The answer: extra response per extra unit of snow-DA input, holding background snowiness and year constant |
| $\alpha$, $\varepsilon$ | Intercept and unexplained remainder |

The tilde on three of these terms is doing enormous work. Every tilde means "this tile's own long-run mean has been subtracted."

**Why does $\beta$ have units of response per kg m⁻² input?** Because it is a slope: (kg m⁻² of runoff) per (kg m⁻² of input). Both numerator and denominator are water masses, so $\beta$ is dimensionless. A value of 0.749 reads directly as "about three-quarters of a marginal unit of anomalous snow-DA input shows up as anomalous runoff."

**Why do the four closing coefficients sum to exactly 1?** Because runoff, ET, storage, and residual are *defined* to sum to the input:

$$
\Delta R + \Delta ET + \Delta S + \epsilon = I_{\mathrm{snow}}
$$

That identity holds at every tile-year, so it survives demeaning. Then all four responses are regressed on the *same* predictor using the *same* design matrix and the *same* sample. Regression is linear, so the coefficients inherit the identity: regressing the sum of the four responses on $\tilde{I}$ is the same as regressing $\tilde{I}$ on itself, which gives exactly 1.

**What does this NOT mean?** The summing-to-one is a construction, not a validation. It would hold even if every individual coefficient were badly estimated. Do not present it as evidence that the analysis is correct — it is a consistency check on the arithmetic, nothing more. In the production run: 0.749 + 0.182 + 0.085 − 0.016 = 1.000.

### 2.10 Why 0.749 is not supposed to equal 0.643

This comes up every time someone reads the two numbers side by side, so it needs a clear answer.

**The direct accounting says 64.3% runoff. The controlled regression says 0.749. Is one of them wrong?**

No. They answer different questions.

**64.3% is an average fate.** Take all the water that snow assimilation added, across every positive-input tile-year in six years, and ask where it ended up in aggregate. This is a mass-accounting statement about a specific pool of water, conditional on the tile-years where assimilation added water.

**0.749 is a marginal response.** Take one *additional* unit of anomalous snow-DA input at a given location in a given year — over and above what that location normally gets, in a year that was otherwise average across the domain, at a location whose background OL snowpack was also average — and ask how much additional runoff appears. This uses signed within-tile interannual variation, includes all 288,402 complete signed tile-years (not only positive-input ones), and removes tile means, year effects, and background snowiness.

These differ in the sample used, the contrast that identifies the effect, and the conditioning. There is no reason for them to be numerically equal, and it would be somewhat suspicious if they were.

**So what is the actual corroboration?** That both independently say the same *qualitative* thing: runoff is the dominant fate, ET is second, storage is small. Two methods with different samples, different weighting, and different confounding structures agree on the physical story. That is meaningful.

**What does it NOT mean?** It does not mean we have two estimates of one parameter that we could average, or that we should reconcile the discrepancy. The paper must not present them as competing estimates of the same quantity.

**How would I explain this to a reviewer?** "The 64% is where the water went. The 0.75 is where the next unit would go. Those are different questions and we report them as different questions."

---

## 3. Section 2.7 — long-term trends

### 3.1 Why analyze DA−OL directly?

**What question are we asking?** Is the *influence of assimilation* changing systematically over 24 years?

**Why would a simpler approach be risky?** The intuitive approach is to fit a trend to DA, fit a trend to OL, and subtract. This is worse for several reasons.

Consider the common-mode problem. Suppose the precipitation forcing gets wetter over 24 years. Then OL gets wetter and DA gets wetter, both for the same reason, which has nothing to do with assimilation. If you fit each separately you get two positive trends, and you are relying on them to cancel in the subtraction. But they were estimated from two different noisy trajectories, each with its own robust-median fitting behaviour, so they do not cancel cleanly.

**What do we actually calculate?** We form the monthly paired difference first:

$$
D_i(t) = \mathrm{DA}_i(t) - \mathrm{OL}_i(t)
$$

and fit the trend to $D_i(t)$.

Now the cancellation happens month by month, in the data, before any fitting. A wet month is wet in both runs, and the difference is unaffected. What survives is the assimilation increment's effect. The trend of $D$ asks precisely the intended question: *is the effect of assimilation itself becoming systematically more positive or negative through time?*

There is a second benefit: identical sample support. Missing OL or DA values are cross-masked before differencing, so we are guaranteed to compare like with like every month. Two separately fitted trends could quietly rest on slightly different month sets.

**What does a significant positive DA−OL trend mean?** That assimilation is progressively pushing the state in one direction relative to the free-running model. **What does it NOT mean?** It does not mean the state itself is trending — DA and OL could both be flat while their difference grows.

OL and DA trends are still computed and reported, as context and control fields. The precipitation control is the clearest use: 3,719 OL and 3,726 DA tiles have significant precipitation trends, with OL–DA slope correlation 0.9998, and *no* paired DA−OL precipitation trend survives field significance. That is exactly what a common-forcing null should look like, and it validates the whole approach.

### 3.2 Why remove seasonality?

**What question are we asking?** Is there a multi-decadal tendency, ignoring the annual cycle?

Monthly RZMC, SCF, ET and similar fields have enormous annual cycles. Snow cover in January versus July is a far bigger swing than any 24-year trend. If you fit a trend to the raw series, the seasonal cycle dominates the residual variance and swamps the signal.

**Why would a simpler approach be risky?** The obvious fix is to subtract each calendar month's mean. This subtly damages the trend, and it is worth seeing why.

Suppose a series rises perfectly linearly with no seasonality at all. Compute the mean of all Januaries: those are spread evenly through the record, so their mean sits near the record midpoint. Same for February, and so on — every calendar-month mean lands near the same value. Subtracting them removes the *level* but the trend survives. So far, fine.

Now suppose the record does not contain a whole number of years, or has missing months, or the seasonal cycle interacts with the trend. Then each calendar month's mean sits at a slightly different position along the trend line, and the twelve means differ from each other by trend-induced amounts rather than real seasonality. Subtracting them removes part of the trend and converts a smooth line into an annual sawtooth.

**What do we actually calculate?** Four steps:

1. estimate a preliminary exact Theil–Sen slope on the raw series;
2. temporarily remove that slope;
3. compute each calendar month's climatology from the *detrended* values;
4. subtract only that seasonal climatology from the **original** series.

The trend is removed at step 2 solely so the climatology is not contaminated by it, and is restored at step 4.

**What does this NOT do?** It does not remove the trend before testing. This is the point people get wrong. After step 4 the series still has its full long-term tendency; only the seasonal cycle has been taken out.

The same adjusted paired tile series feed both main temporal summaries. Figure 16 estimates a whole-record slope at each tile. Figure 17 first area-weights those adjusted tile series within each region and then calculates P1–P9 means. Seasonality is therefore removed once, at the tile level, rather than being estimated again from the regional aggregate.

### 3.3 What Theil–Sen does

**What question are we asking?** How big is the trend?

**What do we actually calculate?** Take every pair of points in the series. For each pair, compute the slope of the line joining them. Take the median of all those pairwise slopes. That is the Theil–Sen estimate.

A tiny example. Values 1, 2, 3, 20 at times 1, 2, 3, 4. Pairwise slopes: (2−1)/1 = 1, (3−1)/2 = 1, (20−1)/3 = 6.33, (3−2)/1 = 1, (20−2)/2 = 9, (20−3)/1 = 17. Median of {1, 1, 1, 6.33, 9, 17} = **3.67**. Ordinary least squares on the same data gives about 5.3, pulled hard by the single value of 20.

**Why use it?** Robust to outliers, and makes no Gaussian-residual assumption. Monthly land fields have heavy tails, occasional extreme months, and — for increment diagnostics — large numbers of exact zeros. A median-based estimator handles all of that gracefully.

**What does it NOT do?** Three things, all important:

- It does not provide a significance test. It is an estimator only.
- It does not solve autocorrelation. The pairwise slopes are not independent.
- It does not address multiple testing.

Those are the jobs of §3.4, §3.5, and §3.6 respectively. Theil–Sen only tells you *how big*.

### 3.4 What Mann–Kendall does

**What question are we asking?** Is there a monotonic tendency at all, or could this be noise?

**What do we actually calculate?** For every pair of points, record only whether the later value is higher or lower than the earlier one — a $+1$ or a $-1$. Sum those signs. Under the null hypothesis of no trend, ups and downs should roughly cancel and the sum should be near zero. A large positive sum indicates a persistent upward tendency.

Note what is *not* used: the sizes of the changes. Mann–Kendall works on ranks and directions only. This is why it needs no distributional assumption.

**How does this differ from Theil–Sen?**

| | Theil–Sen | Mann–Kendall |
|---|---|---|
| Question | How big? | Real or noise? |
| Output | A slope, in units per year | A test statistic and p-value |
| Uses magnitudes? | Yes (median of slopes) | No (signs only) |

They are complementary and always reported together. Theil–Sen without Mann–Kendall gives you a number with no idea whether to believe it; Mann–Kendall without Theil–Sen tells you something is happening but not how much.

### 3.5 Why autocorrelation matters

**What question are we asking?** How much independent evidence do 288 monthly values actually contain?

Land states have memory. A wet month is usually followed by a wet month, because soil moisture drains slowly and the atmosphere has persistent regimes. So consecutive months are not two independent pieces of information — they are closer to one and a bit.

**Why would ignoring this be risky?** Because Mann–Kendall's null distribution assumes independence. Under that assumption, a long run of increases is very unlikely. But under strong persistence, long runs happen *routinely* by chance. A persistent random walk wanders convincingly upward for years at a time while having no trend whatsoever. If you test it assuming independence, you find "significant trends" everywhere.

The synthetic calibration shows exactly this. On AR(1) = 0.8 series with genuinely no trend, the pointwise modified Mann–Kendall still flagged 12% of series — and an *unmodified* test would flag many more. This is not a hypothetical risk.

**What do we actually calculate?** The Hamed–Rao correction, in a deliberately conservative form:

- estimate rank autocorrelation from the Sen-detrended residuals at actual monthly lags 1 through 24;
- require at least 24 pairs per lag, and assess significance at the 0.05 level;
- inflate the test variance based on those correlations;
- **only significantly positive** correlations contribute to the inflation;
- the variance factor is **floored at 1.0**.

Those last two rules make the correction one-directional. Positive persistence makes significance harder to achieve. Negative autocorrelation — which would technically justify making the test *easier* — is not permitted to do so. By construction the correction can only ever make a result less significant.

**What does a significant trend mean after this correction?** That the tendency is distinguishable from noise even after accounting for the series' own memory. **What does non-significance mean?** Only that we could not resolve it. Given the deliberate conservatism, non-significance is weak evidence of absence.

### 3.6 What BH FDR actually means

This section exists because false-discovery rate is very commonly misdescribed, including in published papers.

**Start with the p-value.** A p-value is: *if there were truly no effect, how often would I see a result at least this extreme?* p = 0.03 means such a result would arise 3% of the time under the null. It is a statement about data-under-a-hypothesis. It is emphatically **not** the probability that the hypothesis is true.

**Now the multiplicity problem.** We test 112,573 tiles. If nothing is happening anywhere and we flag p < 0.05, we expect about 5,600 false flags. Because neighbouring tiles are correlated, those false flags cluster into coherent-looking patches, and a map of them will look like a real geophysical pattern. This is the single easiest way to produce a convincing but entirely spurious trend map.

**What Benjamini–Hochberg does.** Rank all p-values from smallest to largest and apply a threshold that grows with rank, so that the *expected proportion of false positives among everything you flag* is held at 5%. If BH flags 7,892 tiles at q = 0.05, we expect roughly 5% of that flagged set — about 395 tiles — to be false discoveries.

**Two misconceptions to correct explicitly:**

> ❌ "q = 0.009 means there is a 0.9% probability the null hypothesis is true."

No. Like a p-value, q is not a posterior probability of any hypothesis. It is a property of the selection procedure applied to a collection of tests.

> ❌ "FDR at 5% means each significant tile has a 5% chance of being a false positive."

No. The control is on the *expected proportion across the flagged set as a whole*, not on any individual tile. Individual tiles vary enormously in their evidence — a tile that squeaked past the threshold is far more likely to be false than one with an overwhelming signal. FDR says nothing about which is which.

**The correct phrasing** is something like: "among the discoveries produced by this procedure, the expected false-discovery proportion is controlled at 5% under the procedure's assumptions."

**How would I explain this to a coauthor?** "It's a quality guarantee on the batch, not a warranty on each item."

**Family structure in this paper.** Each complete mapped field is its own FDR family. OL, DA, and DA−OL are separate families even when they share a sample. FDR must be computed over the complete field — a subset run gets a different threshold and is explicitly labeled non-global.

### 3.7 Why a confidence interval and FDR significance can disagree

**What question are we asking?** Why does a tile sometimes show FDR significance while its confidence interval contains zero?

Because they are two related but mathematically distinct procedures, and neither is derived from the other:

- **FDR significance** comes from the Hamed–Rao modified Mann–Kendall p-value, put through Benjamini–Hochberg across the field. It is rank-based.
- **The reported interval** is SciPy's nominal Theil–Sen interval with its half-width expanded by the square root of the Hamed–Rao variance factor. It is a first-order adjustment to a magnitude-based interval.

That expansion is a reasonable approximation, but it is *not* the exact inversion of the modified Mann–Kendall test. Where the two procedures behave differently — particularly for tied and zero-heavy series — they can disagree.

**Where does this actually happen?** Mostly in the increment activity diagnostics, which have many exact zeros and many ties. In the Phase 1 batch:

| Field | FDR significant | CI disagreements | Fraction |
|---|---:|---:|---:|
| TSOIL1 DA−OL | 8,143 | 2 | 0.03% |
| RZMC correction mean | 3,063 | 3 | 0.10% |
| SFMC correction RMS | 102,004 | 12,985 | 12.73% |
| Absolute soil-water activity | 106,009 | 23,454 | 22.12% |

Notice the pattern: the paired physical state fields agree almost perfectly, and disagreement is concentrated in the zero-heavy activity diagnostics. Additionally, 6,203 significant absolute soil-water-activity tiles have an *exact zero* median Sen slope — the rank test detects a tendency but the median pairwise slope is precisely zero because so many pairwise slopes are zero.

**The rule for this paper:** maps use FDR significance. The first-order interval is not plotted as inferential whiskers for the activity fields.

**What does this NOT mean?** It is not a software bug, and it must not be "fixed" by forcing one product to mimic the other. It is a genuine and expected property of using two different inferential procedures. The production output exposes an explicit `fdr_ci_disagreement` flag so the disagreement is visible rather than hidden. If interval-based mapped inference is ever needed, the right answer is a test-inverted or block-bootstrap interval, not a fudge.

---

## 4. Section 2.7 — regional observing-system periods

### 4.1 How Figure 17 reuses the Figure 16 data

**What question are we asking?** When do the regional DA−OL differences seen in the whole-record trend map emerge?

At every tile we already have the paired, trend-retaining, seasonally adjusted DA−OL series used by Figure 16. For Figure 17 we area-weight those same adjusted tile series within six fixed regions, producing one monthly regional series per region. We then calculate its mean separately in P1–P9.

This is simpler than fitting another seasonal model after aggregation. It also lets every tile retain its own calendar-month climatology. The only difference between the figures is the final summary:

| Figure | Shared input | Final summary |
|---|---|---|
| Figure 16 | Seasonally adjusted paired tile DA−OL | Whole-record Theil–Sen slope at each tile |
| Figure 17 | Seasonally adjusted paired tile DA−OL | Area-weighted regional P1–P9 means |

### 4.2 What an adjacent-period difference means

For P6, for example, the estimate is simply

$$
\overline{D}_{P6}-\overline{D}_{P5}.
$$

It says that average assimilation influence differs between the two observing-system periods. It does **not** say the series jumped instantaneously on the boundary date. A gradual change within either period can produce the same difference.

The six regions are tested together at each boundary with BH FDR. This gives eight separate six-test families, one for each P2−P1 through P9−P8 comparison. The regional boxes were motivated by Figure 16, so Figure 17 is an exploratory timing decomposition of an existing spatial result, not an independent discovery that those regions are special.

### 4.3 Why the effective sample size is small

The regional monthly series retain substantial memory. We estimate one lag-1 autocorrelation from residuals about the period means and use

$$
n_{\mathrm{eff}}=n\frac{1-\rho}{1+\rho}.
$$

With regional $\rho$ between 0.63 and 0.80, a period containing many calendar months contains far fewer independent pieces of information. The effective-sample-size intervals are conservative: they produce the same nine RZMC FDR decisions as a fitted-AR(1) parametric bootstrap, while a moving-block bootstrap adds only three marginal cases.

### 4.4 What changed when we unified the seasonal adjustment

Previously the regional aggregate received a separate least-squares seasonal adjustment. Recomputing Figure 17 from the tile-adjusted series changed no RZMC or SFMC FDR decision. The largest adjacent-period estimate changed by only $8.1\times10^{-6}$ m³ m⁻³. The value of the change is therefore conceptual clarity and reproducibility, not a different scientific result.

---

## Appendix A — superseded segmented-transition workflow

The material below documents an earlier analysis that is no longer used in the manuscript. It is retained only to explain archived outputs; it must not be used to describe the current Figure 17 analysis.

### A.1 Trend versus level change versus slope change

Three different shapes, three different coefficients, constantly confused.

**Trend.** The series rises steadily across the whole record. Draw one straight sloping line from 2000 to 2024.

**Level change.** The series is flat, then at one instant jumps to a new baseline and is flat again. Draw a flat line, a vertical step at April 2015, then another flat line higher up.

**Slope change.** The series is flat, then from one date onward begins climbing. Draw a flat line that bends upward at April 2015 — no jump, just a change of direction.

**Why does this distinction matter so much here?** Because it resolves what otherwise looks like a contradiction in our results. We report a large, highly significant April 2015 jump in DA−OL storage (+2.144 kg m⁻², q = 0.009) *and* no significant whole-record storage trend. Both are true. A single step up, followed by a flat period, is not a monotonic 24-year tendency. A Theil–Sen slope fitted through a step function returns something small and poorly resolved, because most pairwise slopes within each flat segment are near zero.

**How would I explain this to a coauthor?** "It's a staircase, not a ramp. Asking for the trend of a staircase gives you an answer, but not a useful one."

This is exactly why the paper runs the transition analysis separately from the trend analysis rather than treating one as a summary of the other.

### A.2 What segmented regression is doing

**What question are we asking?** At each known observing-system boundary, does the series step or bend?

**What do we actually calculate?** One single regression fitted to the whole 288-month record, containing 28 parameters:

| Component | Count | What it absorbs |
|---|---:|---|
| Intercept | 1 | Overall level |
| P1 baseline slope | 1 | The starting trend |
| Calendar-month effects | 11 | The seasonal cycle |
| Level changes at P2–P9 | 8 | Jumps at each boundary |
| Slope hinges at P2–P6, P8, P9 | 7 | Bends at boundaries long enough to support one |

In plain English: *we fit one continuous description of the entire record that contains seasonality, a baseline trend, an allowed jump at each known boundary, and an allowed change of slope wherever the adjacent periods are long enough to estimate it.*

**Why fit all boundaries simultaneously rather than one at a time?** Because each transition is then estimated while accounting for all the others. If you tested April 2015 in isolation, a real change at June 2007 would sit unmodelled in the residuals, distorting both the fitted baseline and the apparent size of the 2015 step. Fitting jointly means the P6 coefficient answers "how much does the series jump at P6, *given* everything else the model already knows about," which is the question we want.

**What assumptions remain?** That the boundaries are the right dates (they come from the observing-system registry, fixed in advance); that the underlying form really is piecewise-linear-plus-seasonal; and that a single AR(1) adequately describes the residual dependence.

### A.3 Why P7 gets no slope hinge

**What question are we asking?** Can a 15-month period support an independent slope estimate?

P7 — the CYGNSS period, August 2018 to November 2019 — is 15 months long. Our declared minimum reliable segment length is 24 months.

Fitting a slope to 15 monthly points that also contain a seasonal cycle and are strongly autocorrelated would produce a number, but a meaningless one: it would be dominated by whatever the annual cycle and short-term persistence happened to be doing across that particular window.

So P7 receives a **level change** — we do ask whether the series steps at August 2018 — but **no independent slope hinge**. The P6 slope simply continues until the P8 boundary.

**The crucial point:** this was decided from the period registry *before* looking at any results. It is a design constraint, not a post-hoc reaction to an inconvenient estimate. P7 is also exempt from PELT agreement scoring for the same structural reason — a 24-month minimum segment cannot resolve a 15-month period.

**What does it NOT mean?** It does not mean nothing happened during P7. It means our design deliberately declines to make a slope claim it cannot support.

### A.4 Prais–Winsten in plain language

**What question are we asking?** How do we fit a regression whose residuals remember the previous month?

Ordinary regression assumes residuals are independent. Monthly land-surface residuals are not: if the model is above the data this month, it is probably above the data next month too. Ignoring this does not usually bias the coefficients much, but it badly understates their uncertainty.

**What do we actually calculate?** Iteratively:

1. fit the regression;
2. look at the residuals and estimate how strongly each depends on the previous one — the AR(1) coefficient;
3. transform the data to remove that dependence;
4. refit;
5. repeat until the AR(1) estimate stops moving.

Production settings: at most 25 iterations, convergence tolerance 10⁻⁷, AR(1) magnitude capped at 0.98. The first observation is retained rather than discarded, which is the specific feature distinguishing Prais–Winsten from the simpler Cochrane–Orcutt procedure — with only 288 points, throwing one away for no reason is wasteful.

**What does Prais–Winsten NOT do here?** It does not produce our final uncertainty estimates. This is worth stating clearly because it is the natural assumption. We tried using its standard errors — both OLS-with-Newey–West and Prais–Winsten-with-Newey–West — and *both were anti-conservative* in the fixed-seed 288-month AR(1) no-transition test. They declared too many transitions significant when none existed.

So Prais–Winsten's historical job was narrower: it supplied the fitted serial-dependence model that the bootstrap in Appendix A.5 then used. The Newey–West standard errors remain in archived output as diagnostics and play no role in current manuscript inference.

### A.5 Innovation-resampling bootstrap

**What question are we asking?** How much would our transition coefficients move if the noise had come out differently?

**What do we actually calculate?**

1. Fit the complete segmented Prais–Winsten model. This splits the series into a **fitted deterministic part** (seasonality, trends, jumps) and **residuals**.
2. Whiten the residuals into innovations — the genuinely fresh noise arriving each month, once the AR(1) memory is accounted for.
3. Center the innovations.
4. Resample them with replacement — a new random ordering of the same pool of shocks.
5. Propagate them through an AR(1) process to rebuild a noise series with realistic memory, discarding a 120-month burn-in so the simulated series is not contaminated by its arbitrary starting value.
6. Add that synthetic noise to the fitted deterministic model, producing an alternative plausible version of the record.
7. **Refit the entire 28-parameter segmented model** to it.
8. Repeat 1,999 times.

The spread of the P6 coefficient across those 1,999 refits is its bootstrap distribution.

**Why must the whole model be refitted each time?** Because we want the uncertainty in the P6 coefficient *as it is actually estimated* — jointly with 27 other parameters, some of which are correlated with it. If we perturbed the P6 coefficient alone while holding the rest fixed, we would ignore the fact that in a different realization the fitted seasonality and neighbouring boundary terms would also shift, partly absorbing or amplifying the apparent P6 change. Refitting captures that.

**Why is the simulation AR(1) set to the upper 95% bound?** This is a deliberate conservatism. A model with 28 free parameters is flexible enough to absorb some genuine persistence into its deterministic part — the trends and jumps soak up wandering that really belongs to the noise. The residuals therefore look *less* autocorrelated than the truth, and the fitted AR(1) is biased low. If we simulated at that biased-low value, our bootstrap records would be too well-behaved and our intervals too narrow.

Simulating at the upper 95% confidence bound instead (capped at 0.98) means we are asking: *if persistence were at the pessimistic end of what the data permit, how uncertain would the coefficient be?* Both the fitted and simulation values are reported.

**How would I explain this to a coauthor?** "We're deliberately assuming the noise is stickier than it looks, because we know our model has an incentive to make the noise look tidy."

Inference uses centered two-sided empirical bootstrap p-values and basic bootstrap intervals.

### A.6 Why bootstrap CI and FDR can produce the P6 ET result

This is the best worked example in the paper of pointwise versus family-level inference, and it is worth understanding thoroughly because it is the one result most likely to be challenged.

**The numbers.** At P6 (April 2015), for DA−OL evapotranspiration:

- estimate **+0.502 kg m⁻² month⁻¹**
- 95% bootstrap interval **0.094 to 0.883** — excludes zero
- **q = 0.072** — does not meet q < 0.05

**Is this a contradiction?** No.

The bootstrap interval answers a *pointwise* question: considering ET on its own, is the P6 coefficient distinguishable from zero? Answer: yes, though not by a wide margin — the lower bound of 0.094 is not far above zero.

The q-value answers a *family* question: across all nine series tested at this boundary, if we flag everything at least this extreme, what proportion of flagged results would we expect to be false? At 0.072, we cannot hold that expected proportion at 5%.

Both are correct. The interval says "considered alone, this looks real." The q-value says "considered alongside the eight other things we tested at the same boundary, it is not distinctive enough."

**The correct paper language:** *convergent and suggestive positive response; not formally FDR-significant.* Never "significant." Never "no effect."

**The contrast that makes it clear.** At the same boundary, DA−OL total storage:

- estimate **+2.144 kg m⁻²**
- interval **0.932 to 3.266**
- **q = 0.009** — survives comfortably

Same boundary, same family, same method. Storage clears the family bar; ET does not. That is what family correction is for.

**One further point that strengthens the ET case.** PELT, given no dates at all, independently places an ET break exactly in April 2015. And the increment budget across P6 closes to about 5%, requiring an ET response of roughly this magnitude. So there are three convergent lines of evidence — pointwise interval, independent detection, budget closure — and one formal test the result does not pass. The honest report includes all four facts.

---

## Appendix B — superseded independent PELT analysis

This blind-changepoint analysis is likewise historical and is not part of the current manuscript inference.

### B.1 What PELT is asking

**What question are we asking?** Ignoring everything we know about satellites, where does this record appear to break?

No observing-system date is supplied to the algorithm. PELT searches for points at which splitting the record into separate segments substantially improves the statistical description of the data, while charging a penalty for each additional break.

**Why is this worth doing when we already have the known-date test?** Because it can fail in ways the known-date test cannot, and succeed in ways it cannot. It can find breaks we did not anticipate. It can decline to find a break at a date we were confident about. And when we specify April 2015 in advance for physical reasons and an algorithm that has never heard of SMAP independently picks April 2015 out of 288 months, that agreement carries real evidential weight.

### B.2 Why a penalty is necessary

**What question are we asking?** How do we stop the algorithm calling every wiggle a changepoint?

Without a penalty, more breaks always fit better. A separate segment for every single month would fit perfectly and describe nothing. So each break must pay a price, and a break is only accepted if it improves the fit by more than it costs.

The penalty basis is BIC — the Bayesian Information Criterion — which charges roughly (number of parameters) × log(number of observations). The intuition without the derivation: extra parameters must earn their place, and the bar rises with sample size, because in a long record even pure noise offers many superficially attractive places to split.

Higher penalty → fewer, stronger breaks. Lower penalty → more, weaker breaks. Our primary penalty is **1.25 × BIC**, slightly stricter than standard.

### B.3 Why we require consensus

**What question are we asking?** Is this break real, or an artifact of one convenient tuning choice?

Any single penalty and any single cost function might produce an appealing answer by luck. So a break is accepted only if it satisfies **all three** conditions:

1. it appears in the primary method (piecewise AR(1) + linear cost) at 1.25 × BIC;
2. it recurs within three months in at least **half** the penalty grid — multipliers 0.5, 0.75, 1.0, 1.25, 1.5, 2.0 — so it is not an artifact of one penalty;
3. an independently prewhitened linear formulation also finds a break within three months — so it is not an artifact of one treatment of autocorrelation.

Every segment must additionally contain at least 24 months.

This is intentionally demanding, and the numbers show it: across 43 series only 37 breaks are accepted. Twenty match a known boundary within ±3 months, two more within ±6 months, and **fifteen remain unmatched**.

**Why do we keep the unmatched ones?** Because discarding them would be cheating. Reporting only the breaks that landed near a satellite date, while quietly dropping fifteen that did not, would badly overstate the correspondence between the record and the observing-system chronology. The unmatched breaks are retained and explicitly not attributed.

### B.4 What the synthetic tests tell us

We planted known effects in synthetic AR(1) series and asked how often the pipeline recovered them.

**Where it does well:**

- seasonal white-noise null: **zero** accepted breaks;
- AR(1) = 0.8 null: **zero** accepted breaks;
- isolated abrupt level shift: **91.7%** recovery within six months;
- two opposing level shifts: **85.4%**;
- level shift near the minimum-segment edge: **91.7%**.

**Where it does badly:**

- gradual slope hinge: **4.2%** recovery within six months, with a median localization error of about **24.5 months**.

**What follows from this?**

A PELT break near P6 is strong corroboration — the detector essentially never fires on null series, so when it fires the evidence is real.

A PELT *absence* near a known date is nearly uninformative about gradual changes. At 4.2% recovery, failing to detect a slope hinge is the expected outcome even when one is definitely present. This is exactly why known-date slope inference remains the job of the interrupted-series model, and why the changepoint validation gate does not require gradual-hinge recovery.

**How would I explain this to a reviewer?** "The detector is a smoke alarm for steps. It's excellent evidence when it goes off. Its silence tells you almost nothing about a slow temperature rise."

---

## 6. Worked examples from this paper

### 6.1 RZMC regional trend

**The result.** Significant paired DA−OL root-zone soil-moisture trends occur at 7.01% of valid-land tiles (7,892 tiles), of which 7,267 are positive and 625 negative. The area with significant trends rises from 4.36% in OL to 9.18% in DA. Yet the area-weighted OL, DA, and DA−OL domain-mean trends are all non-significant.

**Is that a contradiction?** No. Regional effects of opposite sign cancel in a domain mean, and even same-sign regional effects covering 7% of land dilute heavily against the other 93%.

**The lesson.** Regional modification is not global wetting. The correct statement is that assimilation modifies the regional evolution of root-zone soil moisture without producing a resolved global land-mean trend. Never collapse the two into "DA makes the land wetter over time."

### 6.2 April 2015 RZMC transition

**The result.** The P6−P5 RZMC DA−OL difference is resolved globally and in Australia, southern Africa, and North Africa/Middle East. The estimates are +1.52, +3.46, +3.57, and +1.62 ×10⁻³ m³ m⁻³, respectively. Western North America and northern Eurasia are not resolved at P6.

The P2−P1 pattern is different: it is resolved globally and in western North America and northern Eurasia, but not in the three warm, dry regions. Thus the timing of regional change differs geographically rather than following one global sequence.

**The lesson.** Say that P6−P5 coincides with the introduction of SMAP brightness-temperature assimilation and a broader multi-sensor reorganization. Do not call it an instantaneous April 2015 jump, and do not claim that SMAP is the unique cause.

### 6.3 Snow-water partition

**The result.** Direct accounting: 64.3% runoff, 35.9% ET, 4.2% storage, −2.7% peatland free-standing water, −1.7% residual. Controlled regression: 0.749, 0.182, 0.085, −0.017, 0.0007.

**What is the corroboration?** That two methods — with different samples (247,545 positive-input tile-years versus 288,402 signed tile-years), different weighting, different identifying variation, and different confounders — independently agree that runoff dominates, ET is second, and storage is small.

**What is it NOT?** Two estimates of the same parameter. See §2.10. Do not average them, do not reconcile them, and do not describe the difference as a discrepancy.

---

## 7. Common statistical traps in this manuscript

A checklist. Each of these is a mistake it would be easy to make in this specific paper.

**Correlation ≠ causation.** Every transition result is temporal coincidence with a known date. Write "coincides with," not "caused by."

**Significance ≠ importance.** The global RZMC P6−P5 difference is +0.00152 m³ m⁻³ — resolved but physically small. Report both the significance and the magnitude.

**Non-significance ≠ zero effect.** Write "not resolved," not "no change." Short periods and strong monthly persistence make regional period differences difficult to resolve.

**Global mean ≠ regional behaviour.** RZMC: 7% of tiles significant, domain mean not significant. Both true.

**Period-mean difference ≠ instantaneous level change.** A P6−P5 contrast can arise from gradual evolution within either period; it is not evidence of a jump on 1 April 2015.

**Trend of difference ≠ difference of trends.** We fit the trend of DA−OL. Do not describe it as the difference between the DA trend and the OL trend.

**Pointwise p ≠ FDR q.** Regional tests are corrected in six-test families at each boundary. A small unadjusted p-value does not by itself define a resolved transition.

**CI excludes zero ≠ automatically FDR significant.** Also the reverse: for the zero-heavy activity diagnostics, FDR significance with a zero-containing interval occurs at up to 22% of significant tiles.

**Thousands of neighbouring tiles ≠ thousands of independent experiments.** This is why we block-bootstrap and why we use FDR rather than pointwise stippling.

**Physical consistency ≠ independent validation.** The P6 budget closing to 5% shows the terms are mutually consistent. It uses the same model output as everything else and is not an external check against observations.

---

## 8. One-page cheat sheet

| Method | Question | Why we use it | Main caution | Example in paper |
|---|---|---|---|---|
| Direct mass accounting | Where did added water go? | Mass is conserved, so the budget can be checked for closure | Ratio of sums, not mean of ratios; residual must be reported | 64.3% of snow-DA water becomes runoff |
| Spatial-block bootstrap | How much would this move with different regions sampled? | Neighbouring tiles are not independent | Covers spatial sampling only, not model bias | 61.1–67.4% interval on the runoff fraction |
| Within-tile regression | Does anomalous input predict anomalous response at one place? | Removes permanent geographic contrasts | Corroboration, not causal proof | β_runoff = 0.749 |
| Theil–Sen | How big is the trend? | Robust to outliers and non-Gaussian noise | Estimator only — no significance, no autocorrelation fix | SCF declines 0.000554 yr⁻¹ in OL |
| Modified Mann–Kendall | Is the trend distinguishable from noise? | Rank-based, no distributional assumption; Hamed–Rao handles memory | Deliberately conservative; can only reduce significance | 12% → 0% false positives under AR(1)=0.8 |
| BH FDR | How many of my flagged tiles are false? | Thousands of tests need multiplicity control | Controls the batch, not the individual tile | 7,892 significant RZMC DA−OL tiles |
| Regional period means | When do mapped RZMC differences emerge? | Reuses the same adjusted tile series as the trend map | Difference between period averages, not a boundary jump | P6−P5 is resolved in four warm/dry regions plus global land |
| AR(1) effective sample size | How much independent monthly information is present? | Regional soil moisture has strong memory | Conservative approximation, checked against two bootstraps | 9 of 48 RZMC comparisons resolved |

### The 30-second explanation of our statistical strategy

For the snow-water analysis, the primary result is direct mass accounting, with uncertainty estimated by resampling spatial blocks; a within-tile controlled regression independently corroborates the partition pattern after removing persistent geographic differences, common year effects, and variations in background snowpack. For temporal consistency, we form paired DA-minus-OL tile series, remove seasonality once with a trend-preserving Theil–Sen adjustment, and use those same adjusted tile series in both figures. Figure 16 estimates tile trends with autocorrelation-corrected Mann–Kendall inference and spatial FDR control; Figure 17 area-weights the tile series regionally, compares adjacent P1–P9 means with AR(1)-aware uncertainty, and controls FDR across six regions at each boundary.

### If you only remember one thing per section

- **§2.6:** We followed conserved mass, so the budget closes and the closure is checkable.
- **§2.7 trends:** We fit the trend of the *difference*, and we corrected for both memory and multiplicity, in both cases in the conservative direction.
- **§2.7 periods:** Figure 17 reuses the adjusted Figure 16 tile series and asks when their regional averages differ between observing-system periods; it does not estimate instantaneous jumps.
- **Everywhere:** Where we found nothing, that means *not resolved*, not *no effect*.
