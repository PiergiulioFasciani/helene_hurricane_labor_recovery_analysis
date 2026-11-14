# Methodology

This document describes the **active implementation**, including the existing 14 November 2026 revision. The repository reorganization does not introduce a new estimator, change any numerical result, or reopen model selection. Historical implementations are stored under `archive/`; they are not executed by the current pipeline.

## 1. Research question and units

The project compares labor-market developments in 39 frozen Hurricane Helene-affected North Carolina counties with a geographically restricted synthetic comparison. The treated aggregate is an **equal-county mean**, not a population-weighted state aggregate. The observational treatment definition and geographic exclusions are fixed in reference rosters and source code.

The main panel contains 3,221 counties over January 2022–September 2026: 57 months and 183,597 county-month rows. This is a balanced calendar, not a claim that every outcome is observed every month. The main matched labor analysis ends in March 2026 because of QCEW coverage. Extended LAUS results continue through September 2026.

## 2. Sources and construction

| Input | Construction and role | Active implementation |
|---|---|---|
| NOAA nClimGrid-Daily county precipitation | Sum September 25–30, 2024; require all six daily values. NOAA state identifiers are mapped to Census FIPS; county-name checks accompany geography matching. | `src/preprocessing/01_noaa.do` |
| Census 2024 county PEP | Population denominator and county geography. | `src/preprocessing/02_census.do` |
| ACS 2023 five-year estimates | Household count, median home value and household income fields. | `src/preprocessing/02_census.do` |
| FEMA disasters 4827–4832 | Registration-level verified real- and personal-property losses, aggregated by county. These are verified applicant losses, not total economic storm damage. | `src/preprocessing/03_fema.do` |
| LAUS | Resident employment, unemployment, labor force and unemployment rate. Default is the explicitly identified economicdata.dev distribution of BLS series; direct official flat-file mode is optional. | `src/preprocessing/05_labor.do` |
| QCEW | Official BLS county industry-10 CSVs; monthly workplace jobs from quarterly files, through 2026 Q1. | `src/preprocessing/05_labor.do` |
| Census Vintage 2025 annual components | July 1 population and annual migration components, constructed separately from the monthly labor panel. | `src/migration/01_census_v2025.do` |

Cached raw files live under `data/raw/`. Intermediate transformations live under `data/interim/`. The full run always constructs its own panel and copies the validated result to `data/input/Helene_NC_ES_preinference_panel.dta`. It does not use an old inferred panel as a substitute for source construction.

Missing source values remain missing. October 2025 LAUS outcomes are deliberately excluded from observed labor outcomes; QCEW is not altered merely because LAUS is unavailable. Source audit tables record coverage and excluded source values.

### Exposure measures

County rainfall is the six-day total in millimetres. FEMA losses are combined as real-property plus personal-property verified loss. Derived descriptive measures include loss per resident and a housing-value proxy equal to household count times median home value.

For counties with complete rainfall and economic-burden information, the code constructs percentile-rank scores and the composite HDSI:

`HDSI = 0.4 × rainfall rank score + 0.6 × economic-burden rank score`.

Economic burden is 100 times verified loss divided by the housing-value proxy. Rank scores run from 0 to 100 using the relevant complete-case sample. Missing FEMA records are not silently equated to zero loss. The composite is descriptive and does not define a new treated sample.

## 3. Treatment and donor eligibility

The 39-county treatment roster is frozen in `references/preprocessing/treated_NC_county_roster_2026-11-13.csv` and independently checked against `references/labor/frozen_treated_nc_39.csv`. Source construction compares it with the regenerated FEMA-positive North Carolina definition and can stop if the roster differs.

Donor exclusions are determined from county state FIPS:

- Impact-region exclusions: AL, FL, GA, IL, IN, KY, NC, OH, SC, TN, VA and WV.
- Buffer exclusions: AR, DC, IA, MD, MI, MS, MO, PA and WI.
- Other 2024 hurricane exclusions: LA and TX.

Eligible donors must also have all 24 January 2022–December 2023 observations for the four fitting outcomes, with positive employment, labor force and jobs and nonmissing unemployment rate. Later outcomes do not select donor eligibility. The pipeline requires at least 500 eligible donors and exactly 39 treated counties. Geographic separation is a design restriction, not a proof of no spillovers.

Implementation: `src/preprocessing/06_geography.do` and `src/labor/steps/01_frozen_sample.do`.

## 4. Transformations and the shared synthetic comparison

The four fitting outcomes are:

1. `100 × log(resident employment)`;
2. `100 × log(resident labor force)`;
3. `100 × log(workplace jobs)`;
4. unemployment rate in percentage points.

Each county series is initially centered on its own 2022 mean. Outcome scaling uses the treated aggregate's 2022 temporal standard deviation, with floors of 0.25 for the three log outcomes and 0.15 percentage points for unemployment rate. All four outcomes contribute to one shared set of county weights.

### Objective

For a fitting interval of T months and K = 4 outcomes, let `y[k,t]` be the treated mean and `X[k,t,j]` donor j's transformed outcome. Let `s[k]` denote its fixed scale and `u[j] = 1/J` the uniform donor vector. The final fit is equivalent to:

```text
minimize over w and outcome intercepts a:
  (1 / (2 K T)) × sum[k,t] ((y[k,t] - a[k] - sum[j] X[k,t,j] w[j]) / s[k])²
  + (lambda / 2) × sum[j] (w[j] - u[j])²

subject to w[j] >= 0 and sum[j] w[j] = 1.
```

The intercepts are implemented by demeaning each outcome block over the fitting interval. This fits pre-event movements while allowing a fixed treated/control level difference. Mata uses accelerated projected gradient with a simplex projection and an explicit optimality-gap check.

### Penalty selection and final fitting

The existing penalty-selection protocol is retained: fit 2022 and score 2023 over the grid `{0.03, 0.1, 0.3, 1, 3, 10}`. At this selection stage, the already 2022-centered outcome blocks are used without an additional interval-specific demeaning. A score within 1% of the minimum is considered near-best; ties choose the stronger regularization.

After selecting lambda, the definitive weights are fitted with outcome intercepts over **January 2022–August 2024**. A separate backdated diagnostic fit uses **January 2022–December 2023**, the same selected lambda and the intercept objective. It is not an alternative primary estimator.

The active method uses synthetic-DiD-style **unit weights**, outcome intercepts and uniform pre-event time averaging. It does not implement a separately estimated synthetic-DiD time-weight vector in the active pipeline. Calling it the full canonical synthetic-DiD estimator would therefore overstate the implementation.

Key outputs are `outputs/labor/tables/donor_weights.csv`, `inner_2023_cv.csv`, `training_2022_scales.csv`, `solver_audit.csv` and `monthly_paths.csv`.

## 5. Event study and period effects

The main estimation sample is January 2022–March 2026, excluding September 2024 (transition month) and October 2025 (missing LAUS). All six matched outcomes use the common observed sample. Treated county weights are 1/39; controls receive the shared fitted donor weights. Donors with numerical weights no greater than 1e-10 are excluded from the regression sample.

The regression is a county fixed-effect model with calendar-month effects and treated-by-event-month indicators, estimated using Stata `areg` with analytic weights and standard errors clustered by county. August 2024 is the omitted indicator for parameterization, but **reported coefficients are recentered to the mean of all 32 January 2022–August 2024 months**, using linear combinations and their joint covariance. The reported reference is therefore not August alone.

Six outcomes are reported:

| Outcome | Regression scale | Display scale |
|---|---|---|
| Resident employment | log count | 100 × log effect |
| Labor force | log count | 100 × log effect |
| Workplace jobs | log count | 100 × log effect |
| Unemployment rate | rate | percentage points |
| Jobs / resident employment | difference of log jobs and log employment | 100 × log ratio effect |
| Employment / labor force | difference of log employment and log labor force | 100 × log ratio effect |

The shared weights make the two ratio contrasts internally consistent. LAUS residents and QCEW workplace jobs nevertheless represent different populations; their levels are not interchangeable.

Three primary windows average the event coefficients:

- Immediate: October–December 2024 (3 months).
- 2025 observed window: January–September 2025 (9 months).
- Late matched window: November 2025–March 2026 (5 months).

The six outcomes times three windows form the **18-test Holm family**. Raw p-values and Holm-adjusted p-values are both retained. Displayed 95% intervals use county-clustered standard errors and a t critical value; they are individual intervals, not simultaneous Holm-adjusted confidence intervals.

For labor force, the saved tables additionally include one-sided lower confidence limits and whether those limits exclude losses of 1% or 2%, using `log(0.99)` and `log(0.98)` thresholds. These are distinct from declaring a zero effect.

Implementation: `src/labor/steps/04_joint_event_study.do` and `05_tables_figures.do`.

## 6. Diagnostics supplied by the current specification

The project reports full-pre-event fit, pre-event segment errors, maximum monthly gaps, linear drift contrasts, half-year-bin contrasts, strict joint monthly tests, counts of individually significant pre-event coefficients and the backdated 2024 prediction audit.

The definitive fit includes 2024, so its 2024 fit is in-sample. The backdated weights provide the implemented out-of-sample comparison for January–August 2024. Both sets of quantities are retained with distinct column names (`fit_*` versus `bd_*`). This document records their definition without changing or re-evaluating the accepted specification.

Implementation: `src/labor/steps/03_prefit_and_backdated_audit.do` and the event-study scripts.

## 7. Geographic placebo refits

The default is 199 draws with seed 20261114. Each draw samples 39 eligible donor counties as a pseudo-treated group and removes them from its own donor pool. The fit uses all four outcomes, the same full-pre-event intercept objective, and the primary selected lambda. Default convergence settings are 30000 iterations and tolerance 0.00005; strict mode rejects a stage with unconverged fits.

The observed and placebo shocks compare October–December 2024 with the pre-event mean gap. Output includes draws, convergence information, sampled county identities, and ranks for all converged draws and fixed relative pre-fit-error screens. Lower-tail and absolute-extremeness ranks use a +1 numerator/denominator correction.

These are geographic benchmarking ranks. Treatment was not randomized across these groups, so the ranks are not randomized causal p-values. No old employment-only placebo procedure is called.

Implementation and outputs: `src/labor/steps/07_joint_placebo_refits.do`; `outputs/labor/placebos/`.

## 8. Exploratory labor extensions and exposure curves

### Longer LAUS horizon

The extended LAUS analysis continues through September 2026 and reports additional windows, including November 2025–September 2026 and April–September 2026. It retains the shared primary comparison. Its p-values are exploratory and outside the 18-test main Holm family.

### Exposure associations

Within the treated counties, the code compares changes relative to the common synthetic path. Post-event county changes are rebased to each county's own pre-event mean gap. A January–August 2024 placebo window is measured against that county's 2022–2023 mean gap.

Exposure variables are total rainfall and HDSI. Outcomes include employment, labor force, jobs, unemployment rate and unemployed residents. Windows are the pre-event placebo, October–December 2024, January–September 2025 and November 2025–September 2026; workplace jobs in the last window end in March 2026. Incomplete county-window outcomes are marked missing rather than averaged over silently shorter periods.

The curves use local-linear Gaussian smoothing and 299 county bootstrap draws, with pointwise intervals and effective-sample-size support diagnostics. The code also reports linear exposure slopes per one sample standard deviation of the exposure, using HC3 uncertainty, when at least 15 observations are available.

These are descriptive within-treated-county exposure associations. The shared synthetic baseline changes levels, not the identifying source of cross-county exposure variation. Neither the curves nor the linear slopes establish a causal dose-response relationship. Pointwise bands and exploratory slope tests are not included in the main Holm family.

Implementation: `src/labor/steps/04b_extended_laus.do` and `06_kernel_rebased.do`.

## 9. Annual migration and population comparison

Census Vintage 2025 components produce a separate county-year panel for 2020–2025. Population growth is 100 times the annual population change divided by prior-year population. Net domestic migration per 1,000 residents uses prior-year population as its denominator.

The annual regression uses 2021–2025, county and year effects, and treated-by-year indicators relative to 2024. Standard errors are clustered by county. It compares the 39 treated counties with geographically eligible counties with available annual observations, **without applying the labor synthetic weights**. Group summaries are equal-county means.

The 2025 component spans July 2024–July 2025 and includes both pre- and post-Helene months. Net migration cannot distinguish increased departures from reduced arrivals or track the original residents individually. Annual estimates and diagnostics are exploratory and outside the main multiplicity family.

Implementation and outputs: `src/migration/01_census_v2025.do`; `outputs/migration/`.

## 10. Reproducibility, provenance and reporting

- Raw records, intermediate data, generated input panels, results and source code occupy separate directories.
- The full run checks county-month identifiers, expected panel size, frozen treatment and missingness before estimation.
- `outputs/audit/STAGE_STATUS.csv` gates result families in the report. File existence alone does not establish a successful stage.
- Full source builds and estimation-only reruns have different provenance; the latter also records `RERUN_STATUS.txt` and does not rebuild migration or sources.
- Historical logs retain historical paths. The reorganization map documents their current destinations instead of rewriting the historical record.
- Source-code checksums and raw-source hashes are distinct. Reporting additionally records the hashes of its table inputs in `outputs/reports/graphs/manifest.json`.
- The readable report and chart atlas transform presentation only. They do not recalculate weights, regressions or placebo fits.

## 11. Scope of inference and versioning

The design is observational and was revised after earlier outcomes had been inspected. It is not a preregistered design. Fixed rosters and reproducible code do not by themselves establish causal identification.

County clustering does not fully address correlated regional shocks or uncertainty from estimated donor weights. Geographic exclusions do not guarantee absence of spillovers. Source releases may revise observations; preserve the actual raw snapshots for replication. Migration, exposure curves and longer-horizon estimates are separate exploratory analyses.

The reorganization preserves the active specification and saved numerical results. This methodology document is descriptive: methodological changes belong in a separately documented analysis revision, not in a file-layout change.
