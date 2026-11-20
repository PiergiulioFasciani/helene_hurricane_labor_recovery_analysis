# Methodology

This document describes the **active implementation**. Historical versions are
stored under `archive/` and are not executed by the current pipeline.

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

For a fitting interval of T months and K = 4 outcomes, let `y[k,t]` be the
treated mean, `X[k,t,j]` donor j's transformed outcome, `s[k]` its fixed scale
and `u[j] = 1/J` the uniform donor vector. The fit is:

```text
minimize over w and outcome intercepts a:
  (1 / (2 K T)) x sum[k,t] ((y[k,t] - a[k] - sum[j] X[k,t,j] w[j]) / s[k])^2
  + (lambda / 2) x sum[j] (w[j] - u[j])^2

subject to w[j] >= 0 and sum[j] w[j] = 1.
```

The intercepts are implemented by demeaning each outcome block over the fitting
interval, which fits pre-event movements while allowing a fixed treated/control
level difference. The ridge term shrinks toward uniform weights rather than
toward zero, so it spreads weight across donors instead of discarding them.
Mata solves it with accelerated projected gradient, a simplex projection and an
explicit optimality-gap check.

Lambda is selected by fitting 2022 and scoring 2023 over the grid
`{0.03, 0.1, 0.3, 1, 3, 10}`, with ties resolved toward stronger
regularization. The definitive weights are then fitted with intercepts over
**January 2022-August 2024**. A backdated diagnostic fit uses **January
2022-December 2023** with the same lambda and objective; it is a check, not an
alternative primary estimator.

The method uses synthetic-DiD-style **unit weights**, outcome intercepts and
uniform pre-event time averaging. It does not estimate a separate
synthetic-DiD time-weight vector, so calling it the full canonical synthetic
DiD estimator would overstate the implementation.

Derivations, the duality-gap stopping rule and the full selection table are in
`docs/methods_note/helene_methods.tex`.

Key outputs: `outputs/labor/tables/donor_weights.csv`, `inner_2023_cv.csv`,
`training_2022_scales.csv`, `solver_audit.csv` and `monthly_paths.csv`.

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

Implementation: `src/labor/steps/04_joint_event_study.do` and `05_tables_figures.do`.

## 6. Pre-event diagnostics

Two tests are reported per outcome, and they answer different questions.

- **Drift.** A linear contrast across the 32 pre-event months. This is the test
  that bears on the design: a comparison that was already trending apart will
  attribute that trend to the storm.
- **Strict joint zero.** A joint test that all 31 free pre-event coefficients
  are exactly zero. This rejects for every outcome, at month-to-month
  deviations of roughly 0.05-0.15 units against estimated effects an order of
  magnitude larger. It is reported for completeness, not as a screen. The count
  of individually significant pre-event coefficients and the largest absolute
  pre-event coefficient are reported alongside it, to give the scale of the
  noise it is rejecting.

The definitive weights are fitted over the whole pre-event window, so flatness
there is **in-sample and is not evidence**. The one out-of-sample check in the
design is the backdated audit: weights fitted on 2022-2023 only are used to
predict January-August 2024, which they never saw, and the ratio of forecast to
training RMSE is reported.

Implementation: `src/labor/steps/03_prefit_and_backdated_audit.do` and the
event-study scripts.

## 7. Geographic placebo refits

The default is 199 draws with seed 20261114. Each draw samples 39 eligible donor counties as a pseudo-treated group and removes them from its own donor pool. The fit uses all four outcomes, the same full-pre-event intercept objective, and the primary selected lambda. Default convergence settings are 30000 iterations and tolerance 0.00005; strict mode rejects a stage with unconverged fits.

The observed and placebo shocks compare October-December 2024 with the
pre-event mean gap. Output includes the draws, a convergence flag, the sampled
county identities and ranks over all converged draws. Lower-tail and
absolute-extremeness ranks use a +1 numerator/denominator correction.

These are geographic benchmarking ranks. Treatment was not randomized across these groups, so the ranks are not randomized causal p-values. No old employment-only placebo procedure is called.

Implementation and outputs: `src/labor/steps/07_joint_placebo_refits.do`; `outputs/labor/placebos/`.

## 8. Exploratory longer-horizon labor extension

The extended LAUS analysis continues through September 2026 and reports additional windows, including November 2025–September 2026 and April–September 2026. It retains the shared primary comparison. Its p-values are exploratory and outside the 18-test main Holm family.

Implementation: `src/labor/steps/04b_extended_laus.do`.

## 8b. Damage-tercile event study

The damage gradient is estimated **inside** the event study, with the same
sample, analytic weights, county and month effects, county clustering and
pre-event-mean reference as Section 5. Nothing in Section 5 is replaced. The 39
treated counties are split into terciles of exposure (13 each) and each tercile
gets its own event-study series against the same synthetic comparison; the
high-minus-low window contrast is reported with its joint covariance.

- **Primary:** `hdsi_tercile`, terciles of the composite HDSI.
- **Robustness:** `rain_tercile`, terciles of six-day rainfall alone, i.e.
  meteorology with no FEMA claim behaviour in it. Reported as a footnote
  comparison, not as a second set of headline results.

Each tercile effect is an ATT for that exposure group against never-treated
donors, so it requires only the parallel-trends assumption of Section 5 applied
within the group. **Exposure is not randomly assigned** across the treated
counties: damage correlates with terrain, income and industry mix, so the
high-minus-low contrast compares differently exposed counties and is not a
randomised dose.

Multiplicity: Holm is applied within this family, which sits outside the main
18-test family.

Implementation and outputs: `src/labor/steps/04c_dose_event_study.do`;
`outputs/labor/tables/dose_*.csv`; atlas pages `03b`-`03c`.

## 9. Annual migration and population comparison

Census Vintage 2025 components produce a separate county-year panel for 2020–2025. Population growth is 100 times the annual population change divided by prior-year population. Net domestic migration per 1,000 residents uses prior-year population as its denominator.

The annual regression uses 2021–2025, county and year effects, and treated-by-year indicators relative to 2024. Standard errors are clustered by county. It compares the 39 treated counties with geographically eligible counties with available annual observations, **without applying the labor synthetic weights**. Group summaries are equal-county means.

The 2025 component spans July 2024-July 2025 and includes both pre- and
post-Helene months. Net migration cannot distinguish increased departures from
reduced arrivals, and does not track original residents individually. Annual
estimates are exploratory and outside the main multiplicity family.

**A negative contrast is slower net in-migration, not population loss.** Net
migration stays positive for the treated group throughout: the high tercile
falls from about +11.1 to +4.8 per 1,000. Two of the 39 counties show true net
outflow and six show falling population. The damaged counties attracted fewer
new residents than they otherwise would have; they did not depopulate, and the
paper must not say they did.

The same annual comparison is also run **by damage tercile**, using the HDSI groups of Section 8b, so that migration enters the design as an outcome of the dose-response framework rather than as a control. Migration is realized after treatment: conditioning on it would be a bad control and would remove the channel under study. Two comparison schemes are reported, unweighted eligible counties and donors carrying the labor synthetic-control weights; the pre-tests on the 2021–2023 treated-by-year coefficients discriminate between them.

Implementation and outputs: `src/migration/01_census_v2025.do`, `src/migration/02_migration_dose.do`; `outputs/migration/`.

## 10. Reproducibility and provenance

- Raw records, intermediate data, generated panels, results and source code
  occupy separate directories.
- The full run checks county-month identifiers, panel size, the frozen
  treatment roster and missingness before estimation.
- `outputs/audit/STAGE_STATUS.csv` gates result families in the report. A file
  existing does not establish that its stage succeeded.
- Full source builds and estimation-only reruns have different provenance; the
  latter records `RERUN_STATUS.txt` and does not rebuild migration or sources.
- Source-code checksums and raw-source hashes are distinct. Reporting records
  the hashes of its table inputs in `outputs/reports/graphs/manifest.json`.
- The readable report and chart atlas transform presentation only. They do not
  recalculate weights, regressions or placebo fits.
- Historical logs retain historical paths and are not rewritten.

## 11. Scope of inference

The design is observational and was revised after earlier outcomes had been
inspected. It is not preregistered. Fixed rosters and reproducible code do not
by themselves establish causal identification.

County clustering does not fully address correlated regional shocks or the
uncertainty from estimated donor weights, and the geographic exclusions do not
guarantee the absence of spillovers. Source releases may revise observations,
so the raw snapshots should be preserved for replication. Migration,
longer-horizon and damage-tercile estimates are exploratory relative to the
main 18-test family.

This document is descriptive. Methodological changes belong in a documented
analysis revision, recorded below with a date.

## Appendix. Removed from earlier versions

- **Gaussian exposure curves** (stage `06_kernel_rebased.do`, 18 Nov 2026). The
  damage gradient is estimated inside the event study instead (Section 8b), so
  it needs no separate assumptions. Preserved under
  `archive/kernel_dropped_2026-11-18/`.
- **Continuous per-SD dose families** (18 Nov 2026). A per-SD slope is a causal
  response only under a stronger parallel-trends assumption than the binned
  tercile contrasts require.
- **Placebo fit screens** (18 Nov 2026). Ranks are now computed over all
  converged draws.
- **Duplicate pre-event diagnostics** (19 Nov 2026): in-sample fit columns,
  half-year bin contrasts, the last-eight-month contrast, pre-event RMSE and
  the maximum pre-event t-statistic. Each duplicated one of the two tests in
  Section 6 or described the in-sample fit.
- **One-sided labor-force limits** and the 1%/2% rule-out flags (19 Nov 2026).
  Not reported in any table or figure.
