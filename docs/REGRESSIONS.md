# How the regressions actually run

Every estimate in this project comes from one of **five** regressions. This
file is the map; `outputs/audit/REGRESSIONS.txt` is the receipt, written by the
estimation stages themselves on each run, so it reports what Stata executed
rather than what a comment claims. If the two ever disagree, the receipt is
right.

Nothing else in the pipeline estimates anything. Stages 01-03 build the sample
and the donor weights, stage 05 formats, stage 07 refits the weights on placebo
groups, and the reporting layer only reads CSVs.

## The one specification, in words

All three labor regressions are the same model with the treated indicators
split differently:

```
areg  y   <treated-by-event-month indicators>   i.mdate
      if <sample>   [aw=<weight>],
      absorb(county_id)  vce(cluster county_id)
```

| Element | What it is |
|---|---|
| `y` | one of six outcomes, each already in display units (100 x log, or percentage points) |
| event indicators | treated county x event-month dummies, event time `k = mdate - ym(2024,9)` |
| `i.mdate` | calendar-month effects, common to treated and donors |
| `absorb(county_id)` | county fixed effects |
| `[aw=]` | treated county `1/39`; donor county its fitted synthetic weight `mc_w_joint` |
| `if` | treated counties, plus donors whose weight exceeds `1e-10` |
| `vce` | clustered on county: 154 clusters, `df_r = 153` |

**Three event months carry no indicator**, and this trips people up:

- `k = -1` (**August 2024**) is the omitted category that parameterises the model.
- `k = 0` (**September 2024**) is the storm month itself, dropped from the sample.
- `k = 13` (**October 2025**) has no LAUS data.

**August 2024 is not the reference.** The omitted month only parameterises the
regression. Every reported coefficient is then recentred to the **mean of all
32 pre-event months** (January 2022 - August 2024) with `lincom`:

```
beta_k(reported) = beta_k - (1/32) * sum over the 32 pre-event months
```

so a reported coefficient answers "relative to the typical pre-storm gap", not
"relative to August 2024". The same recentring is applied inside each damage
tercile. This is why the code builds a `premean` macro and never reports a raw
`_b[]` for an event month.

The **linear pre-trend test** is one `lincom` over the pre-event coefficients
weighted by `k + 16.5` (the centre of the 32-month window), so it is a test of
drift, not of level.

## The five regressions

| # | What | File | Estimates | Writes |
|---|---|---|---|---|
| 1 | **Main event study** | `src/labor/steps/04_joint_event_study.do` | 6 outcomes x 1 model | `event_coefficients.csv`, `period_effects_raw.csv`, `period_effects_WITH_HOLM.csv` |
| 2 | **Damage terciles** | `src/labor/steps/04c_dose_event_study.do` | 6 outcomes x 2 families | `dose_event_coefficients.csv`, `dose_period_effects.csv` |
| 3 | **Extended LAUS** | `src/labor/steps/04b_extended_laus.do` | 5 outcomes, exploratory | `LAUS_extended_*_EXPLORATORY.csv` |
| 4 | **Migration, pooled** | `src/migration/02_migration_dose.do` | 2 outcomes x 2 schemes | `migration_by_damage.csv` (`group = all`) |
| 5 | **Migration, by tercile** | `src/migration/02_migration_dose.do` | 2 outcomes x 2 schemes | `migration_by_damage.csv` (`group = low/mid/high/hi-lo`) |

### 1. Main event study

The headline result. Six outcomes share one set of donor weights, which is what
makes the two ratio outcomes (`jobs/emp`, `emp/lf`) internally consistent.
Period effects average the event coefficients over three windows: immediate
(Oct-Dec 2024), 2025 (Jan-Sep 2025) and late matched (Nov 2025-Mar 2026). The
six outcomes x three windows are the **18-test Holm family**.

### 2. Damage terciles

Identical sample, weights, fixed effects and clustering to #1. The only change
is that each treated-by-event-month indicator is split three ways by exposure
tercile, 13 treated counties each, so each tercile gets its own event path
against the same donors. Because all three terciles are in **one** regression,
the high-minus-low contrast is taken with the joint covariance rather than by
differencing two separate models.

`hdsi_tercile` (composite damage) is primary; `rain_tercile` (six-day rainfall
alone, no FEMA claim behaviour in it) is the robustness footnote.

Each tercile effect is an ATT for that exposure group against never-treated
donors, so it needs only the parallel-trends assumption of #1 applied within
the group. Exposure is **not** randomly assigned across the 39 counties.

### 3. Extended LAUS

Same design pushed to September 2026 on LAUS alone, past where QCEW can match
it. Exploratory, outside the Holm family.

### 4 and 5. Migration

A different panel: annual, July-to-July, clustered on `demo_id`. 2024 is the
omitted last unexposed year and the reported effect is the 2025 coefficient;
2021-2023 are pre-tests. Two weighting schemes are run, `unwt` and `syn`;
`syn` is primary because the unweighted scheme fails its pre-test for the two
less-damaged terciles.

Migration is an **outcome**, never a control. It is realised after the storm,
so conditioning on it would be a bad control and would remove the channel under
study.

## Reading the receipt

`outputs/audit/REGRESSIONS.txt` is regenerated on every pipeline run. Each
block reports the command, the regressors, the sample, the weights, the
clustering, the reference, and the `N`/clusters/`df_r` that Stata actually
returned. To see what is being estimated without reading any Stata:

```
open outputs/audit/REGRESSIONS.txt
```
