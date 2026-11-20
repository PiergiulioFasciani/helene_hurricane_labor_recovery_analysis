# Output index

All active generated products are under this directory. Source code is under `../src/`; source observations and working panels are under `../data/`. Start with the readable report and figure atlas rather than a derived dataset.

## Reader-facing results

| File | Contents |
|---|---|
| [reports/HELENE_COMPLETE_RESULTS.txt](reports/HELENE_COMPLETE_RESULTS.txt) | Interpretive summary, units, effect tables, diagnostics and scope |
| [reports/HELENE_GRAPH_ATLAS.pdf](reports/HELENE_GRAPH_ATLAS.pdf) | Printable atlas of the current saved results |
| [reports/graphs/index.html](reports/graphs/index.html) | Local gallery of atlas figures |
| [reports/graphs/manifest.json](reports/graphs/manifest.json) | Figure inventory, source-table hashes and unavailable families |

Regenerate these with `do "run_report_only.do"` from the project root. This does not rerun estimation.

## Audit and provenance (`audit/`)

- `STAGE_STATUS.csv`: full-run stage outcomes and return codes.
- `SOURCE_BUILD_AUDIT.txt`: source route and panel lineage recorded during construction.
- `SOURCE_SHA256.txt`: raw-file hashes recorded at the source-build hashing step.
- `CENSUS_STATUS.txt`: annual migration source availability and coverage.
- `RERUN_STATUS.txt`: return codes from an estimation-only rerun, when present.
- `REPORT_BUILD_OK.txt`: successful presentation-build marker.
- `REORGANIZATION_VALIDATION.json`: file-preservation and layout checks for the reorganization.

Earlier audit records are not rewritten to pretend they were produced at new paths. An old rerun audit is not the status of a later report-only refresh.

## Source construction (`preprocessing/`)

- `audits/`: source geography, county assignment, exposure, donor eligibility, missingness and coverage CSVs.
- `panels/Helene_NC_ES_preinference_panel.dta`: generated national panel; the estimation input copy is `data/input/Helene_NC_ES_preinference_panel.dta` relative to the project root.
- `logs/source_reconstruction.log`: source-build execution record.
- `graphs/`: reserved for source-stage graphics.

## Labor analysis (`labor/`)

| Location | Key products |
|---|---|
| `tables/donor_weights.csv` | Shared donor weights by county FIPS |
| `tables/solver_audit.csv` | Penalty, convergence and weight concentration |
| `tables/inner_2023_cv.csv` | Penalty-selection scores |
| `tables/monthly_paths.csv` | Treated, definitive and backdated comparison paths |
| `tables/pretrend_train_holdout_ALL.csv` | Current fit and backdated diagnostic fields |
| `tables/twfe_pretrend_tests.csv` | Saved event-study diagnostics |
| `tables/event_coefficients.csv` | Monthly main event estimates |
| `tables/period_effects_WITH_HOLM.csv` | Main 18-test effect family, raw and adjusted p-values |
| `tables/LAUS_extended_*` | Exploratory longer-horizon labor estimates |
| `tables/dose_event_coefficients.csv` | Exposure-interacted event study: monthly level and per-SD coefficients (HDSI, rainfall) and HDSI-tercile series |
| `tables/dose_period_effects.csv` | Window averages for the same terms, high-minus-low tercile contrast, Holm within each exploratory family |
| `tables/dose_pretrend_tests.csv` | Pre-event tests on the exposure terms |
| `../migration/migration_by_damage.csv` | Annual migration and population by HDSI tercile, two comparison schemes |
| `tables/dose_treated_exposure.csv` | Standardised exposures and HDSI tercile of each treated county |
| `placebos/joint_placebo_draws.csv` | Each pseudo-treated fit and its diagnostics |
| `placebos/joint_placebo_summary.csv` | Outcome ranks and fit-screen summaries |
| `placebos/joint_placebo_counties.csv` | County membership for every draw |
| `placebos/joint_placebos.log` | Placebo execution record |
| `data/` | Derived `.dta` panels and `.ster` estimation objects |
| `figures/` | Original stage-generated PNGs and `.gph` figures |
| `logs/` | Labor estimation logs |

The original stage graphics are distinct from the curated presentation gallery under `reports/graphs/`.

## Annual Census analysis (`migration/`)

- `census_county_year.csv` / `.dta`: annual population and migration panel.
- `census_group_comparison.csv`: unweighted treated/comparison group means.
- `census_annual_effects.csv` / `.dta`: annual regression contrasts and diagnostics.

## Interpretation and retention

CSV files are machine-readable numerical results. Stata datasets are working analysis products. Logs and status files document execution, not effects. Missing outputs may indicate failed or skipped stages; consult status before interpretation.

Historical outputs are stored in `../archive/`, outside the active result tree. Do not copy archived results into active folders to fill missing stages. Keep raw data and reference rosters when cleaning regenerable outputs.
