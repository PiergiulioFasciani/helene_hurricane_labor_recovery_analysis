# Hurricane Helene research pipeline

A source-to-results Stata project for studying county labor-market outcomes after Hurricane Helene, with a shared four-outcome synthetic comparison, geographic placebo refits, exploratory exposure-response curves, and annual Census migration comparisons.

The active analysis is the existing **14 November 2026 specification**. This repository reorganization changes file locations and documentation, not the estimator, treatment roster, source observations, or saved numerical results. See [METHODOLOGY.md](METHODOLOGY.md) for the implemented design.

## Start here

| What you want | Open or run |
|---|---|
| Read the results | [Readable results](outputs/reports/HELENE_COMPLETE_RESULTS.txt) |
| View all current figures | [PDF atlas](outputs/reports/HELENE_GRAPH_ATLAS.pdf) or [browser gallery](outputs/reports/graphs/index.html) |
| Understand the design | [Methodology](METHODOLOGY.md) |
| Find a specific output | [Output index](outputs/README.md) |
| Inspect execution and source provenance | [Audit folder](outputs/audit/) |
| Rebuild everything from raw sources | `do "run_all.do"` from the project root |
| Re-estimate from the existing generated panel | `do "rerun_estimation_only.do"` from the root |
| Refresh the report and figures only | `do "run_report_only.do"` from the root |

The PDF, gallery and result tables are generated files. They exist in this working copy but are not bundled automatically in a clean source-code checkout.

## Repository layout

```text
Helene_Definitive_Joint_Research/
├── README.md                       # Setup, execution and navigation
├── METHODOLOGY.md                   # Implemented statistical design
├── run_all.do                      # Full source-to-results entry point
├── rerun_estimation_only.do        # Re-estimation using existing panel
├── run_report_only.do              # Presentation refresh only
├── config.example.do               # Public configuration template
├── config.local.do                 # Private overrides; never commit
├── src/
│   ├── build_panel_from_sources.do  # Mandatory builder for a full run
│   ├── preprocessing/              # NOAA, Census, FEMA, labor, geography, panel
│   ├── labor/
│   │   ├── run_*.do                # Labor-stage entry points
│   │   ├── config.example.do       # Internal solver defaults
│   │   └── steps/                  # Fitting, diagnostics, inference, curves, placebos
│   ├── migration/                  # Annual Census construction and comparisons
│   └── reporting/                  # Readable TXT, PDF atlas, PNG/HTML gallery
├── references/
│   ├── preprocessing/              # Frozen treatment roster used by source build
│   └── labor/                      # Frozen treatment roster used by estimator
├── data/
│   ├── raw/                        # Cached external source snapshots
│   ├── interim/                    # Regenerable intermediate datasets
│   └── input/                      # Generated national panel used by estimation
├── outputs/
│   ├── README.md                   # Detailed output index
│   ├── audit/                      # Status, source lineage and source hashes
│   ├── preprocessing/              # Source audits, panels and build logs
│   ├── labor/                      # Labor tables, datasets, figures, logs, placebos
│   ├── migration/                  # Annual county panel and migration tables
│   └── reports/                    # Final readable results and curated graph atlas
├── docs/                           # Source setup, protocol and relocation map
├── scripts/                        # Raw-source checksum utility
├── tests/                          # Structure, source-contract and reporting checks
├── archive/                        # Historical snapshots; never executed by active runners
├── SHA256SUMS.txt                   # Current distributable source-file checksums
└── LICENSE                         # Code license
```

Source code never writes analysis output into `src/`. Frozen reference rosters are separate from generated sample-assignment tables. Historical run logs may retain the paths used when they were created; [REORGANIZATION_MAP.json](docs/REORGANIZATION_MAP.json) maps old locations to current ones.

## Requirements

- **Stata:** source files declare version 16.0; the working project has been executed with Stata 19.5. Event-study stages use `areg ... vce(cluster county_id)` with analytic weights; the synthetic-control solver is Mata (FISTA with simplex projection).
- **Python 3:** NumPy, pandas and Matplotlib for reporting. Install into the interpreter you want Stata to use:

  ```sh
  python3 -m pip install -r src/reporting/requirements.txt
  ```

- **Network access:** needed when required raw files are not already cached.
- **Shell utilities on macOS:** `/usr/bin/curl`, `/bin/sh` and `shasum` support fallback downloads and raw-source hashing.
- **Storage:** keep raw snapshots, intermediate panels and outputs locally; the raw downloads are the largest part of the project.

No external Stata estimation package is required by the active fitting code; the optimizer is implemented in Mata.

## Private configuration

For a new checkout, copy `config.example.do` to `config.local.do`. Do not overwrite an existing private configuration.

| Setting | Purpose / default |
|---|---|
| `H_CENSUS_KEY` | Census API key for source acquisition |
| `H_LAUS_SOURCE_MODE` | `mirror` by default; `official` requests BLS flat files |
| `H_JOINT_PLACEBO_REPS` | 199 pseudo-treated group refits |
| `H_JOINT_PLACEBO_SEED` | 20261114 |
| `H_JOINT_PLACEBO_ITER` | 30000 iterations maximum per placebo fit |
| `H_JOINT_PLACEBO_TOL` | 0.00005 convergence tolerance |
| `H_JOINT_PLACEBO_STRICT` | 1; unsuccessful convergence fails the stage |
| `H_RUN_CENSUS` | 1; include annual migration module |
| `H_CENSUS_FILE` | Optional local path to the Census annual source CSV |
| `H_PYTHON` | Optional absolute path to Python with reporting dependencies installed |

The primary optimizer's existing defaults are in `src/labor/config.example.do`: maximum 35000 iterations and KKT tolerance 0.00002. Changing a design or solver setting is a substantive configuration change, not part of this reorganization.

`config.local.do` is ignored by Git and excluded from the source checksum manifest and shareable source package. Keep API keys out of logs, reports and published files.

## Running the project

### Full source-to-results run

In Stata:

```stata
cd "~/Downloads/Helene_Definitive_Joint_Research"
do "run_all.do"
```

The runner:

1. Rebuilds the national panel from raw source records inside this repository.
2. Validates county-month uniqueness, expected size, treatment roster and required missingness.
3. Fits the active shared four-outcome model and produces labor estimates and exploratory curves.
4. Refits geographic placebo groups and processes annual Census migration when enabled.
5. Regenerates readable results and figures from the stages recorded as passed.

Existing raw files are cached inputs. The full runner does not import a prior combined analysis panel from a sibling project. It replaces its generated input panel after successful construction. If source construction or panel validation fails, downstream estimation is blocked.

### Estimation-only run

```stata
do "rerun_estimation_only.do"
```

Requires `data/input/Helene_NC_ES_preinference_panel.dta`. Re-estimates labor and placebos and then requests reporting. It does not download sources, rebuild the source panel or re-estimate migration. Its return-code audit is `outputs/audit/RERUN_STATUS.txt`; retain the full-run `STAGE_STATUS.csv` as the earlier source-stage record and inspect both when assessing provenance.

### Report-only run

```stata
do "run_report_only.do"
```

Or from the project root in a terminal:

```sh
python3 src/reporting/make_report.py
```

This reads saved CSVs, writes the readable text, PDF and gallery, and never fits an estimator. The report checks `outputs/audit/STAGE_STATUS.csv`; an unpassed stage does not supply numerical results merely because old files exist.

If Stata selects a Python without the required packages, set `H_PYTHON` to the absolute path of the correctly configured interpreter. The reporting wrapper also tries standard macOS Python locations. `outputs/audit/REPORT_BUILD_OK.txt` records successful report generation.

### Individual labor stages

For advanced use, the module entry points are in `src/labor/`: `run_multicontrol.do` and `run_joint_placebos.do`. They expect that directory as the working directory and the required earlier outputs to exist. Prefer the root entry points for a coherent pipeline run; manual module runs do not independently refresh the master stage-status table.

## Inputs and provenance

See [RAW_DATA_SETUP.md](docs/RAW_DATA_SETUP.md) for exact source URLs and cached filenames. The active inputs are NOAA county precipitation, Census PEP and ACS, FEMA registration-level verified losses, LAUS labor records, official QCEW county files, and Census annual migration components.

The default LAUS route is **economicdata.dev, an independent distributor of BLS series**, not a direct BLS bulk download. The route is explicitly documented in the build audit. Preserve source snapshots and hashes because vintages can change. The alternative official route may require manual caching when automated BLS access is unavailable.

`SHA256SUMS.txt` describes current source files. `outputs/audit/SOURCE_SHA256.txt` describes raw snapshots hashed during the source build. These serve different purposes. The source hash step occurs before the optional migration stage, so it should not be assumed to cover a newly downloaded annual migration file unless that file was already present at hashing time.

## Output conventions

- `outputs/labor/tables/`: analysis CSVs, suitable for independent inspection or plotting.
- `outputs/labor/data/`: derived Stata datasets and saved estimation objects.
- `outputs/labor/figures/`: original stage-generated PNG and Stata graph files.
- `outputs/labor/placebos/`: draws, summaries, sampled counties and placebo log.
- `outputs/reports/graphs/`: curated presentation PNGs, local HTML gallery and input-hash manifest.
- `outputs/migration/`: annual county data and comparison tables.
- `outputs/audit/`: execution and provenance records; these are not effect estimates.

The report's `graphs/manifest.json` lists inputs and their hashes, generated pages and unavailable families. The input CSVs are retained rather than pasted wholesale into the readable TXT.

## Validation and troubleshooting

From the root:

```sh
python3 tests/check_repository.py
python3 tests/test_source_first.py
python3 tests/test_reporting.py
```

The first two inspect source structure and the source-first contract. The reporting test verifies that stale estimates cannot appear when their stage failed. These checks do not replace an end-to-end source and estimation run.

| Symptom | Inspect |
|---|---|
| Missing source file / `r(601)` | `outputs/preprocessing/logs/source_reconstruction.log`, cached raw files, source audit |
| Failed panel assertion | Source geography, coverage and treatment audits in `outputs/preprocessing/audits/` |
| Labor-stage error | `outputs/labor/logs/` and solver audit tables |
| Placebo-stage error | `outputs/labor/placebos/joint_placebos.log` and convergence columns in the draw table |
| Reporting failure | Python dependencies, `H_PYTHON`, and `outputs/audit/REPORT_BUILD_OK.txt` |
| Expected output absent | Stage status first; a skipped or failed stage must not be replaced with historical output |

A `PASSED` stage records execution success. Interpretation follows the implemented methodology and reported uncertainty; the reorganization does not reinterpret the analysis.

## Archives, sharing and maintenance

`archive/` holds historical snapshots and discontinued working material. Active runners do not reference or execute them. Archived code retains historical paths and is not a current entry point.

Raw data, private configuration, generated datasets and large outputs are ignored by Git. Share source files and reference rosters separately from source snapshots; provide results as a distinct bundle when needed. Never infer that a source-only ZIP is a complete data snapshot.

The code is distributed under [LICENSE](LICENSE). Source datasets retain their providers' terms and attribution requirements.

## Uploading to GitHub

The repository includes a `.gitignore`, LF line-ending rules, source checksums, and a GitHub Actions workflow. Empty working directories are retained by `.gitkeep` files. Raw data, generated outputs, historical snapshots, local environments and private configuration are excluded.

For browser upload, extract `Helene_GitHub_Ready.zip` and upload the **contents** of its project folder. Include the hidden `.github/`, `.gitignore`, `.gitattributes` and `.gitkeep` files; do not upload the ZIP itself as the repository's only file. The local `.git/` directory is not included in the ZIP.

For Git, run from the project root. Initialize first when starting from the extracted ZIP; initialization is already done in the working repository. Then review `git status` before staging:

```sh
git init -b main
git status
git add .
git commit -m "Prepare Helene research pipeline"
git remote add origin YOUR_GITHUB_REPOSITORY_URL
git push -u origin main
```

The project has a local `main` branch initialized; no remote, commit or upload is created automatically. Never use `git add -f` to override the exclusions for private configuration or data.

GitHub Actions runs structure, source-contract, reporting-isolation and packaging checks on pushes and pull requests. It does **not** download source data or run Stata estimation; a licensed Stata runtime and local source snapshots are needed for that. The source checksum manifest includes every upload file except itself. After source changes, run `python3 scripts/update_source_manifest.py`, then `python3 tests/test_github_ready.py`.
