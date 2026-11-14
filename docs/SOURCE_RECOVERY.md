# Automatic LAUS source recovery and provenance (2026-11-14)

## Why the preceding source-first run stopped

In the supplied `source_reconstruction.log`, NOAA, Census 2024 PEP, ACS 2023,
FEMA 4827--4832 and the exposure merge completed. The next stage attempted
`https://download.bls.gov/pub/time.series/la/la.area`; BLS returned HTTP 403.
The repository stopped with Stata `r(601)` **before** constructing the national
labor panel. This was a code regression introduced when an official-only
restriction replaced a previously working automatic fallback.

## Exact source route used in the original successful scripts

The original `Helene_LAUS_QCEW_403_FALLBACK.do` and the patched pre-inference
`05_labor.do` obtained LAUS county-month observations from:

`https://api.economicdata.dev/v1/laus/download/laus.csv`

This is an **independent data distributor** that republishes BLS LAUS records.
It is not the official BLS bulk-download server. The original script filtered
`series == "nsa"` and `location` beginning `us-county-`, parsed 5-digit FIPS,
and imported resident employment, unemployment, labor force and the rate.
QCEW came directly from official BLS quarterly industry=10 API CSVs.

## What the corrected repository now does

- Default `H_LAUS_SOURCE_MODE "mirror"` restores the original automatic source.
- Downloads the distributor's full CSV to
  `data/raw/LAUS/LAUS_mirror_all_counties.csv` (and attempts to save metadata).
- Checks that the downloaded file has the expected source-column header and
  that imported fields and county-month identifiers are valid.
- Keeps original CSVs in `data/raw/`, checks SHA-256 in
  `outputs/audit/SOURCE_SHA256.txt`, and records the LAUS source URL,
  provider and mode in `SOURCE_BUILD_AUDIT.txt`.
- Rebuilds the 2022m1--2026m9 calendar and the merged national panel **inside
  the repository**, never from a previous `.dta`.
- Preserves October 2025 as a calendar month with missing LAUS outcomes;
  QCEW is not imputed or altered in that month.
- Checks the 3,221-by-57 county-month size, all 39 frozen treated NC counties,
  and existing source-geography and donor rules. It does not alter the sample
  to fit results.
- Optional `H_LAUS_SOURCE_MODE "official"` uses actual BLS flat files, but
  can require manual downloads due to HTTP 403. This is **not** the default.

## Important replication limits

The independent distributor updates its data. Rebuilding today may not give
bitwise identical labor estimates to a CSV obtained on a prior date. Preserve
raw CSV snapshots and their checksums and inspect prior-vs-new reference
comparisons. Do not claim the direct BLS bulk-data endpoint was used when the
mirror was used. The supplied log additionally flagged all 1,021 NOAA/Census
county **name** comparisons as disagreements despite FIPS merges; this
remains an open validation issue until `noaa_census_name_disagreements.csv`
is inspected. A successful build alone cannot certify this geography check.

## Repair the existing working folder (keeps cached downloaded sources)

Overwrite `src/build_panel_from_sources.do`, `run_all.do`,
`src/preprocessing/05_labor.do`, and the supporting README/config/test
files with the attached patch. **Do not delete the existing `data/raw/`**.
Then in Stata:

```
cd "~/Downloads/Helene_Unified_Research_GitHub_SOURCE_FIRST"
do "run_all.do"
```

The stage must fetch/validate LAUS source records before downstream inference.
If a new source revision changes treatment, county coverage or model fit, the
builder must stop or report the difference, not silently force agreement.

### Provenance

- Third-party distributor: https://economicdata.dev/datasets/laus/
- Distributor bulk CSV: https://api.economicdata.dev/v1/laus/download/laus.csv
- Original government source: https://www.bls.gov/lau/
- Official QCEW county CSV API: https://data.bls.gov/cew/data/api/
