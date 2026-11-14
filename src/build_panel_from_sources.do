/* SOURCE-FIRST PRE-INFERENCE BUILD. Run by root run_all.do ONLY.
   Inputs: raw public NOAA/Census/FEMA/BLS sources, stored under data/raw/.
   Outputs: in-repo generated county-month panel and intermediates.
   No preconstructed national labor or severity data may be loaded. */
version 16.0
if "${H_ROOT}"=="" {
    di as error "H_ROOT is unset. Run the TOP-LEVEL run_all.do."
    exit 198
}

* Do not read src/preprocessing/config.local.do: private overrides could
* quietly reactivate historical pre-built datasets.
global ROOT "${H_ROOT}/src/preprocessing"
global RAW "${H_ROOT}/data/raw"
global INT "${H_ROOT}/data/interim"
global OUT "${H_ROOT}/outputs/preprocessing"
global REF "${H_ROOT}/references/preprocessing"
global LOG "${OUT}/logs"

* Critical safeguards: source records ONLY, no derived-data cache. The
* third-party BLS-data mirror is explicit, not mislabeled as an official file.
global FEMA_SOURCE "api"
global LABOR_SOURCE "rebuild"
global LAUS_SOURCE_MODE "${H_LAUS_SOURCE_MODE}"
global STOP_IF_TREATMENT_CHANGED 1
global CENSUS_KEY "${H_CENSUS_KEY}"
global QCEW_2026_MAXQ 1

global EXPECTED_END_MDATE = ym(2026,9)

foreach d in "data" "data/raw" "data/interim" "data/input" "outputs" "outputs/audit" ///
    "outputs/preprocessing" "outputs/preprocessing/panels" ///
    "outputs/preprocessing/audits" ///
    "outputs/preprocessing/graphs" "outputs/preprocessing/logs" ///
    "outputs" "outputs/reports" {
    capture mkdir "${H_ROOT}/`d'"
}
capture log close _all
log using "${LOG}/source_reconstruction.log", text replace name(srcbuild)
di as result "SOURCE-FIRST BUILD, raw=${RAW}"
di as text "BLS LAUS CSV route: ${LAUS_SOURCE_MODE}; FEMA v2 raw API rows."

do "${ROOT}/01_noaa.do"
do "${ROOT}/02_census.do"
do "${ROOT}/03_fema.do"
do "${ROOT}/04_exposure_index.do"
do "${ROOT}/05_labor.do"
do "${ROOT}/06_geography.do"
do "${ROOT}/07_final_es_panel.do"

* The master input path is ONLY populated from the in-repository source build.
local generated "${OUT}/panels/Helene_NC_ES_preinference_panel.dta"
capture confirm file "`generated'"
if _rc {
    di as error "Source build finished without producing `generated'."
    exit 601
}
use "`generated'", clear
isid fips mdate
assert _N==3221*57
quietly count if mdate==ym(2024,8) & treated_nc==1
assert r(N)==39
quietly count if mdate==ym(2025,10) & !missing(laus_employed)
assert r(N)==0
copy "`generated'" "${H_ROOT}/data/input/Helene_NC_ES_preinference_panel.dta", replace

* Stata-only source contract: URL/definition strings not obtained from a cache.
tempname audit
file open `audit' using "${H_ROOT}/outputs/audit/SOURCE_BUILD_AUDIT.txt", write text replace
file write `audit' "SOURCE-FIRST BUILD; run date: `c(current_date)'" _n
file write `audit' "Root: ${H_ROOT}" _n
file write `audit' "Panel generated: `generated'" _n
file write `audit' "Panel input convenience copy: ${H_ROOT}/data/input/Helene_NC_ES_preinference_panel.dta" _n
file write `audit' "Census key is NOT recorded (private)." _n
file write `audit' "NOAA: Sep25-Sep30 2024 NClimGrid-Daily county totals; raw CSV under data/raw/" _n
file write `audit' "Census: official 2024 county PEP and 2023 ACS5 housing data" _n
file write `audit' "FEMA: official v2 IHP valid registrations; disasters 4827-4832" _n
file write `audit' "BLS LAUS source mode: ${LAUS_SOURCE_MODE}" _n
file write `audit' "BLS LAUS source: ${LAUS_SOURCE_DESC}" _n
file write `audit' "BLS LAUS raw file: ${LAUS_SOURCE_FILE}" _n
file write `audit' "BLS LAUS URL: ${LAUS_SOURCE_URL}" _n
file write `audit' "LAUS mirror is a THIRD PARTY transformation of BLS series; not direct BLS raw." _n
file write `audit' "BLS QCEW: official quarterly industry 10 county CSV; 2022Q1-2026Q1" _n
file write `audit' "Treated FIPS: 39 frozen counties, verified against regenerated FEMA results; only joint geography and 2022-23 four-outcome coverage" _n
file write `audit' "No prebuilt labor, severity or pre-inference .dta accepted as source; NO employment-only SCM stage." _n
file close `audit'
log close srcbuild
