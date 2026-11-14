/* HELENE DEFINITIVE JOINT-ONLY pipeline: source-built panel -> four-outcome
   synthetic control -> diagnostics, matched TWFE, kernel, new joint placebos
   -> annual migration -> one TXT report. No historical estimators. */
version 16.0
clear all
set more off
set varabbrev off
set graphics off
capture log close _all
local hroot "`c(pwd)'"
global H_ROOT "`hroot'"
foreach dir in "data" "data/input" "data/raw" "data/raw/census" ///
    "outputs" "outputs/audit" "outputs/reports" "outputs/labor" ///
    "outputs/preprocessing" "outputs/migration" {
    capture mkdir "`hroot'/`dir'"
}
* Immutable defaults; private local config overrides only source/options.
global H_CENSUS_KEY ""
global H_LAUS_SOURCE_MODE "mirror"
global H_JOINT_PLACEBO_REPS 199
global H_JOINT_PLACEBO_SEED 20261114
global H_JOINT_PLACEBO_ITER 30000
global H_JOINT_PLACEBO_TOL 0.00005
global H_JOINT_PLACEBO_STRICT 1
global H_RUN_CENSUS 1
global H_CENSUS_FILE ""
capture confirm file "`hroot'/config.local.do"
if !_rc do "`hroot'/config.local.do"
local panel "`hroot'/data/input/Helene_NC_ES_preinference_panel.dta"
local stagecsv "`hroot'/outputs/audit/STAGE_STATUS.csv"
file open hs using "`stagecsv'", write text replace
file write hs "stage,status,return_code" _n
file close hs
local any_fail 0
local haspanel 0
local primary_ok 0

* Always rebuild in-repository, NEVER import a prior inferred panel.
capture erase "`panel'"
di as result "STAGE 1/5: source records -> national panel (no synthetic fitting)"
capture noisily do "`hroot'/src/build_panel_from_sources.do"
local pre_rc=_rc
file open hs using "`stagecsv'", write text append
if `pre_rc' {
    file write hs "preinference,FAILED,`pre_rc'" _n
    local any_fail 1
    di as error "SOURCE BUILD FAILED r(`pre_rc'): downstream inference blocked."
}
else file write hs "preinference,PASSED,0" _n
file close hs
if `pre_rc'==0 {
    capture noisily shell /bin/sh "`hroot'/scripts/raw_sha256.sh" ///
        "`hroot'/data/raw" "`hroot'/outputs/audit/SOURCE_SHA256.txt"
}
if `pre_rc'==0 {
    capture noisily use "`panel'", clear
    local panel_rc=_rc
    if `panel_rc'==0 {
        capture isid fips mdate
        local panel_rc=_rc
    }
    if `panel_rc'==0 {
        capture assert _N==3221*57
        local panel_rc=_rc
    }
    if `panel_rc'==0 {
        quietly count if treated_nc==1 & mdate==ym(2024,8)
        if r(N)!=39 local panel_rc=459
    }
    if `panel_rc'==0 {
        capture assert missing(laus_employed) & missing(laus_unemployed) & ///
            missing(laus_labor_force) & missing(laus_unemployment_rate) ///
            if mdate==ym(2025,10)
        local panel_rc=_rc
    }
    if `panel_rc'==0 {
        quietly count if joint_eligible_donor==1 & mdate==ym(2024,8)
        if r(N)<500 local panel_rc=459
    }
    file open hs using "`stagecsv'", write text append
    if `panel_rc' {
        file write hs "panel_integrity,FAILED,`panel_rc'" _n
        local any_fail 1
    }
    else {
        file write hs "panel_integrity,PASSED,0" _n
        local haspanel 1
    }
    file close hs
}
else {
    file open hs using "`stagecsv'", write text append
    file write hs "panel_integrity,SKIPPED_SOURCE_BUILD_FAILED,`pre_rc'" _n
    file close hs
}

if `haspanel' {
    di as result "STAGE 2/5: ONLY four-outcome joint SCM + recovery inference"
    cd "`hroot'/src/labor"
    capture noisily do "run_multicontrol.do"
    local rc=_rc
    cd "`hroot'"
    file open hs using "`stagecsv'", write text append
    if `rc' {
        file write hs "joint_multioutcome,FAILED,`rc'" _n
        local any_fail 1
    }
    else {
        file write hs "joint_multioutcome,PASSED,0" _n
        local primary_ok 1
    }
    file close hs
    if `primary_ok' {
        di as result "STAGE 3/5: refitted four-outcome geographic placebos"
        cd "`hroot'/src/labor"
        capture noisily do "run_joint_placebos.do"
        local rc=_rc
        cd "`hroot'"
        file open hs using "`stagecsv'", write text append
        if `rc' {
            file write hs "joint_placebos,FAILED,`rc'" _n
            local any_fail 1
        }
        else file write hs "joint_placebos,PASSED,0" _n
        file close hs
    }
    else {
        file open hs using "`stagecsv'", write text append
        file write hs "joint_placebos,SKIPPED_PRIMARY_FAILED,0" _n
        file close hs
    }
    if ${H_RUN_CENSUS}==1 {
        di as result "STAGE 4/5: annual Census migration (exploratory)"
        capture noisily do "`hroot'/src/migration/01_census_v2025.do"
        local rc=_rc
        file open hs using "`stagecsv'", write text append
        if `rc' {
            file write hs "census_annual_migration,FAILED_OR_UNAVAILABLE,`rc'" _n
            local any_fail 1
        }
        else file write hs "census_annual_migration,PASSED,0" _n
        file close hs
    }
    else {
        file open hs using "`stagecsv'", write text append
        file write hs "census_annual_migration,SKIPPED_BY_CONFIG,0" _n
        file close hs
    }

}
else {
    file open hs using "`stagecsv'", write text append
    file write hs "downstream,SKIPPED_MISSING_OR_INVALID_PANEL,601" _n
    file close hs
    local any_fail 1
}
cd "`hroot'"
di as result "STAGE 5/5: regenerate consolidated TXT and stage audit"
capture noisily do "src/reporting/generate_report.do"
local report_rc=_rc
if `report_rc' {
    di as error "Text report failed: r(`report_rc')."
    local any_fail 1
}
set graphics on
if `any_fail' {
    di as error "SOME STAGES FAILED. Review outputs/audit/STAGE_STATUS.csv and HELENE_COMPLETE_RESULTS.txt"
    exit 459
}
di as result "COMPLETE: outputs/reports/HELENE_COMPLETE_RESULTS.txt"
