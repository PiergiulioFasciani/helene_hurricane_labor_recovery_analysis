/* Re-estimate from the EXISTING source-built panel: stages 2-3 of run_all.do
   (joint synthetic-DiD + event study, then placebo refits), the annual
   migration-by-damage stage, and the report. Does not rebuild migration
   source data: run src/migration/01_census_v2025.do or run_all.do for that.
   report. Does not rebuild data/input from raw sources; use run_all.do for a
   full source-to-results run. Writes outputs/audit/RERUN_STATUS.txt. */
version 16.0
clear all
set more off
local hroot "`c(pwd)'"
global H_ROOT "`hroot'"

* Start a fresh regression log: every estimation command appends its own
* expanded specification to this file as it runs (see src/labor/steps/_logspec.doh).
capture mkdir "${H_ROOT}/outputs/audit"
tempname _rl
file open `_rl' using "${H_ROOT}/outputs/audit/REGRESSIONS.txt", write replace text
file write `_rl' "HELENE PIPELINE -- REGRESSIONS AS ACTUALLY ESTIMATED" _n
file write `_rl' "Written by the estimation stages themselves on $S_DATE $S_TIME." _n
file write `_rl' "Map and interpretation: docs/REGRESSIONS.md" _n
file write `_rl' "" _n
file close `_rl'

capture confirm file "`hroot'/data/input/Helene_NC_ES_preinference_panel.dta"
if _rc {
    di as error "No source-built panel. Run run_all.do first."
    exit 601
}
global H_JOINT_PLACEBO_REPS 199
global H_JOINT_PLACEBO_SEED 20261114
global H_JOINT_PLACEBO_ITER 30000
global H_JOINT_PLACEBO_TOL 0.00005
global H_JOINT_PLACEBO_STRICT 1
capture confirm file "`hroot'/config.local.do"
if !_rc do "`hroot'/config.local.do"
local statusfile "`hroot'/outputs/audit/RERUN_STATUS.txt"
capture erase "`statusfile'"
cd "`hroot'/src/labor"
capture noisily do "run_multicontrol.do"
local rc1=_rc
capture log close _all
local rc2=.
if `rc1'==0 {
    cd "`hroot'/src/labor"
    capture noisily do "run_joint_placebos.do"
    local rc2=_rc
    capture log close _all
}
cd "`hroot'"
* Migration by damage tercile: consumes the tercile assignment written by 04c,
* so it must follow the labor stages and precede the report.
local rc3=.
if `rc1'==0 {
    capture noisily do "`hroot'/src/migration/02_migration_dose.do"
    local rc3=_rc
    capture log close _all
    if `rc3'==0 {
        * Register the stage so the report will render it.
        tempname sc
        file open `sc' using "`hroot'/outputs/audit/STAGE_STATUS.csv", read text
        local hasrow 0
        file read `sc' line
        while r(eof)==0 {
            if strpos("`line'","migration_by_damage")>0 local hasrow 1
            file read `sc' line
        }
        file close `sc'
        if `hasrow'==0 {
            file open `sc' using "`hroot'/outputs/audit/STAGE_STATUS.csv", write text append
            file write `sc' "migration_by_damage,PASSED,0" _n
            file close `sc'
        }
    }
}
local rc4=.
if `rc1'==0 & `rc2'==0 {
    capture noisily do "src/reporting/generate_report.do"
    local rc4=_rc
}
file open rs using "`statusfile'", write text replace
file write rs "joint_multioutcome rc=`rc1'" _n "joint_placebos rc=`rc2'" _n ///
    "migration_by_damage rc=`rc3'" _n "report rc=`rc4'" _n
file close rs
set graphics on
di as result "RERUN FINISHED: multioutcome rc=`rc1', placebos rc=`rc2', migration rc=`rc3', report rc=`rc4'"
