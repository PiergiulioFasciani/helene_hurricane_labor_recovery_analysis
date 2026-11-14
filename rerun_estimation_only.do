/* Re-estimate from the EXISTING source-built panel: stages 2-3 of run_all.do
   (joint synthetic-DiD + event study + kernel, then placebo refits) and the
   report. Does not rebuild data/input from raw sources; use run_all.do for a
   full source-to-results run. Writes outputs/audit/RERUN_STATUS.txt. */
version 16.0
clear all
set more off
local hroot "`c(pwd)'"
global H_ROOT "`hroot'"
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
local rc3=.
if `rc1'==0 & `rc2'==0 {
    capture noisily do "src/reporting/generate_report.do"
    local rc3=_rc
}
file open rs using "`statusfile'", write text replace
file write rs "joint_multioutcome rc=`rc1'" _n "joint_placebos rc=`rc2'" _n "report rc=`rc3'" _n
file close rs
set graphics on
di as result "RERUN FINISHED: multioutcome rc=`rc1', placebos rc=`rc2', report rc=`rc3'"
