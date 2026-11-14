/* Resume after successful 02_fit_weights.do, without repeating optimization. */
version 16.0
clear all
set more off
set graphics off
global MC_ROOT "`c(pwd)'"
if "${H_ROOT}"=="" global H_ROOT "${MC_ROOT}/../.."
global MC_OUT "${H_ROOT}/outputs/labor"
capture confirm file "config.local.do"
if !_rc do "config.local.do"
else do "config.example.do"
capture confirm file "${MC_OUT}/data/02_weighted_county_month.dta"
if _rc {
    di as error "Run run_multicontrol.do first to produce donor weights."
    exit 601
}
* Stata can retain a log from an earlier run; close it before starting this one.
capture log close _all
log using "${MC_OUT}/logs/multicontrol_postfit.log", replace text
do "steps/03_prefit_and_backdated_audit.do"
do "steps/04_joint_event_study.do"
do "steps/04b_extended_laus.do"
do "steps/05_tables_figures.do"
do "steps/06_kernel_rebased.do"
log close
set graphics on
