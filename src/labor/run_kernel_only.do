/* Resume just final Gaussian-kernel stage after successful stages 01-05.
   Does NOT re-fit donors or rerun primary event studies. Stata 16+. */
version 16.0
clear all
set more off
set graphics off
set linesize 120

global MC_ROOT "`c(pwd)'"
if "${H_ROOT}"=="" global H_ROOT "${MC_ROOT}/../.."
global MC_OUT "${H_ROOT}/outputs/labor"
capture confirm file "config.local.do"
if !_rc do "config.local.do"
else do "config.example.do"

capture confirm file "${MC_OUT}/data/02_weighted_county_month.dta"
if _rc {
    di as error "Missing fitted county-month data. Run run_multicontrol.do first."
    exit 601
}
capture confirm file "${MC_INPUT}"
if _rc {
    di as error "Source input file does not exist: ${MC_INPUT}"
    exit 601
}
capture mkdir "${MC_OUT}"
capture mkdir "${MC_OUT}/logs"
capture mkdir "${MC_OUT}/data"
capture mkdir "${MC_OUT}/tables"
capture mkdir "${MC_OUT}/figures"
capture log close _all
log using "${MC_OUT}/logs/multicontrol_kernel_only.log", text replace

do "steps/06_kernel_rebased.do"

di as result "KERNEL COMPLETE. Results under ${MC_OUT}/"
log close
set graphics on
