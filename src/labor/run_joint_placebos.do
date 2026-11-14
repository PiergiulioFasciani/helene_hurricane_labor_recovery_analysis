/* REFITTED multi-outcome spatial group placebos, independent of legacy placebo.
   Requires stage 02 data, does not modify weight training or any primary result. */
version 16.0
set more off
set graphics off
capture log close _all
global MC_ROOT "`c(pwd)'"
if "${H_ROOT}"=="" global H_ROOT "${MC_ROOT}/../.."
global MC_OUT "${H_ROOT}/outputs/labor"
capture confirm file "${MC_OUT}/data/02_weighted_county_month.dta"
if _rc {
    di as error "Missing joint fitted panel. Run run_multicontrol.do first."
    exit 601
}
capture mkdir "${MC_OUT}/placebos"
local default_reps 49
local default_seed 20261114
local default_iter 30000
local default_tol 0.00005
local default_strict 1
if "${H_JOINT_PLACEBO_REPS}"=="" global H_JOINT_PLACEBO_REPS `default_reps'
if "${H_JOINT_PLACEBO_SEED}"=="" global H_JOINT_PLACEBO_SEED `default_seed'
if "${H_JOINT_PLACEBO_ITER}"=="" global H_JOINT_PLACEBO_ITER `default_iter'
if "${H_JOINT_PLACEBO_TOL}"=="" global H_JOINT_PLACEBO_TOL `default_tol'
if "${H_JOINT_PLACEBO_STRICT}"=="" global H_JOINT_PLACEBO_STRICT `default_strict'
log using "${MC_OUT}/placebos/joint_placebos.log", text replace name(jplog)
di as result "MULTI-OUTCOME PLACEBOS: refit on four Jan22-Aug24 objectives with intercept; B=${H_JOINT_PLACEBO_REPS}"
do "steps/07_joint_placebo_refits.do"
log close jplog
set graphics on
