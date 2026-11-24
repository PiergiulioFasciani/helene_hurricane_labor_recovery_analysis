/* Reproducible multi-outcome control-group redesign, Stata 16+.
   Run from the unzipped repository root. This never edits its source input. */
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

* Override any old local paths in a unified source-first run.
if "${H_ROOT}"!="" global MC_INPUT "${H_ROOT}/data/input/Helene_NC_ES_preinference_panel.dta"
capture mkdir "${MC_OUT}"
capture mkdir "${MC_OUT}/data"
capture mkdir "${MC_OUT}/tables"
capture mkdir "${MC_OUT}/figures"
capture mkdir "${MC_OUT}/logs"
* Stata can retain a log from an earlier run; close it before starting this one.
capture log close _all
log using "${MC_OUT}/logs/multicontrol_run.log", text replace

di as result "Definitive four-outcome joint estimation: ${MC_INPUT}"
capture confirm file "${MC_INPUT}"
if _rc {
    di as error "INPUT NOT FOUND: ${MC_INPUT}"
    di as error "Copy config.example.do to config.local.do and set MC_INPUT."
    exit 601
}

do "steps/01_frozen_sample.do"
do "steps/02_fit_weights.do"
do "steps/03_prefit_and_backdated_audit.do"
do "steps/04_joint_event_study.do"
do "steps/04b_extended_laus.do"
do "steps/04c_dose_event_study.do"
do "steps/04d_backdated_weights_robustness.do"
do "steps/05_tables_figures.do"

di as result "COMPLETE. Results: ${MC_OUT}/"
log close
set graphics on
