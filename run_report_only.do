/* Regenerate readable results and all figures without refitting models. */
version 16.0
set more off
global H_ROOT "`c(pwd)'"
/* Only for the optional H_PYTHON interpreter override; no estimator setting
   in config.local.do affects reporting. */
capture confirm file "config.local.do"
if !_rc do "config.local.do"
do "src/reporting/generate_report.do"
