/* 03: BACKDATED OUT-OF-SAMPLE AUDIT.
   Slimmed 19 Nov 2026.

   The definitive weights are fitted over the whole pre-event window, so their
   pre-event fit is in-sample and is not evidence. It is therefore no longer
   computed here; the event-study pre-event tests in stage 04 serve that role.

   What remains is the one genuinely out-of-sample check in the design:
   backdated weights (same objective, same lambda, same intercept rule, fitted
   on 2022-2023 only) are used to predict January-August 2024, which they never
   saw. A ratio of forecast to training error near one would mean donors chosen
   on two years kept tracking the treated counties in the third.

   Removed 19 Nov 2026: the in-sample fit columns (fit_*) and the four in-time
   placebo contrasts. The latter compared three-month windows inside 2024 with
   earlier reference months; the Helene contrast exceeded the largest of them
   by only about 1.6x, so the check added little next to the geographic
   placebos and the null in the least-damaged tercile. Previous version:
   archive/03_before_slim_2026-11-19.do.txt */
version 16.0
di as txt "03: backdated out-of-sample audit (2022-23 weights predicting Jan-Aug 2024)"
use "${MC_OUT}/data/02_paths.dta", clear
isid mdate
local tr0=ym(2022,1)
local tr1=ym(2023,12)
local va0=ym(2024,1)
local va1=ym(2024,8)
tempname D
postfile `D' str12 outcome double bd_train_rmse bd_holdout_rmse bd_rmspe_ratio ///
    bd_holdout_max ///
    using "${MC_OUT}/data/03_fit_audit.dta", replace
foreach o in emp lf jobs urate {
    quietly gen double mc_gbd=mc_y_`o'-mc_syn_bd_`o'
    quietly gen double mc_sqb=mc_gbd^2
    quietly summarize mc_sqb if inrange(mdate,`tr0',`tr1'), meanonly
    local tr=sqrt(r(mean))
    quietly summarize mc_sqb if inrange(mdate,`va0',`va1'), meanonly
    local va=sqrt(r(mean))
    quietly summarize mc_gbd if inrange(mdate,`va0',`va1')
    local mx=max(abs(r(min)),abs(r(max)))
    post `D' ("`o'") (`tr') (`va') (`va'/`tr') (`mx')
    drop mc_gbd mc_sqb
}
postclose `D'
use "${MC_OUT}/data/03_fit_audit.dta", clear
assert !missing(bd_train_rmse,bd_holdout_rmse,bd_rmspe_ratio,bd_holdout_max)
assert _N==4
isid outcome
export delimited using "${MC_OUT}/tables/pretrend_train_holdout_ALL.csv", replace
list, noobs
* The ratio is the least interpretable form of this number: it is inflated by a
* very tight training fit. What matters for the design is (a) whether the
* out-of-sample gap TRENDS, which is tested on the event-study coefficients in
* stage 04d, and (b) its absolute size next to the estimated effects. Both are
* reported there; the absolute forecast error is printed here in display units.
quietly summarize bd_holdout_rmse, meanonly
di as txt "Out-of-sample forecast error, Jan-Aug 2024 (display units, log-% or pp):"
list outcome bd_holdout_rmse bd_holdout_max, noobs
di as txt "Compare these with the estimated effects before reading the ratio column."
quietly count if bd_rmspe_ratio>2
if r(N)>0 di as txt "NOTE: ratio exceeds 2 for `r(N)' outcome(s); see stage 04d for the trend-vs-amplitude split."
