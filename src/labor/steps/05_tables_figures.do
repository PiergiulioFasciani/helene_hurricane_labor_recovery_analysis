/* Transparent publication-friendly outputs. 18 prespecified joint-weight
   tests (six outcomes x three windows) are adjusted using Holm's step-down
   family-wise-error procedure; all raw p-values are retained. */
version 16.0
di as txt "05: Prespecified multiplicity correction, graphs and national panel"
use "${MC_OUT}/data/04_period_effects.dta", clear
preserve
    keep if model=="joint"
    quietly count
    assert r(N)==18
    sort pvalue outcome window
    gen double holm_step=min(1,(19-_n)*pvalue)
    gen double holm_p=holm_step
    forvalues j=2/18 {
        quietly replace holm_p=max(holm_p,holm_p[_n-1]) in `j'
    }
    keep model outcome window holm_p
    isid model outcome window
    tempfile holm
    save `holm'
restore
merge 1:1 model outcome window using `holm', keep(master match) nogen

gen double display_beta=cond(outcome=="mc_yurate", beta,100*beta)
gen double display_lo=cond(outcome=="mc_yurate",lower95,100*lower95)
gen double display_hi=cond(outcome=="mc_yurate",upper95,100*upper95)
gen str7 display_units="logpct"
replace display_units="pp" if outcome=="mc_yurate"
label var holm_p "Holm-adjusted p, joint-weight primary family only (18 tests)"
order model outcome window months display_beta display_lo display_hi ///
    display_units pvalue holm_p beta se lower95 upper95
save "${MC_OUT}/data/05_reported_period_effects.dta", replace
export delimited using "${MC_OUT}/tables/period_effects_WITH_HOLM.csv", replace

* Only these 18 joint-primary hypotheses form the multiplicity family.
list outcome window display_beta display_lo display_hi pvalue holm_p ///
    if model=="joint", noobs sepby(outcome)

* Figures are rendered by src/reporting/make_report.py from the exported CSVs.
* Stage 05 no longer writes Stata .png/.gph duplicates (removed 18 Nov 2026):
* nothing read them and they were regenerated on every run.

use "${MC_OUT}/data/04_event_coefficients.dta", clear
* The merged national panel written here until 18 Nov 2026 was never read by
* any stage, test or report, and cost 47MB per run. Donor weights are already
* distributed as outputs/labor/tables/donor_weights.csv; merge them onto the
* source panel directly if a national file is ever needed.

di as result "05 complete: multiplicity-adjusted tables (figures are rendered by the reporting stage)."
