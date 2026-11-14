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

* Only the definitive joint gap is displayed; no other control estimators.
use "${MC_OUT}/data/02_paths.dta", clear
foreach o in emp lf jobs urate {
    gen double joint_gap=mc_y_`o'-mc_syn_joint_`o'
    twoway (line joint_gap mdate, lcolor(navy) lwidth(medthick)), ///
        xline(`=ym(2024,9)', lcolor(maroon) lpattern(dash)) ///
        yline(0, lcolor(gs10)) xlabel(, angle(45)) ///
        title("Treated NC minus definitive joint synthetic: `o'") ///
        subtitle("One common four-outcome weight vector") ///
        xtitle("Month") ytitle("Gap after fitted intercept (log-% / pp)") ///
        note("Weights fitted Jan22-Aug24 with intercept; pre-event gap averages zero") ///
        graphregion(color(white))
    graph export "${MC_OUT}/figures/pretrend_and_recovery_`o'.png", replace width(2100)
    graph save "${MC_OUT}/figures/pretrend_and_recovery_`o'.gph", replace
    drop joint_gap
}

use "${MC_OUT}/data/04_event_coefficients.dta", clear
gen double db=cond(outcome=="mc_yurate",beta,100*beta)
gen double dlo=cond(outcome=="mc_yurate",lower95,100*lower95)
gen double dhi=cond(outcome=="mc_yurate",upper95,100*upper95)
foreach y in mc_yemp mc_ylf mc_yjobs mc_yurate mc_ydiff mc_yshare {
    twoway (rcap dhi dlo event_month if model=="joint" & outcome=="`y'", ///
        lcolor(gs10)) ///
        (connected db event_month if model=="joint" & outcome=="`y'", ///
        mcolor(navy) lcolor(navy) msymbol(O) msize(vsmall)), ///
        xline(0, lcolor(maroon) lpattern(dash)) yline(0,lcolor(gs8)) ///
        title("Multi-outcome joint-weight TWFE: `y'") ///
        subtitle("Relative to pre-event mean (Jan22-Aug24); 95% county-clustered CIs") ///
        xtitle("Months from September 2024") ///
        ytitle("Effect (log points x100; rate in pp)") ///
        legend(off) graphregion(color(white)) ///
        note("Synthetic-DiD weights fitted on Jan22-Aug24. Oct 2025 LAUS unavailable.")
    graph export "${MC_OUT}/figures/twfe_joint_`y'.png", replace width(2100)
    graph save "${MC_OUT}/figures/twfe_joint_`y'.gph", replace
}

* Create national panel with new donor assignments WITHOUT altering the input.
* Remaining US counties receive zero new donor weight, not missing labor data.
use "${MC_INPUT}", clear
isid fips mdate
merge m:1 fips using "${MC_OUT}/data/02_donor_weights.dta", ///
    keep(master match) gen(mc_merge)
gen byte mc_new_donor=(mc_merge==3)
assert mc_new_donor==0 if treated_nc==1
foreach w in mc_w_joint {
    replace `w'=0 if missing(`w')
}
assert !(treated_nc & mc_new_donor)
drop mc_merge
label var mc_w_joint "Joint emp/LF/QCEW/UR synthetic-DiD weight, Jan22-Aug24 fit"
label var mc_new_donor "Strict geographic donor with complete four-outcome 2022-23"
format mdate %tm
isid fips mdate
assert _N==3221*57
compress
save "${MC_OUT}/data/Helene_NC_ES_MULTIOUTCOME_national_panel.dta", replace

di as result "05 complete: multiplicity-adjusted tables, event-study figures, and extended national inference panel."
