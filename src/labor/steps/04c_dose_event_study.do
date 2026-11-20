/* 04c: DAMAGE-TERCILE EVENT STUDY.
   Added 16 Nov 2026; reduced to a minimal design 18 Nov 2026.

   The 39 treated counties are split into terciles of storm exposure and each
   tercile gets its own event-study series against the same synthetic control.
   Sample, weights, fixed effects, clustering and the pre-event-mean reference
   are identical to stage 04; only the treated indicators are split.

   Two exposure measures:
     hdsi_tercile  PRIMARY. Composite damage index (rainfall rank and verified
                   loss relative to housing value).
     rain_tercile  ROBUSTNESS. Six-day rainfall alone, i.e. meteorology with no
                   FEMA claim behaviour in it. Reported as a single footnote
                   comparison, not as a second set of headline results.

   Each tercile effect is an ATT for that exposure group relative to never-
   treated donors, so it needs only the parallel-trends assumption of stage 04
   applied within the group. Exposure is not randomly assigned across the
   treated counties: damage correlates with terrain, income and industry mix,
   so the high-minus-low contrast is a comparison of differently exposed
   counties, not a randomised dose.

   Removed 18 Nov 2026: the continuous per-SD interaction families for HDSI
   and rainfall, and their average-exposure level terms. A per-SD slope is a
   causal response only under a stronger parallel-trends assumption than the
   binned contrasts need, and nothing downstream depended on it. The previous
   version is at archive/04c_before_minimal_2026-11-18.do.txt. */
version 16.0
di as txt "04c: damage-tercile event study (HDSI primary, rainfall robustness)"

tempfile dose
use "${MC_INPUT}", clear
keep if treated_nc==1 & mdate==ym(2024,8)
keep fips hdsi prcp_total_mm
isid fips
quietly count
assert r(N)==39
assert !missing(hdsi,prcp_total_mm)
xtile md_ter=hdsi, nq(3)
xtile md_rter=prcp_total_mm, nq(3)
label var md_ter "HDSI tercile, 1 least to 3 most damaged"
label var md_rter "Six-day rainfall tercile, 1 least to 3 most"
tab md_ter md_rter
export delimited fips hdsi prcp_total_mm md_ter md_rter ///
    using "${MC_OUT}/tables/dose_treated_exposure.csv", replace
save `dose'

use "${MC_OUT}/data/02_weighted_county_month.dta", clear
keep if inrange(mdate,ym(2022,1),ym(2026,3)) & ///
    mdate!=ym(2024,9) & mdate!=ym(2025,10)
assert !missing(mc_raw_emp,mc_raw_lf,mc_raw_jobs,mc_raw_ur)
by fips, sort: assert _N==49
isid fips mdate
merge m:1 fips using `dose', keep(master match) nogen
assert !missing(md_ter,md_rter) if treated_nc==1
foreach v in md_ter md_rter {
    replace `v'=0 if treated_nc==0
}

gen double mc_yemp=mc_raw_emp/100
gen double mc_ylf=mc_raw_lf/100
gen double mc_yjobs=mc_raw_jobs/100
gen double mc_yurate=mc_raw_ur
gen double mc_ydiff=mc_yjobs-mc_yemp
gen double mc_yshare=mc_yemp-mc_ylf
gen int md_k=mdate-ym(2024,9)

* Six indicator families: three HDSI terciles and three rainfall terciles.
local prefixes "md_t1 md_t2 md_t3 md_r1 md_r2 md_r3"
foreach pf of local prefixes {
    local all_`pf' ""
    local pre_`pf' ""
    local presum_`pf' ""
    local slope_`pf' ""
}
forvalues k=-32/18 {
    if `k'!=-1 & `k'!=0 & `k'!=13 {
        if `k'<0 local s="m`=abs(`k')'"
        else local s="p`k'"
        gen byte md_ev_`s'=(treated_nc==1 & md_k==`k')
        forvalues g=1/3 {
            gen byte md_t`g'_`s'=md_ev_`s'*(md_ter==`g')
            gen byte md_r`g'_`s'=md_ev_`s'*(md_rter==`g')
        }
        drop md_ev_`s'
        foreach pf of local prefixes {
            local all_`pf' "`all_`pf'' `pf'_`s'"
            if `k'<-1 {
                local pre_`pf' "`pre_`pf'' `pf'_`s'"
                if "`presum_`pf''"=="" local presum_`pf' "`pf'_`s'"
                else local presum_`pf' "`presum_`pf'' + `pf'_`s'"
                local a=abs(`k'+16.5)
                if `k'+16.5>0 {
                    if "`slope_`pf''"=="" local slope_`pf' "`a'*`pf'_`s'"
                    else local slope_`pf' "`slope_`pf'' + `a'*`pf'_`s'"
                }
                else local slope_`pf' "`slope_`pf'' - `a'*`pf'_`s'"
            }
        }
    }
}
foreach pf of local prefixes {
    local premean_`pf' "(`presum_`pf'')/32"
}

tempname E P D
postfile `E' str16 dose str12 outcome str16 term int event_month ///
    double beta se lower95 upper95 pvalue long N_clust N_obs ///
    using "${MC_OUT}/data/04c_dose_event_coefficients.dta", replace
postfile `P' str16 dose str12 outcome str16 term str12 window byte months ///
    double beta se lower95 upper95 pvalue long N_clust N_obs ///
    using "${MC_OUT}/data/04c_dose_period_effects.dta", replace
postfile `D' str16 dose str12 outcome str16 term double pretrend_p slope_p ///
    max_abs_pre int n_sig05 long N_clust N_obs ///
    using "${MC_OUT}/data/04c_dose_pretrend_tests.dta", replace

tempvar mw mins
gen double `mw'=cond(treated_nc==1,1/39,mc_w_joint)
gen byte `mins'=treated_nc==1 | (mc_w_joint>1e-10 & mc_donor==1)
assert `mw'>0 if `mins'

foreach y in mc_yemp mc_ylf mc_yjobs mc_yurate mc_ydiff mc_yshare {
    foreach fam in hdsi_tercile rain_tercile {
        if "`fam'"=="hdsi_tercile" {
            local rhs "`all_md_t1' `all_md_t2' `all_md_t3'"
            local terms "md_t1 md_t2 md_t3"
        }
        else {
            local rhs "`all_md_r1' `all_md_r2' `all_md_r3'"
            local terms "md_r1 md_r2 md_r3"
        }
        local labels "low mid high"
        di as txt "04c `fam': `y'"
        quietly areg `y' `rhs' i.mdate if `mins' [aw=`mw'], ///
            absorb(county_id) vce(cluster county_id)
        local spec_tag  "04c `fam' / `y'"
        local nev : word count `rhs'
        local spec_cmd  "areg `y' <`nev' tercile-by-event-month indicators> i.mdate if <sample> [aw=<weight>], absorb(county_id) vce(cluster county_id)"
        local spec_rhs  "`nev' regressors: the stage-04 event indicators split three ways by exposure tercile (low/mid/high, 13 treated counties each), so each tercile gets its own event path against the same donors; plus i.mdate"
        local spec_sample "identical to stage 04: Jan 2022-Mar 2026, excluding Sep 2024 and Oct 2025"
        local spec_wt   "identical to stage 04: treated=1/39, donor=mc_w_joint; analytic"
        local spec_ref  "each tercile series recentred to its own 32-month pre-event mean; high-minus-low contrast taken with joint covariance"
        local spec_out  "04c_dose_*.dta -> dose_event_coefficients.csv, dose_period_effects.csv"
        include "${MC_ROOT}/steps/_logspec.doh"
        local n=e(N)
        local g=e(N_clust)
        local df=e(df_r)
        local tcrit=invttail(`df',.025)
        forvalues j=1/3 {
            local pf : word `j' of `terms'
            local lab : word `j' of `labels'
            quietly test `pre_`pf''
            local pp=r(p)
            quietly lincom `slope_`pf''
            local sp=r(p)
            local maxb=0
            local nsig=0
            forvalues k=-32/18 {
                if `k'!=0 & `k'!=13 {
                    if `k'<0 local s="m`=abs(`k')'"
                    else local s="p`k'"
                    if `k'==-1 quietly lincom 0*`pf'_m2 - `premean_`pf''
                    else quietly lincom `pf'_`s' - `premean_`pf''
                    local b=r(estimate)
                    local se=r(se)
                    local pv=2*ttail(`df',abs(`b'/`se'))
                    post `E' ("`fam'") ("`y'") ("`lab'") (`k') (`b') (`se') ///
                        (`b'-`tcrit'*`se') (`b'+`tcrit'*`se') (`pv') (`g') (`n')
                    if `k'<0 {
                        if abs(`b')>`maxb' local maxb=abs(`b')
                        if `pv'<.05 local nsig=`nsig'+1
                    }
                }
            }
            post `D' ("`fam'") ("`y'") ("`lab'") (`pp') (`sp') ///
                (`maxb') (`nsig') (`g') (`n')
            foreach window in immediate year2025 late {
                if "`window'"=="immediate" local ks "1 2 3"
                if "`window'"=="year2025" local ks "4 5 6 7 8 9 10 11 12"
                if "`window'"=="late" local ks "14 15 16 17 18"
                local len : word count `ks'
                local expr ""
                foreach kk of local ks {
                    if "`expr'"=="" local expr "`pf'_p`kk'"
                    else local expr "`expr' + `pf'_p`kk'"
                }
                quietly lincom (`expr')/`len' - `premean_`pf''
                post `P' ("`fam'") ("`y'") ("`lab'") ("`window'") (`len') ///
                    (r(estimate)) (r(se)) (r(estimate)-`tcrit'*r(se)) ///
                    (r(estimate)+`tcrit'*r(se)) (r(p)) (`g') (`n')
            }
        }
        * High minus low, same window, from the joint covariance.
        local hi : word 3 of `terms'
        local lo : word 1 of `terms'
        foreach window in immediate year2025 late {
            if "`window'"=="immediate" local ks "1 2 3"
            if "`window'"=="year2025" local ks "4 5 6 7 8 9 10 11 12"
            if "`window'"=="late" local ks "14 15 16 17 18"
            local len : word count `ks'
            local eh ""
            local el ""
            foreach kk of local ks {
                if "`eh'"=="" {
                    local eh "`hi'_p`kk'"
                    local el "`lo'_p`kk'"
                }
                else {
                    local eh "`eh' + `hi'_p`kk'"
                    local el "`el' + `lo'_p`kk'"
                }
            }
            quietly lincom ((`eh')/`len' - `premean_`hi'') - ((`el')/`len' - `premean_`lo'')
            post `P' ("`fam'") ("`y'") ("high_minus_low") ("`window'") (`len') ///
                (r(estimate)) (r(se)) (r(estimate)-`tcrit'*r(se)) ///
                (r(estimate)+`tcrit'*r(se)) (r(p)) (`g') (`n')
        }
    }
}
postclose `E'
postclose `P'
postclose `D'

use "${MC_OUT}/data/04c_dose_event_coefficients.dta", clear
isid dose outcome term event_month
export delimited using "${MC_OUT}/tables/dose_event_coefficients.csv", replace

use "${MC_OUT}/data/04c_dose_period_effects.dta", clear
isid dose outcome term window
* Holm retained in the saved table for completeness; the report shows raw p.
gen double holm_p=.
foreach fam in hdsi_tercile rain_tercile {
    quietly count if dose=="`fam'" & term=="high_minus_low"
    assert r(N)==18
    gen double md_ord=pvalue if dose=="`fam'" & term=="high_minus_low"
    sort md_ord outcome window
    gen double md_step=min(1,(19-_n)*pvalue) if !missing(md_ord)
    forvalues j=2/18 {
        quietly replace md_step=max(md_step,md_step[_n-1]) in `j'
    }
    replace holm_p=md_step if !missing(md_ord)
    drop md_ord md_step
}
gen double display_beta=cond(outcome=="mc_yurate",beta,100*beta)
gen double display_lo=cond(outcome=="mc_yurate",lower95,100*lower95)
gen double display_hi=cond(outcome=="mc_yurate",upper95,100*upper95)
sort dose outcome term window
export delimited using "${MC_OUT}/tables/dose_period_effects.csv", replace
list dose outcome term window display_beta display_lo display_hi pvalue ///
    if dose=="hdsi_tercile", noobs sepby(outcome)

use "${MC_OUT}/data/04c_dose_pretrend_tests.dta", clear
isid dose outcome term
export delimited using "${MC_OUT}/tables/dose_pretrend_tests.csv", replace
list dose outcome term slope_p n_sig05 max_abs_pre if dose=="hdsi_tercile", noobs

di as result "04c complete: damage-tercile event study (HDSI primary, rainfall robustness)."
