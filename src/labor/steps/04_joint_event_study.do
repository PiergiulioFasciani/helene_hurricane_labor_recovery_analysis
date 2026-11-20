/* Synthetic-DiD event study (revision 14 Nov 2026).
   The regression is unchanged (county and month effects, treated x month
   indicators, August 2024 omitted, donor weights as analytic weights), but
   every reported coefficient and window effect is now expressed RELATIVE TO
   THE PRE-EVENT MEAN (Jan 2022-Aug 2024, uniform time weights), not to the
   single month of August 2024. Reported pre-event coefficients therefore
   average zero by construction and August 2024 gets its own estimate.
   Joint Wald tests are invariant to this re-expression.
   Main: common FOUR-OUTCOME weights.
   Same Jan22-Mar26 outcomes/units, Oct25 LAUS missing, Sep24 excluded.
   All six matched outcomes use identical rows for common-weight inference.
   Distinct 'specific' weights NEVER used to build workplace/resident contrast. */
version 16.0
di as txt "04: Matched-sample TWFE event study for jointly trained donor weights"
use "${MC_OUT}/data/02_weighted_county_month.dta", clear

* Exactly the common observed months for all four labor outcomes.
keep if inrange(mdate,ym(2022,1),ym(2026,3)) & ///
    mdate!=ym(2024,9) & mdate!=ym(2025,10)
assert !missing(mc_raw_emp,mc_raw_lf,mc_raw_jobs,mc_raw_ur)
by fips, sort: assert _N==49
isid fips mdate

* Raw outcomes. Differences are direct regressions with joint covariance.
gen double mc_yemp=mc_raw_emp/100
gen double mc_ylf=mc_raw_lf/100
gen double mc_yjobs=mc_raw_jobs/100
gen double mc_yurate=mc_raw_ur
gen double mc_ydiff=mc_yjobs-mc_yemp
gen double mc_yshare=mc_yemp-mc_ylf

gen int mc_k=mdate-ym(2024,9)
local events ""
local pre ""
local presum ""
forvalues k=-32/18 {
    if `k'!=-1 & `k'!=0 & `k'!=13 {
        if `k'<0 local s="m`=abs(`k')'"
        else local s="p`k'"
        gen byte mc_ev_`s'=(treated_nc==1 & mc_k==`k')
        local events "`events' mc_ev_`s'"
        if `k'<-1 local pre "`pre' mc_ev_`s'"
        if `k'<-1 {
            if "`presum'"=="" local presum "mc_ev_`s'"
            else local presum "`presum' + mc_ev_`s'"
        }
    }
}
* Pre-event mean of the 32 monthly gaps (Aug 2024 is the omitted zero).
local premean "(`presum')/32"

* Linear pre-trend contrast: sum over k=-32..-1 of (k+16.5)*beta_k.
local slopeexpr ""
forvalues k=-2(-1)-32 {
    local a=abs(`k'+16.5)
    local s="m`=abs(`k')'"
    if `k'+16.5>0 {
        if "`slopeexpr'"=="" local slopeexpr "`a'*mc_ev_`s'"
        else local slopeexpr "`slopeexpr' + `a'*mc_ev_`s'"
    }
    else local slopeexpr "`slopeexpr' - `a'*mc_ev_`s'"
}
tempname E P D
postfile `E' str12 model str12 outcome int event_month double beta se ///
    lower95 upper95 double pvalue long N_clust N_obs ///
    using "${MC_OUT}/data/04_event_coefficients.dta", replace
postfile `P' str12 model str12 outcome str12 window byte months ///
    double beta se lower95 upper95 pvalue long N_clust N_obs ///
    using "${MC_OUT}/data/04_period_effects.dta", replace
postfile `D' str12 model str12 outcome double pretrend_p slope_p ///
    max_abs_pre int n_sig05 long N_clust N_obs using "${MC_OUT}/data/04_twfe_diagnostics.dta", replace

local model "joint"
* Weights and sample are identical for every outcome, so build them once.
*   treated county  -> 1/39 (equal-county treated mean)
*   donor county    -> its fitted synthetic weight mc_w_joint
*   in sample       -> treated, or a donor carrying non-negligible weight
tempvar mw mins
gen double `mw'=cond(treated_nc==1,1/39,mc_w_joint)
gen byte `mins'=treated_nc==1 | (mc_w_joint>1e-10 & mc_donor==1)
assert `mw'>0 if `mins'

foreach y in mc_yemp mc_ylf mc_yjobs mc_yurate mc_ydiff mc_yshare {
        di as txt "TWFE `model': `y'"
        quietly areg `y' `events' i.mdate if `mins' [aw=`mw'], ///
            absorb(county_id) vce(cluster county_id)
        local spec_tag  "04 joint / `y'"
        local nev : word count `events'
        local spec_cmd  "areg `y' <`nev' event indicators> i.mdate if <sample> [aw=<weight>], absorb(county_id) vce(cluster county_id)"
        local spec_rhs  "`nev' treated-by-event-month indicators, k=-32..18 omitting k=-1 (Aug 2024, omitted reference), k=0 (Sep 2024, transition month) and k=13 (Oct 2025, LAUS missing); plus i.mdate calendar-month effects"
        local spec_sample "Jan 2022-Mar 2026, excluding Sep 2024 and Oct 2025; 39 treated counties plus donors with mc_w_joint>1e-10"
        local spec_wt   "treated=1/39 (equal-county mean); donor=mc_w_joint (fitted synthetic weight); analytic"
        local spec_ref  "coefficients recentred to the mean of the 32 pre-event months Jan2022-Aug2024, NOT to Aug 2024 alone"
        local spec_out  "04_event_coefficients.dta / 04_period_effects.dta -> event_coefficients.csv, period_effects_raw.csv"
        include "${MC_ROOT}/steps/_logspec.doh"
        local n=e(N)
        local g=e(N_clust)
        local df=e(df_r)
        local tcrit=invttail(`df',.025)
        local pp=.
        capture quietly test `pre'
        if !_rc local pp=r(p)
        quietly lincom `slopeexpr'
        local sp=r(p)
        local maxb=0
        local nsig=0
        forvalues k=-32/18 {
            if `k'!=0 & `k'!=13 {
                if `k'<0 local s="m`=abs(`k')'"
                else local s="p`k'"
                if `k'==-1 quietly lincom 0*mc_ev_m2 - `premean'
                else quietly lincom mc_ev_`s' - `premean'
                local b=r(estimate)
                local se=r(se)
                local lb=`b'-`tcrit'*`se'
                local ub=`b'+`tcrit'*`se'
                local pv=2*ttail(`df',abs(`b'/`se'))
                post `E' ("`model'") ("`y'") (`k') (`b') (`se') ///
                    (`lb') (`ub') (`pv') (`g') (`n')
                if `k'<0 {
                    if abs(`b')>`maxb' local maxb=abs(`b')
                    if `pv'<.05 local nsig=`nsig'+1
                }
            }
        }
        post `D' ("`model'") ("`y'") (`pp') (`sp') (`maxb') (`nsig') (`g') (`n')
        foreach window in immediate year2025 late {
            if "`window'"=="immediate" local terms "mc_ev_p1 mc_ev_p2 mc_ev_p3"
            if "`window'"=="year2025" local terms "mc_ev_p4 mc_ev_p5 mc_ev_p6 mc_ev_p7 mc_ev_p8 mc_ev_p9 mc_ev_p10 mc_ev_p11 mc_ev_p12"
            if "`window'"=="late" local terms "mc_ev_p14 mc_ev_p15 mc_ev_p16 mc_ev_p17 mc_ev_p18"
            local len : word count `terms'
            local expr ""
            foreach term of local terms {
                if "`expr'"=="" local expr "`term'"
                else local expr "`expr' + `term'"
            }
            quietly lincom (`expr')/`len' - `premean'
            local b=r(estimate)
            local s=r(se)
            local pv=r(p)
            local lb=`b'-`tcrit'*`s'
            local ub=`b'+`tcrit'*`s'
            post `P' ("`model'") ("`y'") ("`window'") (`len') ///
                (`b') (`s') (`lb') (`ub') (`pv') (`g') (`n')
        }
}
postclose `E'
postclose `P'
postclose `D'

use "${MC_OUT}/data/04_event_coefficients.dta", clear
isid model outcome event_month
export delimited using "${MC_OUT}/tables/event_coefficients.csv", replace
use "${MC_OUT}/data/04_period_effects.dta", clear
isid model outcome window
export delimited using "${MC_OUT}/tables/period_effects_raw.csv", replace
use "${MC_OUT}/data/04_twfe_diagnostics.dta", clear
isid model outcome
export delimited using "${MC_OUT}/tables/twfe_pretrend_tests.csv", replace

di as result "04 complete: joint primary ONLY (six outcomes), pre-event-mean reference."
