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
local va ""
local presum ""
forvalues k=-32/18 {
    if `k'!=-1 & `k'!=0 & `k'!=13 {
        if `k'<0 local s="m`=abs(`k')'"
        else local s="p`k'"
        gen byte mc_ev_`s'=(treated_nc==1 & mc_k==`k')
        local events "`events' mc_ev_`s'"
        if `k'<-1 local pre "`pre' mc_ev_`s'"
        if inrange(`k',-8,-2) local va "`va' mc_ev_`s'"
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
* Half-year bin means equal to the Jan-Aug 2024 mean (4 restrictions).
* Bin j has 6 months; the 2024 bin has 8 (seven coefficients plus the zero).
local last ""
forvalues k=-8/-2 {
    local last "`last' - 6*mc_ev_m`=abs(`k')'"
}
forvalues j=1/4 {
    local lo=-32+6*(`j'-1)
    local hi=`lo'+5
    local bin`j' ""
    forvalues k=`lo'/`hi' {
        if "`bin`j''"=="" local bin`j' "8*mc_ev_m`=abs(`k')'"
        else local bin`j' "`bin`j'' + 8*mc_ev_m`=abs(`k')'"
    }
    local bin`j' "`bin`j'' `last'"
}

tempname E P D
postfile `E' str12 model str12 outcome int event_month double beta se ///
    lower95 upper95 double pvalue long N_clust N_obs ///
    using "${MC_OUT}/data/04_event_coefficients.dta", replace
postfile `P' str12 model str12 outcome str12 window byte months ///
    double beta se lower95 upper95 pvalue lf_one_sided_lb ///
    byte rule_out_1pct rule_out_2pct long N_clust N_obs ///
    using "${MC_OUT}/data/04_period_effects.dta", replace
postfile `D' str12 model str12 outcome double pretrend_p last8_p ///
    slope_p bin_p pre_rmse max_abs_pre max_abs_t int n_sig05 ///
    long N_clust N_obs using "${MC_OUT}/data/04_twfe_diagnostics.dta", replace

local model "joint"
foreach y in mc_yemp mc_ylf mc_yjobs mc_yurate mc_ydiff mc_yshare {
    local wvar="mc_w_joint"
        tempvar mw mins
        gen double `mw'=cond(treated_nc==1,1/39,`wvar')
        gen byte `mins'=treated_nc==1 | (`wvar'>1e-10 & mc_donor==1)
        assert `mw'>0 if `mins'
        di as txt "TWFE `model': `y'"
        quietly areg `y' `events' i.mdate if `mins' [aw=`mw'], ///
            absorb(county_id) vce(cluster county_id)
        local n=e(N)
        local g=e(N_clust)
        local df=e(df_r)
        local tcrit=invttail(`df',.025)
        local pp=.
        local vp=.
        capture quietly test `pre'
        if !_rc local pp=r(p)
        capture quietly test `va'
        if !_rc local vp=r(p)
        quietly lincom `slopeexpr'
        local sp=r(p)
        quietly test (`bin1'=0) (`bin2'=0) (`bin3'=0) (`bin4'=0)
        local bp=r(p)
        local ssq=0
        local maxb=0
        local maxt=0
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
                    local ssq=`ssq'+(`b')^2
                    if abs(`b')>`maxb' local maxb=abs(`b')
                    if abs(`b'/`se')>`maxt' local maxt=abs(`b'/`se')
                    if `pv'<.05 local nsig=`nsig'+1
                }
            }
        }
        post `D' ("`model'") ("`y'") (`pp') (`vp') (`sp') (`bp') ///
            (sqrt(`ssq'/32)) (`maxb') (`maxt') (`nsig') (`g') (`n')
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
            local lower1=.
            local r1=.
            local r2=.
            if "`y'"=="mc_ylf" {
                local lower1=`b'-invttail(`df',.05)*`s'
                local r1=(`lower1'>ln(.99))
                local r2=(`lower1'>ln(.98))
            }
            post `P' ("`model'") ("`y'") ("`window'") (`len') ///
                (`b') (`s') (`lb') (`ub') (`pv') (`lower1') ///
                (`r1') (`r2') (`g') (`n')
        }
        estimates save "${MC_OUT}/data/04_est_`model'_`y'.ster", replace
        drop `mw' `mins'
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
