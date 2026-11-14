/* Secondary longer followup through September 2026, LAUS outcomes ONLY.
   The workplace-job comparison cannot be carried past March 2026.
   This extended family is explicitly exploratory, not one of the 18
   matched primary tests with Holm adjustment.
   Revision 14 Nov 2026: like stage 04, coefficients and window effects are
   expressed relative to the pre-event mean (Jan 2022-Aug 2024). */
version 16.0
di as txt "04b: Additional LAUS-only followup, April to September 2026"
use "${MC_OUT}/data/02_weighted_county_month.dta", clear
keep if mdate!=ym(2024,9) & mdate!=ym(2025,10)
assert inrange(mdate,ym(2022,1),ym(2026,9))
assert !missing(mc_raw_emp,mc_raw_lf,mc_raw_ur,mc_raw_unemp)
by fips, sort: assert _N==55
isid fips mdate

gen double mx_emp=mc_raw_emp/100
gen double mx_lf=mc_raw_lf/100
gen double mx_urate=mc_raw_ur
gen double mx_share=mx_emp-mx_lf
gen double mx_unemp=mc_raw_unemp
gen int mx_k=mdate-ym(2024,9)
local events ""
local pre ""
local valid ""
local presum ""
forvalues k=-32/24 {
    if `k'!=-1 & `k'!=0 & `k'!=13 {
        if `k'<0 local s="m`=abs(`k')'"
        else local s="p`k'"
        gen byte mx_ev_`s'=(treated_nc==1 & mx_k==`k')
        local events "`events' mx_ev_`s'"
        if `k'<-1 local pre "`pre' mx_ev_`s'"
        if inrange(`k',-8,-2) local valid "`valid' mx_ev_`s'"
        if `k'<-1 {
            if "`presum'"=="" local presum "mx_ev_`s'"
            else local presum "`presum' + mx_ev_`s'"
        }
    }
}
local premean "(`presum')/32"
local slopeexpr ""
forvalues k=-2(-1)-32 {
    local a=abs(`k'+16.5)
    local s="m`=abs(`k')'"
    if `k'+16.5>0 {
        if "`slopeexpr'"=="" local slopeexpr "`a'*mx_ev_`s'"
        else local slopeexpr "`slopeexpr' + `a'*mx_ev_`s'"
    }
    else local slopeexpr "`slopeexpr' - `a'*mx_ev_`s'"
}

tempname E P D
postfile `E' str12 model str12 outcome int event_month ///
    double beta se lower95 upper95 pvalue long N_clust N_obs ///
    using "${MC_OUT}/data/04b_laus_event_coefficients.dta", replace
postfile `P' str12 model str12 outcome str15 window byte months ///
    double beta se lower95 upper95 pvalue lf_one_sided_lb ///
    byte rule_out_1pct rule_out_2pct long N_clust N_obs ///
    using "${MC_OUT}/data/04b_laus_period_effects.dta", replace
postfile `D' str12 model str12 outcome double pretrend_p last8_p slope_p ///
    using "${MC_OUT}/data/04b_laus_pretrends.dta", replace

local model "joint"
foreach y in mx_emp mx_lf mx_urate mx_share mx_unemp {
    local wvar="mc_w_joint"
        tempvar ww included
        gen double `ww'=cond(treated_nc==1,1/39,`wvar')
        gen byte `included'=treated_nc==1 | (mc_donor==1 & `wvar'>1e-10)
        quietly areg `y' `events' i.mdate if `included' [aw=`ww'], ///
            absorb(county_id) vce(cluster county_id)
        local df=e(df_r)
        local crit=invttail(`df',.025)
        local pre_p=.
        local val_p=.
        capture quietly test `pre'
        if !_rc local pre_p=r(p)
        capture quietly test `valid'
        if !_rc local val_p=r(p)
        quietly lincom `slopeexpr'
        local slope_p=r(p)
        post `D' ("`model'") ("`y'") (`pre_p') (`val_p') (`slope_p')
        forvalues k=-32/24 {
            if `k'!=0 & `k'!=13 {
                if `k'<0 local s="m`=abs(`k')'"
                else local s="p`k'"
                if `k'==-1 quietly lincom 0*mx_ev_m2 - `premean'
                else quietly lincom mx_ev_`s' - `premean'
                local b=r(estimate)
                local se=r(se)
                local pp=2*ttail(`df',abs(`b'/`se'))
                post `E' ("`model'") ("`y'") (`k') (`b') (`se') ///
                    (`b'-`crit'*`se') (`b'+`crit'*`se') (`pp') (e(N_clust)) (e(N))
            }
        }
        foreach window in late_extended aprsep2026 {
            local terms ""
            if "`window'"=="late_extended" {
                forvalues k=14/24 {
                    local terms "`terms' mx_ev_p`k'"
                }
            }
            else {
                forvalues k=19/24 {
                    local terms "`terms' mx_ev_p`k'"
                }
            }
            local nterms : word count `terms'
            local expr=""
            foreach v of local terms {
                if "`expr'"=="" local expr="`v'"
                else local expr="`expr' + `v'"
            }
            quietly lincom (`expr')/`nterms' - `premean'
            local b=r(estimate)
            local se=r(se)
            local pp=r(p)
            local lo=`b'-`crit'*`se'
            local hi=`b'+`crit'*`se'
            local lower1=.
            local rule1=.
            local rule2=.
            if "`y'"=="mx_lf" {
                local lower1=`b'-invttail(`df',.05)*`se'
                local rule1=(`lower1'>ln(.99))
                local rule2=(`lower1'>ln(.98))
            }
            post `P' ("`model'") ("`y'") ("`window'") (`nterms') ///
                (`b') (`se') (`lo') (`hi') (`pp') (`lower1') (`rule1') (`rule2') ///
                (e(N_clust)) (e(N))
        }
        drop `ww' `included'
}
postclose `E'
postclose `P'
postclose `D'

use "${MC_OUT}/data/04b_laus_event_coefficients.dta", clear
isid model outcome event_month
export delimited using "${MC_OUT}/tables/LAUS_extended_event_coefficients_EXPLORATORY.csv", replace
use "${MC_OUT}/data/04b_laus_period_effects.dta", clear
isid model outcome window
export delimited using "${MC_OUT}/tables/LAUS_extended_period_effects_EXPLORATORY.csv", replace
use "${MC_OUT}/data/04b_laus_pretrends.dta", clear
export delimited using "${MC_OUT}/tables/LAUS_extended_pretrends.csv", replace

di as result "04b complete: LAUS alone through Sep 2026 (exploratory extension)."
