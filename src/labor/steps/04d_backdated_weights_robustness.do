/* 04d: THE MAIN EVENT STUDY REFITTED ON BACKDATED WEIGHTS.
   Added 22 Nov 2026.

   The definitive weights are fitted over the whole pre-event window, so their
   pre-event fit is in-sample and the backdated audit in stage 03 can only test
   the fitting PROCEDURE, not the weights actually used. This stage closes that
   gap directly: it re-estimates stage 04 and the HDSI-tercile contrast of
   stage 04c using mc_w_joint_bd, the weight vector fitted on 2022-2023 ONLY,
   which never saw January-August 2024.

   Everything else is identical to stages 04 and 04c: same sample window, same
   excluded months, same county and calendar-month effects, same clustering,
   same pre-event-mean recentring, same three windows. Only the donor weights
   and the donor inclusion rule change, both now keyed on the backdated vector.

   Read it as a robustness column, not as a replacement. The backdated weights
   are fitted on 24 months instead of 32, so they are handicapped by
   construction; if the damage gradient survives them it survives a control
   that could not have been tuned to the pre-storm year. */
version 16.0
di as txt "04d: robustness refit on 2022-23 (backdated) weights"

use "${MC_OUT}/data/02_weighted_county_month.dta", clear
capture confirm variable mc_w_joint_bd
if _rc {
    di as error "mc_w_joint_bd not found: rerun stage 02 so backdated weights are stored."
    exit 111
}
keep if inrange(mdate,ym(2022,1),ym(2026,3)) & ///
    mdate!=ym(2024,9) & mdate!=ym(2025,10)
assert !missing(mc_raw_emp,mc_raw_lf,mc_raw_jobs,mc_raw_ur)
isid fips mdate

* Damage terciles, defined exactly as in stage 04c.
tempfile dose
preserve
    use "${MC_INPUT}", clear
    keep if treated_nc==1 & mdate==ym(2024,8)
    keep fips hdsi
    isid fips
    quietly count
    assert r(N)==39
    assert !missing(hdsi)
    xtile md_ter=hdsi, nq(3)
    keep fips md_ter
    save `dose'
restore
merge m:1 fips using `dose', keep(master match) nogen
replace md_ter=0 if treated_nc==0
assert inrange(md_ter,1,3) if treated_nc==1

gen double mc_yemp=mc_raw_emp/100
gen double mc_ylf=mc_raw_lf/100
gen double mc_yjobs=mc_raw_jobs/100
gen double mc_yurate=mc_raw_ur
gen double mc_ydiff=mc_yjobs-mc_yemp
gen double mc_yshare=mc_yemp-mc_ylf
gen int mc_k=mdate-ym(2024,9)

* Aggregate indicators, and the same indicators split by tercile.
local events ""
local presum ""
local prefixes "bt1 bt2 bt3"
foreach pf of local prefixes {
    local all_`pf' ""
    local presum_`pf' ""
}
forvalues k=-32/18 {
    if `k'!=-1 & `k'!=0 & `k'!=13 {
        if `k'<0 local s="m`=abs(`k')'"
        else local s="p`k'"
        gen byte bd_ev_`s'=(treated_nc==1 & mc_k==`k')
        local events "`events' bd_ev_`s'"
        forvalues g=1/3 {
            gen byte bt`g'_`s'=bd_ev_`s'*(md_ter==`g')
            local all_bt`g' "`all_bt`g'' bt`g'_`s'"
        }
        if `k'<-1 {
            if "`presum'"=="" local presum "bd_ev_`s'"
            else local presum "`presum' + bd_ev_`s'"
            foreach pf of local prefixes {
                if "`presum_`pf''"=="" local presum_`pf' "`pf'_`s'"
                else local presum_`pf' "`presum_`pf'' + `pf'_`s'"
            }
        }
    }
}
local premean "(`presum')/32"

* Linear pre-trend contrast, as in stage 04: sum over k=-32..-2 of (k+16.5)*beta_k.
local slopeexpr ""
forvalues k=-2(-1)-32 {
    local a=abs(`k'+16.5)
    local s="m`=abs(`k')'"
    if `k'+16.5>0 {
        if "`slopeexpr'"=="" local slopeexpr "`a'*bd_ev_`s'"
        else local slopeexpr "`slopeexpr' + `a'*bd_ev_`s'"
    }
    else local slopeexpr "`slopeexpr' - `a'*bd_ev_`s'"
}
* Two pre-event sub-windows. The weights were fitted on 2022-23 ONLY, so
* k=-32..-9 is in-sample and k=-8..-2 (Jan-Jul 2024) is a genuine holdout.
* Comparing coefficient amplitude across the two is the out-of-sample test.
local insamp ""
local holdout ""
local holdtest ""
forvalues k=-32/-2 {
    local s="m`=abs(`k')'"
    if `k'<=-9 local insamp "`insamp' `s'"
    else {
        local holdout "`holdout' `s'"
        local holdtest "`holdtest' (bd_ev_`s' - `premean' = 0)"
    }
}

foreach pf of local prefixes {
    local premean_`pf' "(`presum_`pf'')/32"
}

* Weights and sample keyed on the BACKDATED vector.
tempvar mw mins
gen double `mw'=cond(treated_nc==1,1/39,mc_w_joint_bd)
gen byte `mins'=treated_nc==1 | (mc_w_joint_bd>1e-10 & mc_donor==1)
assert `mw'>0 if `mins'

tempname P Q
postfile `Q' str12 outcome double drift_p holdout_p ///
    max_abs_insample rms_insample max_abs_holdout rms_holdout ///
    long N_clust N_obs using "${MC_OUT}/data/04d_backdated_pretrends.dta", replace
postfile `P' str12 spec str12 outcome str16 term str12 window byte months ///
    double beta se lower95 upper95 pvalue long N_clust N_obs ///
    using "${MC_OUT}/data/04d_backdated_period_effects.dta", replace

foreach y in mc_yemp mc_ylf mc_yjobs mc_yurate mc_ydiff mc_yshare {

    * ---- aggregate ----
    di as txt "04d aggregate: `y'"
    quietly areg `y' `events' i.mdate if `mins' [aw=`mw'], ///
        absorb(county_id) vce(cluster county_id)
    local nev : word count `events'
    local spec_tag  "04d backdated / `y' (aggregate)"
    local spec_cmd  "areg `y' <`nev' event indicators> i.mdate if <sample> [aw=<backdated weight>], absorb(county_id) vce(cluster county_id)"
    local spec_rhs  "`nev' treated-by-event-month indicators, as in stage 04"
    local spec_sample "identical to stage 04, but donors are those with mc_w_joint_bd>1e-10"
    local spec_wt   "treated=1/39; donor=mc_w_joint_bd, fitted on 2022-23 ONLY (never saw Jan-Aug 2024)"
    local spec_ref  "recentred to the 32-month pre-event mean, as in stage 04"
    local spec_out  "04d_backdated_period_effects.dta -> backdated_weights_robustness.csv"
    include "${MC_ROOT}/steps/_logspec.doh"
    local df=e(df_r)
    local tcrit=invttail(`df',.025)
    local g=e(N_clust)
    local n=e(N)

    * ---- out-of-sample pre-event diagnostics ----
    * No differential TREND is what parallel trends requires; amplitude is
    * reported beside it because the two can move independently.
    quietly lincom `slopeexpr'
    local dp=r(p)
    quietly test `holdtest'
    local hp=r(p)
    foreach part in insamp holdout {
        local mx=0
        local ss=0
        local cnt=0
        foreach s of local `part' {
            quietly lincom bd_ev_`s' - `premean'
            local bb=r(estimate)
            if abs(`bb')>`mx' local mx=abs(`bb')
            local ss=`ss'+(`bb')^2
            local cnt=`cnt'+1
        }
        local max_`part'=`mx'
        local rms_`part'=sqrt(`ss'/`cnt')
    }
    post `Q' ("`y'") (`dp') (`hp') (`max_insamp') (`rms_insamp') ///
        (`max_holdout') (`rms_holdout') (`g') (`n')

    foreach window in immediate year2025 late {
        if "`window'"=="immediate" local ks "1 2 3"
        if "`window'"=="year2025"  local ks "4 5 6 7 8 9 10 11 12"
        if "`window'"=="late"      local ks "14 15 16 17 18"
        local expr ""
        local len=0
        foreach k of local ks {
            local len=`len'+1
            if "`expr'"=="" local expr "bd_ev_p`k' - `premean'"
            else local expr "`expr' + bd_ev_p`k' - `premean'"
        }
        quietly lincom (`expr')/`len'
        post `P' ("aggregate") ("`y'") ("all") ("`window'") (`len') ///
            (r(estimate)) (r(se)) (r(estimate)-`tcrit'*r(se)) ///
            (r(estimate)+`tcrit'*r(se)) (r(p)) (`g') (`n')
    }

    * ---- by HDSI tercile, one regression ----
    di as txt "04d terciles: `y'"
    quietly areg `y' `all_bt1' `all_bt2' `all_bt3' i.mdate if `mins' [aw=`mw'], ///
        absorb(county_id) vce(cluster county_id)
    local nev2 : word count `all_bt1' `all_bt2' `all_bt3'
    local spec_tag  "04d backdated / `y' (hdsi terciles)"
    local spec_cmd  "areg `y' <`nev2' tercile-by-event-month indicators> i.mdate if <sample> [aw=<backdated weight>], absorb(county_id) vce(cluster county_id)"
    local spec_rhs  "`nev2' regressors: stage-04d indicators split three ways by HDSI tercile"
    local spec_sample "identical to the aggregate model above"
    local spec_wt   "treated=1/39; donor=mc_w_joint_bd, fitted on 2022-23 ONLY"
    local spec_ref  "each tercile recentred to its own 32-month pre-event mean; high-minus-low with joint covariance"
    local spec_out  "04d_backdated_period_effects.dta -> backdated_weights_robustness.csv"
    include "${MC_ROOT}/steps/_logspec.doh"
    local df=e(df_r)
    local tcrit=invttail(`df',.025)
    local g=e(N_clust)
    local n=e(N)
    local labels "low mid high"
    foreach window in immediate year2025 late {
        if "`window'"=="immediate" local ks "1 2 3"
        if "`window'"=="year2025"  local ks "4 5 6 7 8 9 10 11 12"
        if "`window'"=="late"      local ks "14 15 16 17 18"
        forvalues j=1/3 {
            local lab : word `j' of `labels'
            local expr ""
            local len=0
            foreach k of local ks {
                local len=`len'+1
                if "`expr'"=="" local expr "bt`j'_p`k' - `premean_bt`j''"
                else local expr "`expr' + bt`j'_p`k' - `premean_bt`j''"
            }
            quietly lincom (`expr')/`len'
            post `P' ("tercile") ("`y'") ("`lab'") ("`window'") (`len') ///
                (r(estimate)) (r(se)) (r(estimate)-`tcrit'*r(se)) ///
                (r(estimate)+`tcrit'*r(se)) (r(p)) (`g') (`n')
        }
        * high minus low, taken inside the single regression
        local ehi ""
        local elo ""
        local len=0
        foreach k of local ks {
            local len=`len'+1
            if "`ehi'"=="" local ehi "bt3_p`k' - `premean_bt3'"
            else local ehi "`ehi' + bt3_p`k' - `premean_bt3'"
            if "`elo'"=="" local elo "bt1_p`k' - `premean_bt1'"
            else local elo "`elo' + bt1_p`k' - `premean_bt1'"
        }
        quietly lincom ((`ehi')-(`elo'))/`len'
        post `P' ("tercile") ("`y'") ("high_minus_low") ("`window'") (`len') ///
            (r(estimate)) (r(se)) (r(estimate)-`tcrit'*r(se)) ///
            (r(estimate)+`tcrit'*r(se)) (r(p)) (`g') (`n')
    }
}
postclose `P'
postclose `Q'

use "${MC_OUT}/data/04d_backdated_pretrends.dta", clear
assert !missing(drift_p,holdout_p)
isid outcome
gen double f=cond(outcome=="mc_yurate",1,100)
foreach v in max_abs_insample rms_insample max_abs_holdout rms_holdout {
    replace `v'=`v'*f
}
drop f
label var max_abs_holdout "Largest pre-event coefficient in the Jan-Jul 2024 holdout"
export delimited using "${MC_OUT}/tables/backdated_pretrends.csv", replace
di as txt "04d: out-of-sample pre-event diagnostics (weights never saw 2024H1)"
list outcome drift_p holdout_p rms_insample rms_holdout, noobs

use "${MC_OUT}/data/04d_backdated_period_effects.dta", clear
assert !missing(beta,se,pvalue)
isid spec outcome term window
gen double display_beta=cond(outcome=="mc_yurate",beta,100*beta)
gen double display_lo=cond(outcome=="mc_yurate",lower95,100*lower95)
gen double display_hi=cond(outcome=="mc_yurate",upper95,100*upper95)
gen str12 display_units=cond(outcome=="mc_yurate","pp","log-%")
export delimited using "${MC_OUT}/tables/backdated_weights_robustness.csv", replace
di as result "04d: robustness estimates on backdated weights written."
list outcome term window display_beta pvalue if spec=="tercile" & ///
    window=="immediate" & inlist(term,"low","high","high_minus_low"), noobs sepby(outcome)
