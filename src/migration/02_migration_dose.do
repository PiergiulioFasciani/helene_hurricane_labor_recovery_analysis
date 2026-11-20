/* 02: ANNUAL MIGRATION AND POPULATION BY DAMAGE TERCILE (added 18 Nov 2026).
   Migration enters the design as an OUTCOME, never as a control: it is
   realized after treatment, so conditioning on it would be a bad control and
   would remove the very channel under study.

   Same estimand shape as the labor event study, at annual July-to-July
   frequency: county and year effects, treated-by-year indicators, 2024 the
   omitted (last unexposed) year, standard errors clustered by county.
   HDSI terciles are read from the labor stage so the groups are IDENTICAL to
   those in src/labor/steps/04c_dose_event_study.do.

   Two comparison schemes are estimated and both are reported:
     unweighted : every geographically eligible county counts equally
                  (the scheme used by stage 01).
     synthetic  : donors carry the fitted labor synthetic-control weights,
                  treated counties carry 1/39, exactly as in stage 04/04c.
   The synthetic scheme makes migration comparable with the labor outcomes;
   the unweighted scheme is retained because the weights were fitted to labor
   dynamics, not to population. Choose on the reported pre-tests, not on the
   2025 estimate. */
version 16.0
di as txt "02: annual migration and population by HDSI tercile"

local panel "${H_ROOT}/outputs/migration/census_county_year.dta"
capture confirm file "`panel'"
if _rc {
    di as error "Missing county-year panel. Run src/migration/01_census_v2025.do first."
    exit 601
}
local terc "${H_ROOT}/outputs/labor/tables/dose_treated_exposure.csv"
capture confirm file "`terc'"
if _rc {
    di as error "Missing tercile assignment. Run src/labor/steps/04c_dose_event_study.do first."
    exit 601
}
local wts "${H_ROOT}/outputs/labor/tables/donor_weights.csv"
capture confirm file "`wts'"
if _rc {
    di as error "Missing donor weights. Run src/labor/run_multicontrol.do first."
    exit 601
}

tempfile assign terciles weights
use "${H_ROOT}/data/input/Helene_NC_ES_preinference_panel.dta", clear
keep if mdate==ym(2024,8)
keep fips treated_nc es_donor_pool
isid fips
save `assign'

import delimited using "`terc'", clear varnames(1) stringcols(1)
keep fips md_ter
isid fips
quietly count
assert r(N)==39
save `terciles'

import delimited using "`wts'", clear varnames(1) stringcols(1)
keep fips mc_w_joint
isid fips
save `weights'

use "`panel'", clear
isid fips year
merge m:1 fips using `assign', keep(match) nogen
merge m:1 fips using `terciles', keep(master match) nogen
merge m:1 fips using `weights', keep(master match) nogen
replace md_ter=0 if treated_nc==0
replace mc_w_joint=0 if missing(mc_w_joint)
assert !missing(md_ter) if treated_nc==1
assert inlist(md_ter,1,2,3) if treated_nc==1

keep if (treated_nc==1 | es_donor_pool==1) & inrange(year,2021,2025)
egen long demo_id=group(fips)

* Weighting schemes. Treated counties always carry 1/39, as in stage 04.
gen double dm_w_unwt=1
gen double dm_w_syn=cond(treated_nc==1,1/39,mc_w_joint)
gen byte dm_in_unwt=1
gen byte dm_in_syn=(treated_nc==1 | (mc_w_joint>1e-10 & treated_nc==0))

* Treated-by-year indicators: pooled and one set per damage tercile.
foreach yr in 2021 2022 2023 2025 {
    gen byte dm_all_`yr'=(treated_nc==1 & year==`yr')
    forvalues g=1/3 {
        gen byte dm_t`g'_`yr'=(treated_nc==1 & md_ter==`g' & year==`yr')
    }
}

tempname R
postfile `R' str10 scheme str18 outcome str6 group double b2025 se lower95 ///
    upper95 pvalue pretest_p long n_treated N_clust N_obs ///
    using "${H_ROOT}/outputs/migration/02_migration_dose.dta", replace

foreach scheme in unwt syn {
    foreach y in netdom_per1000 pop_growth_pct {
        * ---- pooled: all 39 treated counties ----
        capture quietly areg `y' dm_all_2021 dm_all_2022 dm_all_2023 dm_all_2025 ///
            i.year if dm_in_`scheme'==1 [aw=dm_w_`scheme'], ///
            absorb(demo_id) vce(cluster demo_id)
        if !_rc {
            local df=e(df_r)
            local crit=invttail(`df',.025)
            local b=_b[dm_all_2025]
            local se=_se[dm_all_2025]
            local pv=2*ttail(`df',abs(`b'/`se'))
            local pp=.
            capture quietly test dm_all_2021 dm_all_2022 dm_all_2023
            if !_rc local pp=r(p)
            local _c=e(N_clust)
            local _n=e(N)
            tempname _ls
            file open `_ls' using "${H_ROOT}/outputs/audit/REGRESSIONS.txt", write append text
            file write `_ls' "MODEL  02 migration / `scheme' / `y' / pooled (all 39)" _n
            file write `_ls' "  command    areg `y' dm_all_2021 dm_all_2022 dm_all_2023 dm_all_2025 i.year if dm_in_`scheme'==1 [aw=dm_w_`scheme'], absorb(demo_id) vce(cluster demo_id)" _n
            file write `_ls' "  regressors one treated-by-year indicator per year; 2021-2023 are the pre-tests, 2025 is the effect" _n
            file write `_ls' "  absorbed   demo_id (county fixed effects, annual panel)" _n
            file write `_ls' "  sample     annual July-to-July periods 2021-2025; 2024 is the omitted last unexposed year" _n
            file write `_ls' "  weights    scheme `scheme' (unwt = eligible counties equally; syn = donors carrying the labor synthetic weights)" _n
            file write `_ls' "  vce        cluster(demo_id)" _n
            file write `_ls' "  fitted     N=" %12.0fc (`_n') "  clusters=" %6.0f (`_c') "  df_r=" %6.0f (`df') _n
            file write `_ls' "  reference  2024; the reported effect is the 2025 coefficient" _n
            file write `_ls' "  writes     02_migration_dose.dta -> migration_by_damage.csv" _n
            file write `_ls' "" _n
            file close `_ls'
            post `R' ("`scheme'") ("`y'") ("all") (`b') (`se') ///
                (`b'-`crit'*`se') (`b'+`crit'*`se') (`pv') (`pp') ///
                (39) (e(N_clust)) (e(N))
        }
        * ---- by damage tercile, one regression, joint covariance ----
        capture quietly areg `y' dm_t1_2021 dm_t1_2022 dm_t1_2023 dm_t1_2025 ///
            dm_t2_2021 dm_t2_2022 dm_t2_2023 dm_t2_2025 ///
            dm_t3_2021 dm_t3_2022 dm_t3_2023 dm_t3_2025 ///
            i.year if dm_in_`scheme'==1 [aw=dm_w_`scheme'], ///
            absorb(demo_id) vce(cluster demo_id)
        if !_rc {
            local df=e(df_r)
            local crit=invttail(`df',.025)
            local _c=e(N_clust)
            local _n=e(N)
            tempname _ls
            file open `_ls' using "${H_ROOT}/outputs/audit/REGRESSIONS.txt", write append text
            file write `_ls' "MODEL  02 migration / `scheme' / `y' / by damage tercile" _n
            file write `_ls' "  command    areg `y' dm_t{1,2,3}_{2021,2022,2023,2025} i.year if dm_in_`scheme'==1 [aw=dm_w_`scheme'], absorb(demo_id) vce(cluster demo_id)" _n
            file write `_ls' "  regressors the same treated-by-year indicators split three ways by HDSI tercile, one regression so the high-minus-low contrast uses the joint covariance" _n
            file write `_ls' "  absorbed   demo_id (county fixed effects, annual panel)" _n
            file write `_ls' "  sample     annual July-to-July periods 2021-2025; 2024 is the omitted last unexposed year" _n
            file write `_ls' "  weights    scheme `scheme' (unwt = eligible counties equally; syn = donors carrying the labor synthetic weights)" _n
            file write `_ls' "  vce        cluster(demo_id)" _n
            file write `_ls' "  fitted     N=" %12.0fc (`_n') "  clusters=" %6.0f (`_c') "  df_r=" %6.0f (`df') _n
            file write `_ls' "  reference  2024; the reported effect is the 2025 coefficient" _n
            file write `_ls' "  writes     02_migration_dose.dta -> migration_by_damage.csv" _n
            file write `_ls' "" _n
            file close `_ls'
            forvalues g=1/3 {
                local b=_b[dm_t`g'_2025]
                local se=_se[dm_t`g'_2025]
                local pv=2*ttail(`df',abs(`b'/`se'))
                local pp=.
                capture quietly test dm_t`g'_2021 dm_t`g'_2022 dm_t`g'_2023
                if !_rc local pp=r(p)
                local gname="low"
                if `g'==2 local gname="mid"
                if `g'==3 local gname="high"
                quietly count if treated_nc==1 & md_ter==`g' & year==2025
                post `R' ("`scheme'") ("`y'") ("`gname'") (`b') (`se') ///
                    (`b'-`crit'*`se') (`b'+`crit'*`se') (`pv') (`pp') ///
                    (r(N)) (e(N_clust)) (e(N))
            }
            * high minus low, from the joint covariance
            quietly lincom dm_t3_2025 - dm_t1_2025
            local b=r(estimate)
            local se=r(se)
            post `R' ("`scheme'") ("`y'") ("hi-lo") (`b') (`se') ///
                (`b'-`crit'*`se') (`b'+`crit'*`se') (r(p)) (.) ///
                (26) (e(N_clust)) (e(N))
        }
    }
}
postclose `R'

use "${H_ROOT}/outputs/migration/02_migration_dose.dta", clear
isid scheme outcome group
export delimited using "${H_ROOT}/outputs/migration/migration_by_damage.csv", replace
list scheme outcome group b2025 lower95 upper95 pvalue pretest_p, noobs sepby(scheme outcome)
di as result "02 complete: migration and population by damage tercile (outcome, not control)."
