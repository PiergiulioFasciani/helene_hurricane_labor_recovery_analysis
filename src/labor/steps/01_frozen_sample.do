/* Freeze treatment; require 2022-23 availability of all four objectives.
   Never use future observations to determine donor eligibility. */
version 16.0
di as txt "01: Freeze geography and build common-outcome candidate panel"
use "${MC_INPUT}", clear
foreach v in fips county_id mdate treated_nc es_donor_pool ///
    laus_employed laus_labor_force ///
    laus_unemployed laus_unemployment_rate qcew_jobs {
    capture confirm variable `v'
    if _rc {
        di as error "Missing source variable: `v'"
        exit 111
    }
}
isid fips mdate
assert strlen(fips)==5
assert !treated_nc | substr(fips,1,2)=="37"
assert !(treated_nc & es_donor_pool==0)
assert _N==3221*57

* Verify the exact 39 original FIPS, not merely their total count.
tempfile mc_roster
preserve
    import delimited using "${H_ROOT}/references/labor/frozen_treated_nc_39.csv", ///
        clear varnames(1) stringcols(1)
    keep fips
    isid fips
    quietly count
    assert r(N)==39
    gen byte mc_fixed_treated=1
    save `mc_roster'
restore
merge m:1 fips using `mc_roster', keep(master match) nogen
replace mc_fixed_treated=0 if missing(mc_fixed_treated)
assert treated_nc==mc_fixed_treated
drop mc_fixed_treated

* Keep only frozen treated units and strict-geography, pre-training-complete controls.
keep if treated_nc==1 | es_donor_pool==1
keep fips county_id mdate treated_nc es_donor_pool ///
    laus_employed laus_unemployed laus_labor_force laus_unemployment_rate qcew_jobs
sort fips mdate
by fips: assert _N==57
by fips: gen byte mc_step=_n
assert mdate==ym(2022,1)+mc_step-1

foreach v in laus_employed laus_labor_force qcew_jobs {
    assert `v'>=0 if !missing(`v')
}
gen double mc_raw_emp=100*ln(laus_employed) if laus_employed>0
gen double mc_raw_lf=100*ln(laus_labor_force) if laus_labor_force>0
gen double mc_raw_jobs=100*ln(qcew_jobs) if qcew_jobs>0
gen double mc_raw_ur=laus_unemployment_rate
gen double mc_raw_unemp=ln(laus_unemployed) if laus_unemployed>0

foreach o in emp lf jobs ur {
    bysort fips: egen double mc_base_`o'= ///
        mean(cond(inrange(mdate,ym(2022,1),ym(2022,12)),mc_raw_`o',.))
    gen double mc_`o'=mc_raw_`o'-mc_base_`o' ///
        if !missing(mc_raw_`o',mc_base_`o')
}

* Exactly 24 observations in 2022-23 for ALL FOUR pre-treatment outcomes.
gen byte mc_goodpre=inrange(mdate,ym(2022,1),ym(2023,12)) & ///
    !missing(mc_emp,mc_lf,mc_jobs,mc_ur)
bysort fips: egen byte mc_npre=total(mc_goodpre)
gen byte mc_precomplete=mc_npre==24
assert mc_precomplete if treated_nc==1

* Joint donors are screened solely for 2022-23 coverage.
gen byte mc_donor=(treated_nc==0 & es_donor_pool==1 & mc_precomplete==1)
preserve
    keep if mdate==ym(2024,8)
    quietly count if treated_nc==1
    assert r(N)==39
    quietly count if mc_donor==1
    assert r(N)>=500
    di as result "Prespecified treated counties: 39; eligible all-outcome donors: " r(N)
    keep fips treated_nc es_donor_pool mc_precomplete mc_npre mc_donor
    export delimited using "${MC_OUT}/tables/frozen_geography_and_precoverage.csv", replace
restore

keep if treated_nc==1 | mc_donor==1
sort fips mdate
by fips: gen byte mc_step_check=_n
assert mc_step==mc_step_check
drop mc_step_check
by fips: assert _N==57
by fips: assert treated_nc==treated_nc[1]
by fips: assert mc_donor==mc_donor[1]
egen int mc_uid=group(fips)
sort mc_uid mdate
isid mc_uid mdate

* QCEW and resident LAUS must align in the primary analysis window.
* We check later coverage but do not SELECT donors from it.
assert !missing(mc_emp,mc_lf,mc_jobs,mc_ur) ///
    if mdate<=ym(2026,3) & mdate!=ym(2025,10)
assert missing(mc_emp,mc_lf,mc_ur) if mdate==ym(2025,10)

label var mc_donor "Frozen strict geographic donor, complete all outcomes 2022-23"
label var mc_emp "100 x log employed residents, centered on own 2022 mean"
label var mc_lf "100 x log resident labor force, centered on own 2022 mean"
label var mc_jobs "100 x log county QCEW jobs, centered on own 2022 mean"
label var mc_ur "Unemployment rate (pp), centered on own 2022 mean"
compress
save "${MC_OUT}/data/01_frozen_candidate_panel.dta", replace
