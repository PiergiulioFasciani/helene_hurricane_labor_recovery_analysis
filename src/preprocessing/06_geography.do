/* Frozen 39-county treatment + geographically isolated donor universe.
   NO synthetic weights or outcome estimates are calculated in this stage.
   Donor completeness uses only the 24 training months, never 2024 validation. */
version 16.0
di as txt "06: Frozen treated roster and four-outcome donor geography"
local training_start=ym(2022,1)
local training_end=ym(2023,12)
local validation_end=ym(2024,8)
local helene_month=ym(2024,9)
local end_month=ym(2026,9)

* National panel MUST remain the master of the merge.
use "${INT}/labor_national_2022_2026m09.dta", clear
local N_before=_N
merge m:1 fips using "${INT}/helene_exposure_11states.dta", ///
    keep(master match) gen(exposure_merge)
assert _N==`N_before'
gen byte in_noaa_sample=(exposure_merge==3)
quietly count if in_noaa_sample & mdate==ym(2024,9)
di as result "Noaa county matches: " r(N)
assert in_noaa_sample==1 if !missing(prcp_total_mm)

* Freeze the original treated county roster (39 FIPS identifiers) to guard
* against later FEMA revisions or post-outcome reclassification.
tempfile frozen_treated
preserve
    import delimited using "${REF}/treated_NC_county_roster_2026-11-13.csv", ///
        clear varnames(1) stringcols(1)
    keep fips
    isid fips
    count
    assert r(N)==39
    gen byte treated_nc=1
    save `frozen_treated'
restore
merge m:1 fips using `frozen_treated', keep(master match) nogen
replace treated_nc=0 if missing(treated_nc)
assert substr(fips,1,2)=="37" & in_noaa_sample==1 if treated_nc==1
quietly count if mdate==ym(2024,9) & substr(fips,1,2)=="37" & ///
    treated_nc!=fema_positive_nc
if r(N)>0 {
    di as error "FEMA-positive NC roster differs from frozen 39-county treatment (n=" r(N) ")."
    if ${STOP_IF_TREATMENT_CHANGED} exit 459
}

* Donor restriction uses COUNTY STATE FIPS; geographic isolation cannot
* prove absence of commuting or other inter-county interference.
gen str2 state_fips=substr(fips,1,2)
gen byte impact_state= ///
    inlist(state_fips,"01","12","13","17","18","21","37","39","45") | ///
    inlist(state_fips,"47","51","54")
gen byte buffer_state= ///
    inlist(state_fips,"05","11","19","24","26","28","29","42","55")
gen byte other_2024_hurricane_state=inlist(state_fips,"22","48")
gen byte eligible_donor=!impact_state & !buffer_state & !other_2024_hurricane_state
assert eligible_donor==0 if treated_nc==1
preserve
    keep if mdate==`training_end'
    keep fips state_fips treated_nc eligible_donor impact_state ///
        buffer_state other_2024_hurricane_state
    isid fips
    export delimited using "${OUT}/audits/county_geographic_eligibility.csv", replace
restore

* Save complete unfiltered national+exposure dataset before SCM.
isid fips mdate


* Four-outcome training eligibility. The untouched Jan-Aug24 holdout is NOT a filter.
gen byte joint_training_ok = inrange(mdate,ym(2022,1),ym(2023,12)) & ///
    laus_employed>0 & !missing(laus_employed) & ///
    laus_labor_force>0 & !missing(laus_labor_force) & ///
    qcew_jobs>0 & !missing(qcew_jobs) & ///
    !missing(laus_unemployment_rate)
bysort fips: egen int joint_training_n=total(joint_training_ok)
gen byte joint_eligible_donor=(eligible_donor==1 & treated_nc==0 & joint_training_n==24)
assert joint_training_n==24 if treated_nc==1
preserve
    keep if mdate==ym(2024,8)
    count if treated_nc==1
    assert r(N)==39
    count if joint_eligible_donor==1
    assert r(N)>=500
    di as result "Four-outcome training-complete eligible donors: " r(N)
    keep fips treated_nc state_fips eligible_donor joint_training_n joint_eligible_donor
    isid fips
    export delimited using "${OUT}/audits/joint_eligible_donor_geography.csv", replace
restore
isid fips mdate
save "${INT}/national_labor_exposure_prefit.dta", replace
