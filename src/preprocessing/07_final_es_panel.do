/* National county-month panel, built entirely from raw-record upstream stages.
   Strict control universe defined by geography + 2022-23 four-outcome coverage.
   No original employment SCM or synthetic weights are built or imported. */
version 16.0
di as txt "07: Generate final NATIONAL four-outcome research panel"
use "${INT}/national_labor_exposure_prefit.dta", clear
isid fips mdate
assert _N==3221*57
assert strlen(fips)==5
bysort fips: assert _N==57
quietly count if mdate==ym(2024,8) & treated_nc==1
assert r(N)==39
assert joint_training_n==24 if treated_nc==1
assert !(treated_nc & eligible_donor)

gen byte es_donor_pool=(treated_nc==1 | joint_eligible_donor==1)
label var es_donor_pool "39 treated + geographically isolated four-outcome 2022-23 complete donors"
label var joint_eligible_donor "Geographic eligible, all four pre-training outcomes present 2022-23"
label var treated_nc "Frozen 39 FEMA-positive NC counties (post-hoc exposure-defined)"
label var hdsi "Descriptive exposure score, not exogenous damage"

* Event-time variables retained for downstream analysis; no estimations here.
gen int event_month=mdate-ym(2024,9)
gen byte transition_sep2024=(event_month==0)
gen byte post_helene_full=(event_month>=1)
assert laus_oct2025_unavailable==(mdate==ym(2025,10))
assert missing(laus_employed) & missing(laus_unemployed) & ///
    missing(laus_labor_force) & missing(laus_unemployment_rate) ///
    if mdate==ym(2025,10)
format mdate %tm
xtset county_id mdate
order fips county_id mdate year month treated_nc eligible_donor ///
    joint_eligible_donor es_donor_pool joint_training_n event_month ///
    laus_employed laus_unemployed laus_labor_force laus_unemployment_rate ///
    qcew_jobs prcp_total_mm total_damage_usd hdsi
compress
save "${OUT}/panels/Helene_NC_ES_preinference_panel.dta", replace

* Audits contain actual source-built coverage and full donor assignments.
preserve
    collapse (count) n_laus=laus_employed n_qcew=qcew_jobs ///
        (sum) n_pool=es_donor_pool n_source_oct2025_values=laus_oct2025_source_flag, by(mdate)
    format mdate %tm
    export delimited using "${OUT}/audits/final_panel_coverage_by_month.csv", replace
restore
preserve
    keep if mdate==ym(2024,8)
    keep fips treated_nc eligible_donor joint_eligible_donor es_donor_pool ///
        joint_training_n in_noaa_sample has_fema_record prcp_total_mm hdsi
    isid fips
    export delimited using "${OUT}/audits/county_treatment_donor_assignment.csv", replace
restore
di as result "SOURCE-BUILT 3221 x 57 panel saved; NO old control weights computed."
