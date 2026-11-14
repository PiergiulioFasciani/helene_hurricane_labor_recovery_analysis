/* Preserve ALL NOAA rain counties, even without FEMA. Rebuild HDSI as
   a descriptive matched-sample variable; NC treatment frozen separately. */
version 16.0
di as txt "04: Assemble NOAA, Census and FEMA"
use "${INT}/noaa_11states.dta", clear
merge 1:1 fips using "${INT}/census_pep_2024.dta", keep(master match)
tab _merge
quietly count if _merge==1
if r(N)>0 {
    preserve
        keep if _merge==1
        export delimited using "${OUT}/audits/noaa_unmatched_pep.csv", replace
    restore
    di as error "NOAA FIPS not matched to Census PEP; review audit; stopping."
    exit 459
}
drop _merge
* Cross-check geographic names; never silently assign a county by FIPS alone.
gen str100 noaa_key=lower(itrim(strtrim(noaa_county_name)))
replace noaa_key=regexr(noaa_key,"^[a-z][a-z]: *", "")
gen str100 census_key=lower(itrim(strtrim(census_county_name)))
foreach nm in noaa_key census_key {
    replace `nm'=subinstr(`nm',".","",.)
    replace `nm'=subinstr(`nm',"'","",.)
    replace `nm'=regexr(`nm'," county$","")
}
gen byte name_match=(noaa_key==census_key)
quietly count if name_match==0
if r(N)>0 {
    preserve
        keep if name_match==0
        keep fips state noaa_county_name census_county_name
        export delimited using "${OUT}/audits/noaa_census_name_disagreements.csv", replace
    restore
    di as error "NOAA/Census county names disagree. Review exported audit (spelling could differ)."
}
else {
    capture erase "${OUT}/audits/noaa_census_name_disagreements.csv"
}
drop noaa_key census_key
merge 1:1 fips using "${INT}/census_acs_2023.dta", keep(master match)
gen byte has_acs=(_merge==3)
drop _merge
merge 1:1 fips using "${INT}/fema_helene_county.dta", keep(master match)
gen byte has_fema_record=(_merge==3)
drop _merge
* Missing verified losses remain missing; no FEMA record != zero loss.
gen double damage_per_resident_usd=total_damage_usd/population_2024 ///
    if population_2024>0 & !missing(total_damage_usd)
gen double housing_value_proxy_usd=households_2023*median_home_value_2023 ///
    if households_2023>0 & median_home_value_2023>0
gen double econ_burden_pct=100*total_damage_usd/housing_value_proxy_usd ///
    if housing_value_proxy_usd>0 & !missing(total_damage_usd)
gen byte score_sample=has_fema_record & !missing(prcp_total_mm,econ_burden_pct)
quietly count if score_sample
local n=r(N)
if `n'<2 {
    di as error "Insufficient counties for the complete-case HDSI."
    exit 459
}
egen double rain_rank=rank(prcp_total_mm) if score_sample
egen double economic_rank=rank(econ_burden_pct) if score_sample
gen double rainfall_score=100*(rain_rank-1)/(`n'-1) if score_sample
gen double economic_damage_score=100*(economic_rank-1)/(`n'-1) if score_sample
gen double hdsi=.4*rainfall_score+.6*economic_damage_score if score_sample
egen double rain_rank_all=rank(prcp_total_mm)
quietly count if !missing(prcp_total_mm)
local nr=r(N)
gen double noaa_exposure_score=100*(rain_rank_all-1)/(`nr'-1) ///
    if !missing(prcp_total_mm) & `nr'>1
* Do not define treatment from a post-disaster index for the estimand.
gen byte fema_positive_nc=substr(fips,1,2)=="37" & ///
    !missing(total_damage_usd) & total_damage_usd>0
isid fips
save "${INT}/helene_exposure_11states.dta", replace
export delimited using "${OUT}/audits/helene_exposure_11states.csv", replace
