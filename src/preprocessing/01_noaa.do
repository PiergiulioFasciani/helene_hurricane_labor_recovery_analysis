/* County-average NOAA precipitation, Sep 25-30 2024, 11 states.
   Primary rainfall exposure is meteorological, independent of FEMA. */
version 16.0
di as txt "01: NOAA precipitation"
local csv "${RAW}/prcp-202409-cty-scaled.csv"
capture confirm file "`csv'"
if _rc {
    local url "https://www.ncei.noaa.gov/data/nclimgrid-daily/access/averages/2024/prcp-202409-cty-scaled.csv"
    capture noisily copy "`url'" "`csv'", replace
    if _rc {
        shell /usr/bin/curl -fL --retry 3 -o "`csv'" "`url'"
    }
}
capture confirm file "`csv'"
if _rc {
    di as error "NOAA download missing: `csv'"
    exit 601
}
import delimited using "`csv'", clear varnames(nonames) stringcols(2)
foreach v in v2 v3 v31 v32 v33 v34 v35 v36 {
    capture confirm variable `v'
    if _rc {
        di as error "NOAA input missing expected column `v'; inspect CSV."
        exit 459
    }
}
rename v2 county_code
rename v3 noaa_county_name
replace county_code = ustrregexra(ustrtrim(county_code), "[^0-9]", "")
assert strlen(county_code)==5 & !missing(real(county_code))
gen str2 ncei_state=substr(county_code,1,2)
gen str2 state=""
gen str2 census_state=""
foreach pair in "01 AL 01" "08 FL 12" "09 GA 13" "12 IN 18" ///
 "15 KY 21" "31 NC 37" "33 OH 39" "38 SC 45" "40 TN 47" ///
 "44 VA 51" "46 WV 54" {
    tokenize "`pair'"
    replace state="`2'" if ncei_state=="`1'"
    replace census_state="`3'" if ncei_state=="`1'"
}
keep if state!=""
* NOTE: NOAA geographic codes use a different STATE prefix. County suffix
* is a candidate Census code; validated against Census code+NAME in step 04.
gen str5 fips=census_state+substr(county_code,3,3)
forvalues d=25/30 {
    local c=`d'+6
    rename v`c' prcp_`d'
    capture confirm numeric variable prcp_`d'
    if _rc destring prcp_`d', replace force
    replace prcp_`d'=. if prcp_`d'<=-900
}
egen byte n_rain_days=rownonmiss(prcp_25-prcp_30)
egen double prcp_total_mm=rowtotal(prcp_25-prcp_30)
replace prcp_total_mm=. if n_rain_days!=6
gen double prcp_avg_mm=prcp_total_mm/6
keep fips state census_state county_code noaa_county_name ///
     prcp_total_mm prcp_avg_mm n_rain_days
isid fips
count
di as result "NOAA study region county count: " r(N)
save "${INT}/noaa_11states.dta", replace
export delimited using "${OUT}/audits/noaa_county_fips_candidates.csv", replace
