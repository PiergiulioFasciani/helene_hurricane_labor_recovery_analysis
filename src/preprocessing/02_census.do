/* Census 2024 county population + 2023 ACS pre-storm covariates. */
version 16.0
di as txt "02: Census PEP and ACS"
local pep "${RAW}/co-est2024-alldata.csv"
capture confirm file "`pep'"
if _rc {
    local url "https://www2.census.gov/programs-surveys/popest/datasets/2020-2024/counties/totals/co-est2024-alldata.csv"
    capture noisily copy "`url'" "`pep'", replace
    if _rc shell /usr/bin/curl -fL --retry 3 -o "`pep'" "`url'"
}
capture confirm file "`pep'"
if _rc exit 601
import delimited using "`pep'", clear varnames(1) stringcols(4 5)
keep if sumlev==50
capture confirm numeric variable state
if !_rc tostring state, replace format(%02.0f)
capture confirm numeric variable county
if !_rc tostring county, replace format(%03.0f)
replace state=strtrim(state)
replace county=strtrim(county)
keep if strpos("|01|12|13|18|21|37|39|45|47|51|54|", ///
    "|"+state+"|")>0
gen str5 fips=state+county
rename ctyname census_county_name
rename popestimate2024 population_2024
keep fips census_county_name population_2024
assert strlen(fips)==5 & !missing(real(fips))
isid fips
save "${INT}/census_pep_2024.dta", replace

* Use a CSV returned by Census API. Never silently import an error page.
local acs "${RAW}/ACS_2023_5year_counties.csv"
capture confirm file "`acs'"
if _rc {
    if "${CENSUS_KEY}"=="" | "${CENSUS_KEY}"=="INSERT_YOUR_CENSUS_KEY" {
        di as error "Need a Census API key in untracked config.local.do (or valid local ACS CSV)."
        exit 198
    }
    local url "https://api.census.gov/data/2023/acs/acs5?get=NAME,B25003_001E,B25077_001E,B19013_001E&for=county%3A%2A&outputFormat=csv&key=${CENSUS_KEY}"
    shell /usr/bin/curl -fL --retry 3 --silent --show-error --max-time 240 -o "`acs'" "`url'"
}
capture confirm file "`acs'"
if _rc exit 601
file open acs_header using "`acs'", read text
file read acs_header first_line
file close acs_header
if strpos(upper(`"`first_line'"'),"B25003_001E")==0 | ///
   strpos(upper(`"`first_line'"'),"B25077_001E")==0 {
    di as error "ACS CSV lacks expected header (possibly API error). Delete invalid file and rerun."
    exit 459
}
import delimited using "`acs'", clear varnames(1) stringcols(_all) bindquote(strict)
foreach v in state county b25003_001e b25077_001e b19013_001e {
    capture confirm variable `v'
    if _rc {
        di as error "ACS missing `v'"
        exit 459
    }
}
replace state=strtrim(state)
replace county=strtrim(county)
keep if strpos("|01|12|13|18|21|37|39|45|47|51|54|", ///
    "|"+state+"|")>0
gen str5 fips=state+county
rename b25003_001e households_2023
rename b25077_001e median_home_value_2023
rename b19013_001e median_household_income_2023
foreach v in households_2023 median_home_value_2023 median_household_income_2023 {
    destring `v', replace force
    replace `v'=. if `v'<=0
}
keep fips households_2023 median_home_value_2023 median_household_income_2023
isid fips
save "${INT}/census_acs_2023.dta", replace
