/* FEMA verified property loss county totals (disasters 4827-4832).
   CAUTION: assistance applicants/verified losses, not total storm damage. */
version 16.0
di as txt "03: FEMA verified losses"
* The FEMA API data (or local raw official API CSV snapshots) are mandatory.
* No preaggregated FEMA county data and no previous severity master are accepted.
local source "${FEMA_SOURCE}"
if "`source'"!="api" {
    di as error "FEMA_SOURCE=api is required for source-first reproducibility."
    exit 459
}
if "`source'"=="api" {
    tempfile fema_rows
    local first=1
    foreach num in 4827 4828 4829 4830 4831 4832 {
        local cache "${RAW}/FEMA_`num'_verified.csv"
        capture confirm file "`cache'"
        if _rc {
            local url "https://www.fema.gov/api/open/v2/IndividualsAndHouseholdsProgramValidRegistrations?%24filter=disasterNumber%20eq%20`num'&%24select=disasterNumber,fips,rpfvl,ppfvl&%24format=csv&%24metadata=off&%24allrecords=true"
            capture noisily copy "`url'" "`cache'", replace
            if _rc {
                di as error "FEMA query for disaster `num' failed. Place the official raw FEMA v2 CSV for that disaster under data/raw/ and rerun."
                exit 603
            }
        }
        import delimited using "`cache'", clear varnames(1) stringcols(2)
        foreach req in disasternumber fips rpfvl ppfvl {
            capture confirm variable `req'
            if _rc {
                di as error "FEMA CSV for disaster `num' missing `req'"
                exit 459
            }
        }
        capture confirm numeric variable disasternumber
        if _rc destring disasternumber, replace
        keep if disasternumber==`num'
        keep fips rpfvl ppfvl
        if `first' {
            save `fema_rows', replace
            local first=0
        }
        else {
            append using `fema_rows'
            save `fema_rows', replace
        }
    }
    use `fema_rows', clear
    capture confirm numeric variable fips
    if !_rc tostring fips, replace format(%05.0f)
    replace fips=strtrim(fips)
    foreach v in rpfvl ppfvl {
        capture confirm numeric variable `v'
        if _rc destring `v', replace ignore(",") force
    }
    gen byte rp_obs=!missing(rpfvl)
    gen byte pp_obs=!missing(ppfvl)
    gen double property_damage_usd=rpfvl
    gen double personal_damage_usd=ppfvl
    keep if strlen(fips)==5 & !missing(real(fips))
    collapse (sum) property_damage_usd personal_damage_usd ///
        (sum) rp_obs pp_obs, by(fips)
    replace property_damage_usd=. if rp_obs==0
    replace personal_damage_usd=. if pp_obs==0
    drop rp_obs pp_obs
}
else {
    di as error "Unknown FEMA_SOURCE `source'"
    exit 198
}
foreach v in fips property_damage_usd personal_damage_usd {
    capture confirm variable `v'
    if _rc {
        di as error "FEMA source missing `v'"
        exit 111
    }
}
capture confirm numeric variable fips
if !_rc tostring fips, replace format(%05.0f)
replace fips=strtrim(fips)
foreach v in property_damage_usd personal_damage_usd {
    capture confirm numeric variable `v'
    if _rc destring `v', replace ignore(",") force
}
keep fips property_damage_usd personal_damage_usd
keep if strlen(fips)==5 & !missing(real(fips))
isid fips
* The original pipeline treats both loss fields as required for complete
* dollar-denominated sums. Missing is NEVER turned into zero automatically.
gen double total_damage_usd=property_damage_usd+personal_damage_usd
save "${INT}/fema_helene_county.dta", replace
