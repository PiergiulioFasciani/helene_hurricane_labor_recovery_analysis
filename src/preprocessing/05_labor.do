/* National county-month LAUS + QCEW, reconstructed from source CSVs.
   Default LAUS source: documented BLS-data mirror used in our first working
   automated pipeline. Optional original BLS time-series files.
   Neither path accepts a prebuilt .dta, nor imputes October 2025 LAUS. */
version 16.0
di as txt "05: National labor-market panel"
* Strict source-record build. No prebuilt national panel is accepted.
local labor_source "${LABOR_SOURCE}"
if "`labor_source'"!="rebuild" {
    di as error "This repository requires LABOR_SOURCE=rebuild from raw BLS data."
    exit 459
}
if "`labor_source'"=="rebuild" {
    * BLS often returns HTTP 403 for automated historical LAUS flat files.
    * Manually place official files under data/raw/LAUS/ when possible.
    capture mkdir "${RAW}/LAUS"
    capture mkdir "${RAW}/QCEW"
    local official 1
    foreach file in la.area la.data.0.CurrentU20-24 la.data.0.CurrentU25-29 {
        capture confirm file "${RAW}/LAUS/`file'"
        if _rc local official 0
    }
    * Source mode is explicit and preserved in the source manifest.
    * "mirror" reproduces the previously successful automated LAUS route.
    * "official" requires official BLS flat files (automatic download attempted).
    local mode "${LAUS_SOURCE_MODE}"
    if !inlist("`mode'", "mirror", "official") {
        di as error "LAUS_SOURCE_MODE must be mirror or official."
        exit 198
    }
    global LAUS_SOURCE_DESC ""
    global LAUS_SOURCE_FILE ""
    global LAUS_SOURCE_URL ""
    if "`mode'"=="official" {
        foreach file in la.area la.data.0.CurrentU20-24 la.data.0.CurrentU25-29 {
            capture confirm file "${RAW}/LAUS/`file'"
            if _rc {
                local official_url "https://download.bls.gov/pub/time.series/la/`file'"
                di as txt "Downloading official LAUS `file' ..."
                capture noisily copy "`official_url'" "${RAW}/LAUS/`file'", replace
                if _rc {
                    capture noisily shell /usr/bin/curl -fL --retry 2 ///
                        --max-time 300 -o "${RAW}/LAUS/`file'" "`official_url'"
                }
            }
            capture confirm file "${RAW}/LAUS/`file'"
            if _rc {
                di as error "Official BLS LAUS `file' unavailable (possibly HTTP 403)."
                di as error "Use LAUS_SOURCE_MODE=mirror to replicate the original automated pipeline."
                exit 601
            }
        }
        local official 1
        global LAUS_SOURCE_DESC "Official BLS LAUS flat files"
        global LAUS_SOURCE_FILE "${RAW}/LAUS/la.data.0.CurrentU20-24; ${RAW}/LAUS/la.data.0.CurrentU25-29"
        global LAUS_SOURCE_URL "https://download.bls.gov/pub/time.series/la/"
    }
    else {
        local official 0
        local mirror "${RAW}/LAUS/LAUS_mirror_all_counties.csv"
        local mirror_url "https://api.economicdata.dev/v1/laus/download/laus.csv"
        di as result "LAUS source: BLS data redistributed by economicdata.dev (third party)."
        capture confirm file "`mirror'"
        if _rc {
            local part "`mirror'.part"
            capture erase "`part'"
            di as txt "Downloading mirror LAUS raw CSV into repository data/raw/LAUS/ ..."
            capture noisily copy "`mirror_url'" "`part'", replace
            if _rc {
                capture erase "`part'"
                capture noisily shell /usr/bin/curl -fL --retry 3 ///
                    --connect-timeout 30 --max-time 900 ///
                    -o "`part'" "`mirror_url'"
            }
            capture confirm file "`part'"
            if _rc {
                di as error "LAUS mirror download failed: `mirror_url'"
                exit 601
            }
            tempname mirrorhead
            file open `mirrorhead' using "`part'", read text
            file read `mirrorhead' header
            file close `mirrorhead'
            if strpos(lower(`"`header'"'),"series")==0 | ///
                strpos(lower(`"`header'"'),"location")==0 | ///
                strpos(lower(`"`header'"'),"labor_force")==0 {
                di as error "Downloaded LAUS response is not the expected raw CSV."
                capture erase "`part'"
                exit 459
            }
            copy "`part'" "`mirror'", replace
            capture erase "`part'"
        }
        * Validate an existing file too; an HTTP error page is not a dataset.
        tempname mirrorhead2
        file open `mirrorhead2' using "`mirror'", read text
        file read `mirrorhead2' header
        file close `mirrorhead2'
        if strpos(lower(`"`header'"'),"series")==0 | ///
            strpos(lower(`"`header'"'),"location")==0 | ///
            strpos(lower(`"`header'"'),"labor_force")==0 {
            di as error "Cached LAUS source is invalid. Check `mirror'."
            exit 459
        }
        * Metadata snapshot is informative, but not mandatory for data import.
        local metadata "${RAW}/LAUS/LAUS_mirror_metadata.json"
        capture confirm file "`metadata'"
        if _rc {
            capture noisily shell /usr/bin/curl -fL --retry 1 ///
                --connect-timeout 20 --max-time 60 ///
                -o "`metadata'" "https://api.economicdata.dev/v1/laus"
        }
        global LAUS_SOURCE_DESC "Third-party economicdata.dev LAUS bulk CSV; upstream BLS LAUS"
        global LAUS_SOURCE_FILE "`mirror'"
        global LAUS_SOURCE_URL "`mirror_url'"
    }
    tempfile laus_long laus_data
    if `official' {
        local first 1
        foreach file in la.data.0.CurrentU20-24 la.data.0.CurrentU25-29 {
            import delimited using "${RAW}/LAUS/`file'", clear ///
                varnames(1) delimiters(tab) stringcols(1 3) bindquote(nobind)
            foreach req in series_id year period value {
                capture confirm variable `req'
                if _rc exit 459
            }
            replace series_id=strtrim(series_id)
            keep if substr(series_id,1,5)=="LAUCN"
            gen str5 fips=substr(series_id,6,5)
            gen str2 measure=substr(series_id,strlen(series_id)-1,2)
            keep if inlist(measure,"03","04","05","06")
            capture confirm numeric variable year
            if _rc destring year, replace
            gen byte month=real(substr(strtrim(period),2,2))
            keep if inrange(year,2022,2026) & inrange(month,1,12)
            keep if ym(year,month)<=ym(2026,9)
            capture confirm numeric variable value
            if _rc gen double laus_value=real(subinstr(strtrim(value),",","",.))
            else gen double laus_value=value
            keep fips year month measure laus_value
            if `first' {
                save `laus_long', replace
                local first 0
            }
            else {
                append using `laus_long'
                save `laus_long', replace
            }
        }
        use `laus_long', clear
        isid fips year month measure
        reshape wide laus_value, i(fips year month) j(measure) string
        rename laus_value03 laus_unemployment_rate
        rename laus_value04 laus_unemployed
        rename laus_value05 laus_employed
        rename laus_value06 laus_labor_force
        gen int mdate=ym(year,month)
        keep fips mdate laus_*
        gen str80 county_name=""
    }
    else {
        * This is the original automatic path used by the working pipeline.
        import delimited using "`mirror'", clear varnames(1) bindquote(strict)
        foreach req in series location date labor_force employed unemployed unemployment_rate {
            capture confirm variable `req'
            if _rc {
                di as error "Third-party LAUS raw CSV lacks `req': `mirror'"
                exit 459
            }
        }
        keep if strtrim(series)=="nsa" & ///
            substr(strtrim(location),1,10)=="us-county-"
        gen str5 fips=substr(strtrim(location),11,5)
        gen int year=real(substr(date,1,4))
        gen byte month=real(substr(date,6,2))
        keep if strlen(fips)==5 & !missing(real(fips)) & ///
            inrange(year,2022,2026) & inrange(month,1,12) & ///
            ym(year,month)<=ym(2026,9)
        foreach v in labor_force employed unemployed unemployment_rate {
            capture confirm numeric variable `v'
            if _rc destring `v', replace ignore(",") force
        }
        rename labor_force laus_labor_force
        rename employed laus_employed
        rename unemployed laus_unemployed
        rename unemployment_rate laus_unemployment_rate
        gen int mdate=ym(year,month)
        gen str80 county_name=""
        keep fips mdate county_name laus_*
    }
    isid fips mdate
    save `laus_data', replace

    * QCEW: industry 10, own_code=0, agg level 70, three monthly job counts.
    tempfile all_qcew
    local firstq 1
    forvalues yy=2022/2026 {
        local maxq 4
        if `yy'==2026 local maxq ${QCEW_2026_MAXQ}
        forvalues qq=1/`maxq' {
            local file "${RAW}/QCEW/QCEW_`yy'_q`qq'_industry10.csv"
            capture confirm file "`file'"
            if _rc {
                local url "https://data.bls.gov/cew/data/api/`yy'/`qq'/industry/10.csv"
                capture noisily copy "`url'" "`file'", replace
                if _rc shell /usr/bin/curl -fL --retry 3 -o "`file'" "`url'"
            }
            capture confirm file "`file'"
            if _rc {
                di as error "QCEW `yy'Q`qq' download missing; manually obtain `file'."
                exit 601
            }
            import delimited using "`file'", clear varnames(1) ///
                stringcols(_all) bindquote(strict)
            foreach req in area_fips own_code industry_code agglvl_code ///
                month1_emplvl month2_emplvl month3_emplvl {
                capture confirm variable `req'
                if _rc {
                    di as error "Invalid QCEW `yy'Q`qq': `req' is missing"
                    exit 459
                }
            }
            keep if strtrim(own_code)=="0" & strtrim(industry_code)=="10" ///
                & strtrim(agglvl_code)=="70"
            gen str5 fips=strtrim(area_fips)
            keep if strlen(fips)==5 & !missing(real(fips))
            rename month1_emplvl qcew_jobs1
            rename month2_emplvl qcew_jobs2
            rename month3_emplvl qcew_jobs3
            foreach v in qcew_jobs1 qcew_jobs2 qcew_jobs3 {
                destring `v', replace ignore(",") force
            }
            * IMPORTANT: whitelist just the four canonical QCEW variables.
            * This removes original BLS 'year', 'quarter', and 35+ other fields,
            * making the next generate statements safe on every quarter.
            keep fips qcew_jobs1 qcew_jobs2 qcew_jobs3
            gen int year=`yy'
            gen byte quarter=`qq'
            isid fips year quarter
            di as txt "QCEW period identifiers rebuilt for `yy'Q`qq'"
            if `firstq' {
                save `all_qcew', replace
                local firstq 0
            }
            else {
                append using `all_qcew'
                save `all_qcew', replace
            }
        }
    }
    use `all_qcew', clear
    isid fips year quarter
    reshape long qcew_jobs, i(fips year quarter) j(month_in_quarter)
    gen int mdate=ym(year,3*(quarter-1)+month_in_quarter)
    keep fips mdate qcew_jobs
    isid fips mdate
    tempfile qcew_monthly
    save `qcew_monthly', replace
    use `laus_data', clear
    merge 1:1 fips mdate using `qcew_monthly', keep(master match) nogen
}

* Enforce complete 57-month calendar by county, even when LAUS is missing.
keep if inrange(mdate,ym(2022,1),ym(2026,9))
isid fips mdate
tempfile labor_obs skeleton
save `labor_obs', replace
keep fips
duplicates drop fips, force
expand 57
bysort fips: gen int mdate=ym(2022,1)+_n-1
format mdate %tm
keep fips mdate
isid fips mdate
save `skeleton', replace
merge 1:1 fips mdate using `labor_obs', keep(master match) nogen
foreach v in laus_employed laus_unemployed laus_labor_force laus_unemployment_rate qcew_jobs {
    capture confirm variable `v'
    if _rc {
        di as error "Missing labor variable `v'"
        exit 111
    }
}
* Oct 2025 LAUS: the CPS was not collected during the 2025 shutdown.
* Some cached/mirror sources contain non-missing Oct values. Preserve these
* source-supplied values separately for audit; NEVER treat them as official
* observed LAUS in the main analysis or substitute zeros/interpolation.
* QCEW jobs are independent workplace measures and are NOT touched here.
gen byte laus_oct2025_unavailable=(mdate==ym(2025,10))
gen byte laus_oct2025_source_flag=0
foreach v in laus_employed laus_unemployed laus_labor_force laus_unemployment_rate {
    replace laus_oct2025_source_flag=1 if ///
        laus_oct2025_unavailable & !missing(`v')
}
quietly count if laus_oct2025_source_flag==1
local n_oct2025_source=r(N)
di as txt "October 2025 LAUS rows containing source-supplied values: `n_oct2025_source'"
if `n_oct2025_source'>0 {
    preserve
        keep if laus_oct2025_source_flag==1
        keep fips mdate laus_employed laus_unemployed ///
            laus_labor_force laus_unemployment_rate qcew_jobs
        sort fips mdate
        export delimited using ///
            "${OUT}/audits/LAUS_Oct2025_source_values_excluded.csv", replace
    restore
    di as txt "Source values archived: outputs/audits/LAUS_Oct2025_source_values_excluded.csv"
}
foreach v in laus_employed laus_unemployed laus_labor_force laus_unemployment_rate {
    quietly replace `v'=. if laus_oct2025_unavailable
}
assert missing(laus_employed) & missing(laus_unemployed) & ///
    missing(laus_labor_force) & missing(laus_unemployment_rate) ///
    if laus_oct2025_unavailable
label var laus_oct2025_source_flag ///
    "Input had Oct 2025 LAUS values; archived and excluded from main panel"
gen byte laus_available=!missing(laus_employed,laus_unemployed,laus_labor_force)
gen byte qcew_available=!missing(qcew_jobs)
* Recreate period variables even if a cached source already contained them.
capture drop year month
gen int year=year(dofm(mdate))
gen byte month=month(dofm(mdate))
gen double ln_laus_employed=ln(laus_employed) if laus_employed>0
gen double ln_qcew_jobs=ln(qcew_jobs) if qcew_jobs>0
egen long county_id=group(fips)
format mdate %tm
isid fips mdate
xtset county_id mdate
assert mdate>=ym(2022,1) & mdate<=ym(2026,9)
* Catch changes in geography rather than altering the donor universe silently.
assert _N==57*3221
save "${INT}/labor_national_2022_2026m09.dta", replace
preserve
    collapse (count) n_laus=laus_employed n_qcew=qcew_jobs, by(mdate)
    format mdate %tm
    export delimited using "${OUT}/audits/labor_monthly_coverage.csv", replace
restore
