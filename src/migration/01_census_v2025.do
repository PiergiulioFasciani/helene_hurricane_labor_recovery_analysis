/* Official Census Vintage 2025 county population + net domestic migration.
   County/year panel is intentionally separate from monthly LAUS/QCEW.
   A missing source is flagged and never silently replaced with made-up data. */
version 16.0
local source "${H_CENSUS_FILE}"
if "`source'"=="" local source "${H_ROOT}/data/raw/census/co-est2025-alldata.csv"
capture confirm file "`source'"
if _rc {
    local url "https://www2.census.gov/programs-surveys/popest/datasets/2020-2025/counties/totals/co-est2025-alldata.csv"
    di as text "Downloading Census CO-EST2025-ALLDATA (July 2025 release vintage)."
    capture noisily copy "`url'" "${H_ROOT}/data/raw/census/co-est2025-alldata.csv", replace
    capture confirm file "${H_ROOT}/data/raw/census/co-est2025-alldata.csv"
    if _rc {
        capture noisily shell /usr/bin/curl -fL --retry 2 --connect-timeout 8 --max-time 50 -o "${H_ROOT}/data/raw/census/co-est2025-alldata.csv" "`url'"
    }
    local source "${H_ROOT}/data/raw/census/co-est2025-alldata.csv"
}
capture confirm file "`source'"
if _rc {
    tempname h
    file open `h' using "${H_ROOT}/outputs/audit/CENSUS_STATUS.txt", write text replace
    file write `h' "UNAVAILABLE: Official Census CO-EST2025-ALLDATA was not accessible. Annual migration NOT estimated." _n
    file write `h' "Manual source: https://www.census.gov/data/datasets/time-series/demo/popest/2020s-counties-total.html" _n
    file close `h'
    di as error "Census annual source missing; no migration data can be computed."
    exit 601
}
import delimited using "`source'", clear varnames(1) stringcols(_all) bindquote(strict)
foreach c in sumlev state county popestimate2024 popestimate2025 domesticmig2025 {
    capture confirm variable `c'
    if _rc {
        di as error "Census file does not contain required header `c'. STOP - cannot silently invent migration."
        exit 111
    }
}
keep if real(sumlev)==50
assert !missing(real(state),real(county))
gen str5 fips=string(real(state),"%02.0f")+string(real(county),"%03.0f")
isid fips
forvalues y=2020/2025 {
    capture confirm variable popestimate`y'
    if _rc exit 111
    gen double pop`y'=real(subinstr(popestimate`y',",","",.))
    gen double dom`y'=.
    gen double intl`y'=.
    gen double net`y'=.
    if `y'>=2021 {
        foreach c in domesticmig internationalmig netmig {
            capture confirm variable `c'`y'
            if _rc {
                di as error "Missing required migration component `c'`y'"
                exit 111
            }
        }
        replace dom`y'=real(subinstr(domesticmig`y',",","",.))
        replace intl`y'=real(subinstr(internationalmig`y',",","",.))
        replace net`y'=real(subinstr(netmig`y',",","",.))
    }
}
keep fips pop20* dom20* intl20* net20*
reshape long pop dom intl net, i(fips) j(year)
rename pop population_july1
rename dom domestic_migration_net
rename intl international_mig_net
rename net total_migration_net
bysort fips (year): gen double population_change=population_july1-population_july1[_n-1] if _n>1
bysort fips (year): gen double pop_growth_pct=100*population_change/population_july1[_n-1] if _n>1 & population_july1[_n-1]>0
bysort fips (year): gen double netdom_per1000=1000*domestic_migration_net/population_july1[_n-1] if _n>1 & population_july1[_n-1]>0
label var domestic_migration_net "Census NET domestic migration, July prior year to July current year"
label var population_july1 "Census PEP annual July 1 resident population"
format netdom_per1000 %9.2f
isid fips year
save "${H_ROOT}/outputs/migration/census_county_year.dta", replace
export delimited using "${H_ROOT}/outputs/migration/census_county_year.csv", replace

* Assign the SAME frozen counties/controls using only pre-treatment roster.
tempfile assign
preserve
    use "${H_ROOT}/data/input/Helene_NC_ES_preinference_panel.dta", clear
    keep if mdate==ym(2024,8)
    keep fips treated_nc es_donor_pool
    isid fips
    save `assign'
restore
merge m:1 fips using `assign', keep(match) nogen
preserve
    keep if treated_nc==1 | es_donor_pool==1
    gen byte grp=treated_nc==1
    collapse (count) counties=population_july1 ///
        (mean) mean_pop_growth_pct=pop_growth_pct ///
          mean_domestic_net=domestic_migration_net ///
          mean_netdom_per1000=netdom_per1000, by(grp year)
    label var grp "1 treated NC; 0 geographically eligible control county"
    export delimited using "${H_ROOT}/outputs/migration/census_group_comparison.csv", replace
restore
/* Annual July-to-July migration / population exploratory event contrasts.
   2024 is the last UNEXPOSED July year, 2025 is partly post-Helene.
   These are county fixed-effect contrasts, not evidence of gross exits. */
preserve
    keep if (treated_nc==1 | es_donor_pool==1) & inrange(year,2021,2025)
    egen long demo_id=group(fips)
    foreach yr in 2021 2022 2023 2025 {
        gen byte d`yr'=treated_nc==1 & year==`yr'
    }
    tempname R
    postfile `R' str20 outcome double b2025 se2025 lo95 hi95 p2025 ///
        pretest_p long N_counties N_obs using ///
        "${H_ROOT}/outputs/migration/census_annual_effects.dta", replace
    foreach y in netdom_per1000 pop_growth_pct {
        quietly count if !missing(`y')
        if r(N)>100 {
            capture quietly areg `y' d2021 d2022 d2023 d2025 i.year, ///
                absorb(demo_id) vce(cluster demo_id)
            if !_rc {
                local b=_b[d2025]
                local se=_se[d2025]
                local df=e(df_r)
                local critical=invttail(`df',.025)
                local lo=`b'-`critical'*`se'
                local hi=`b'+`critical'*`se'
                local p=2*ttail(`df',abs(`b'/`se'))
                local nobs=e(N)
                local ncounty=e(N_clust)
                local pp=.
                capture quietly test d2021 d2022 d2023
                if !_rc local pp=r(p)
                post `R' ("`y'") (`b') (`se') (`lo') (`hi') (`p') ///
                    (`pp') (`ncounty') (`nobs')
            }
        }
    }
    postclose `R'
    use "${H_ROOT}/outputs/migration/census_annual_effects.dta", clear
    export delimited using "${H_ROOT}/outputs/migration/census_annual_effects.csv", replace
restore
quietly count if treated_nc==1 & year==2025
local n=r(N)
tempname h
file open `h' using "${H_ROOT}/outputs/audit/CENSUS_STATUS.txt", write text replace
file write `h' "AVAILABLE: Official Census Vintage 2025 county annual population/migration." _n
file write `h' "2025 treated counties with matched FIPS: `n' of 39." _n
file write `h' "Component DOMESTICMIG2025 refers to July 2024-July 2025, including part of Helene aftermath." _n
file write `h' "Net domestic migration cannot distinguish more departures from fewer arrivals." _n
file close `h'
di as result "Census county-year migration output ready; `n' of 39 treated counties covered in 2025."
