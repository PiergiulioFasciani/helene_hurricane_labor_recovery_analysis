/* Descriptive Gaussian-kernel damage response on the joint synthetic-DiD
   weights. SAME 39 treated counties, no post-informed fitting of weights.
   A common synthetic curve shifts outcome levels, NOT HDSI dose slopes.
   Revision 14 Nov 2026 (pre-trend redesign):
   - Each county's abnormal change is its gap to the synthetic path MINUS ITS
     OWN PRE-EVENT MEAN GAP (Jan 2022-Aug 2024), not minus the single month
     of August 2024. One noisy, seasonally low month no longer sets the
     baseline of every county.
   - A pre-event PLACEBO window is added: Jan-Aug 2024 relative to each
     county's 2022-2023 mean gap. If exposure were unrelated to pre-event
     dynamics its curves are flat.
   - kernel_dose_slopes.csv reports the linear slope per SD of exposure
     (HC3) for every window, placebo included.
   Prespecified Gaussian local linear, county bootstrap and support audit. */
version 16.0
di as txt "06: Exploratory Gaussian damage-response with joint donor baseline"
local B=299
local seed=31072024
tempfile dosefile donor_paths county_abnormal all_windows all_curves

use "${MC_INPUT}", clear
keep if treated_nc==1 & mdate==ym(2024,8)
keep fips hdsi prcp_total_mm
isid fips
save `dosefile'

use "${MC_OUT}/data/02_weighted_county_month.dta", clear
* Source name map: stage 01 stores the unemployment rate in mc_raw_ur,
* not an alternate name. Keep the rate in pp and log count variables x 100.
foreach o in emp lf jobs urate unemp {
    if "`o'"=="unemp" {
        gen double ku_`o'=100*mc_raw_unemp
    }
    else if "`o'"=="urate" {
        gen double ku_`o'=mc_raw_ur
    }
    else {
        gen double ku_`o'=mc_raw_`o'
    }
}

* Strict synthetic donor mean, NEVER silently normalize missing donor weights.
preserve
    keep if mc_donor==1 & mc_w_joint>1e-10
    foreach o in emp lf jobs urate unemp {
        gen double ku_num_`o'=mc_w_joint*ku_`o' if !missing(ku_`o')
        gen double ku_den_`o'=mc_w_joint if !missing(ku_`o')
    }
    collapse (sum) ku_num_* ku_den_*, by(mdate)
    foreach o in emp lf jobs urate unemp {
        gen double ku_syn_`o'=ku_num_`o' if abs(ku_den_`o'-1)<1e-5
    }
    keep mdate ku_syn_*
    isid mdate
    save `donor_paths'
restore

keep if treated_nc==1
merge m:1 mdate using `donor_paths', assert(match) nogen
merge m:1 fips using `dosefile', assert(match) nogen
foreach o in emp lf jobs urate unemp {
    * County gap to the common synthetic path, then two county baselines:
    * full pre-event mean (main windows) and 2022-23 mean (placebo window).
    gen double ku_gap_`o'=ku_`o'-ku_syn_`o' if !missing(ku_`o',ku_syn_`o')
    bysort fips: egen double ku_pre_`o'= ///
        mean(cond(inrange(mdate,ym(2022,1),ym(2024,8)),ku_gap_`o',.))
    bysort fips: egen int ku_npre_`o'= ///
        total(inrange(mdate,ym(2022,1),ym(2024,8)) & !missing(ku_gap_`o'))
    bysort fips: egen double ku_pre23_`o'= ///
        mean(cond(inrange(mdate,ym(2022,1),ym(2023,12)),ku_gap_`o',.))
    bysort fips: egen int ku_npre23_`o'= ///
        total(inrange(mdate,ym(2022,1),ym(2023,12)) & !missing(ku_gap_`o'))
    * Baselines require the complete window; never average a partial one.
    gen double ku_abnormal_`o'=ku_gap_`o'-ku_pre_`o' ///
        if !missing(ku_gap_`o') & ku_npre_`o'==32
    gen double ku_abn_pl_`o'=ku_gap_`o'-ku_pre23_`o' ///
        if !missing(ku_gap_`o') & ku_npre23_`o'==24
}
keep fips mdate hdsi prcp_total_mm ku_abnormal_* ku_abn_pl_*
sort fips mdate
isid fips mdate
save "${MC_OUT}/data/06_treated_monthly_abnormal.dta", replace
save `county_abnormal'

local first=1
foreach window in placebo immediate short medium {
    foreach o in emp lf jobs urate unemp {
        use `county_abnormal', clear
        local start=ym(2024,10)
        local finish=ym(2024,12)
        local required=3
        local src "ku_abnormal_`o'"
        if "`window'"=="placebo" {
            local start=ym(2024,1)
            local finish=ym(2024,8)
            local required=8
            local src "ku_abn_pl_`o'"
        }
        if "`window'"=="short" {
            local start=ym(2025,1)
            local finish=ym(2025,9)
            local required=9
        }
        if "`window'"=="medium" {
            local start=ym(2025,11)
            local finish=ym(2026,9)
            local required=11
            if "`o'"=="jobs" {
                local finish=ym(2026,3)
                local required=5
            }
        }
        keep if inrange(mdate,`start',`finish')
        collapse (mean) change=`src' (count) n_months=`src' ///
            (first) hdsi prcp_total_mm, by(fips)
        gen str12 outcome="`o'"
        gen str12 window="`window'"
        gen byte expected_months=`required'
        gen byte complete=(n_months==expected_months)
        replace change=. if complete==0
        if `first' {
            save `all_windows', replace
            local first=0
        }
        else {
            append using `all_windows'
            save `all_windows', replace
        }
    }
}
use `all_windows', clear
isid fips outcome window
save "${MC_OUT}/data/06_treated_damage_windows.dta", replace
export delimited using "${MC_OUT}/tables/kernel_county_windows.csv", replace

* Linear exposure gradients for every window, the pre-event placebo included.
tempname S
postfile `S' str16 dose str12 window str12 outcome int n_counties ///
    double mean_change slope_per_sd se pvalue ///
    using "${MC_OUT}/data/06_dose_slopes.dta", replace
foreach dose in hdsi prcp_total_mm {
    foreach window in placebo immediate short medium {
        foreach o in emp lf jobs urate unemp {
            use "${MC_OUT}/data/06_treated_damage_windows.dta", clear
            keep if window=="`window'" & outcome=="`o'" & ///
                !missing(`dose',change)
            quietly count
            if r(N)<15 continue
            quietly summarize `dose'
            gen double ku_z=(`dose'-r(mean))/r(sd)
            quietly summarize change, meanonly
            local mu=r(mean)
            quietly regress change ku_z, vce(hc3)
            post `S' ("`dose'") ("`window'") ("`o'") (e(N)) (`mu') ///
                (_b[ku_z]) (_se[ku_z]) (2*ttail(e(df_r),abs(_b[ku_z]/_se[ku_z])))
        }
    }
}
postclose `S'
use "${MC_OUT}/data/06_dose_slopes.dta", clear
isid dose window outcome
export delimited using "${MC_OUT}/tables/kernel_dose_slopes.csv", replace
list dose window outcome slope_per_sd se pvalue if window=="placebo", noobs

mata:
real scalar mc_gauss_ll(real colvector x, real colvector y,
    real scalar at, real scalar h)
{
    real colvector dx, w
    real scalar s0,s1,s2,t0,t1,det
    dx = x :- at
    w = exp(-0.5 :* ((dx:/h):^2))
    s0=sum(w)
    s1=sum(w:*dx)
    s2=sum(w:*(dx:^2))
    t0=sum(w:*y)
    t1=sum(w:*dx:*y)
    det=s0*s2-s1*s1
    if (det<=1e-12 | s0<=1e-12) return(.)
    return((s2*t0-s1*t1)/det)
}
void mc_gauss_curve()
{
    real matrix d, out
    real colvector x,y,z,b,idx,sorted,w
    real scalar j,i,n,B,h,lo,hi,seed,nok
    d=st_data(.,("intensity","change"))
    x=d[,1]; y=d[,2]; n=rows(x)
    B=strtoreal(st_local("boots"))
    h=strtoreal(st_local("bandwidth"))
    lo=strtoreal(st_local("lo"))
    hi=strtoreal(st_local("hi"))
    seed=strtoreal(st_local("bootseed"))
    rseed(seed)
    z=lo :+ (hi-lo):*(0::30)/30
    out=J(31,5,.)
    for (j=1;j<=31;j++) {
        out[j,1]=z[j]
        out[j,2]=mc_gauss_ll(x,y,z[j],h)
        w=exp(-0.5:*(((x:-z[j]):/h):^2))
        out[j,5]=(sum(w)^2)/sum(w:^2)
        b=J(B,1,.)
        for (i=1;i<=B;i++) {
            idx=floor(runiform(n,1):*n):+1
            b[i]=mc_gauss_ll(x[idx],y[idx],z[j],h)
        }
        sorted=sort(select(b,b:<.),1)
        nok=rows(sorted)
        if (nok>30) {
            out[j,3]=sorted[max((1,ceil(.025*nok)))]
            out[j,4]=sorted[min((nok,ceil(.975*nok)))]
        }
    }
    st_matrix("MC_KERNEL",out)
}
end

local series=0
local first=1
foreach dose in hdsi prcp_total_mm {
    foreach window in placebo immediate short medium {
        foreach o in emp lf jobs urate unemp {
            use "${MC_OUT}/data/06_treated_damage_windows.dta", clear
            keep if window=="`window'" & outcome=="`o'" & ///
                !missing(`dose',change)
            quietly count
            local n=r(N)
            if `n'<15 {
                di as error "Kernel not estimable: `dose' / `window' / `o' N=`n'"
                continue
            }
            quietly summarize `dose', detail
            local lo=r(p10)
            local hi=r(p90)
            local bandwidth=.9*min(r(sd),(r(p75)-r(p25))/1.34)*`n'^(-.2)
            if `bandwidth'<=0 | missing(`bandwidth') {
                local bandwidth=.9*r(sd)*`n'^(-.2)
            }
            if `bandwidth'<=0 | missing(`bandwidth') | `hi'<=`lo' {
                di as error "Kernel skipped: constant dose or invalid bandwidth"
                continue
            }
            local series=`series'+1
            local boots=`B'
            local bootseed=`seed'+`series'
            rename `dose' intensity
            mata: mc_gauss_curve()
            clear
            matrix colnames MC_KERNEL = intensity fitted lower upper effective_n
            svmat double MC_KERNEL, names(col)
            gen str16 dose="`dose'"
            gen str12 window="`window'"
            gen str12 outcome="`o'"
            gen int n_counties=`n'
            gen int n_boot=`boots'
            gen double bandwidth=`bandwidth'
            gen long seed=`bootseed'
            isid intensity
            twoway (rarea lower upper intensity, color(gs13)) ///
                (line fitted intensity, lcolor(navy) lwidth(medthick)), ///
                title("Joint-weight Gaussian: `o'") ///
                subtitle("`dose' / `window'; treated counties only") ///
                xtitle("Damage intensity") ///
                ytitle("Abnormal change vs own pre-event mean gap") ///
                legend(order(2 "Local-linear" 1 "Pointwise 95% CI")) ///
                graphregion(color(white))
            graph export "${MC_OUT}/figures/kernel_`dose'_`window'_`o'.png", ///
                replace width(2000)
            graph save "${MC_OUT}/figures/kernel_`dose'_`window'_`o'.gph", replace
            if `first' {
                save `all_curves', replace
                local first=0
            }
            else {
                append using `all_curves'
                save `all_curves', replace
            }
        }
    }
}
if `first' {
    di as error "No Gaussian kernel curves: check damage and window source."
    exit 459
}
use `all_curves', clear
sort dose window outcome intensity
save "${MC_OUT}/data/06_gaussian_rebased_curves.dta", replace
export delimited using "${MC_OUT}/tables/gaussian_rebased_curves.csv", replace

di as result "06 complete: Gaussian descriptive damage intensity (pre-event-mean baseline, placebo window)."
