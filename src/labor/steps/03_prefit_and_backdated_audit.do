/* Diagnostics only. Revision 14 Nov 2026.
   (a) Definitive weights: IN-SAMPLE pre-event fit, Jan 2022-Aug 2024.
       The fit window now includes 2024, so these numbers describe how flat
       the pre-event gap is; they are NOT an out-of-sample validation.
   (b) Backdated weights (2022-2023 fit, same lambda, same intercept rule):
       genuinely out-of-sample Jan-Aug 2024 prediction audit and in-time
       placebo contrasts. Reported whatever they show. */
version 16.0
di as txt "03: Pre-event fit of definitive weights and backdated out-of-sample audit"
use "${MC_OUT}/data/02_paths.dta", clear
isid mdate
local tr0=ym(2022,1)
local tr1=ym(2023,12)
local va0=ym(2024,1)
local va1=ym(2024,8)
local aug=ym(2024,8)
tempname D T
postfile `D' str12 model str12 outcome double fit_rmse fit_rmse_2022_23 ///
    fit_rmse_2024 fit_max_abs_gap bd_train_rmse bd_holdout_rmse bd_rmspe_ratio ///
    bd_holdout_mean bd_aug2024_gap bd_holdout_slope bd_holdout_slope_p ///
    using "${MC_OUT}/data/03_fit_audit.dta", replace
postfile `T' str12 model str12 outcome str18 period ///
    double reference_change using "${MC_OUT}/data/03_time_placebos.dta", replace
foreach o in emp lf jobs urate {
    * (a) definitive weights, in-sample over the whole pre-event window
    quietly gen double mc_gap=mc_y_`o'-mc_syn_joint_`o'
    quietly gen double mc_sq=mc_gap^2
    quietly summarize mc_gap if inrange(mdate,`tr0',`va1'), meanonly
    assert abs(r(mean))<1e-6
    quietly summarize mc_sq if inrange(mdate,`tr0',`va1'), meanonly
    local fall=sqrt(r(mean))
    quietly summarize mc_sq if inrange(mdate,`tr0',`tr1'), meanonly
    local f23=sqrt(r(mean))
    quietly summarize mc_sq if inrange(mdate,`va0',`va1'), meanonly
    local f24=sqrt(r(mean))
    quietly gen double mc_abs=abs(mc_gap)
    quietly summarize mc_abs if inrange(mdate,`tr0',`va1'), meanonly
    local fmax=r(max)
    * (b) backdated weights: 2022-23 training, Jan-Aug 2024 never fitted
    quietly gen double mc_gbd=mc_y_`o'-mc_syn_bd_`o'
    quietly gen double mc_sqb=mc_gbd^2
    quietly summarize mc_sqb if inrange(mdate,`tr0',`tr1'), meanonly
    local tr=sqrt(r(mean))
    quietly summarize mc_sqb if inrange(mdate,`va0',`va1'), meanonly
    local va=sqrt(r(mean))
    quietly summarize mc_gbd if inrange(mdate,`va0',`va1'), meanonly
    local vm=r(mean)
    quietly summarize mc_gbd if mdate==`aug', meanonly
    local ag=r(mean)
    quietly gen int mc_t=mdate-`va0'+1
    quietly regress mc_gbd mc_t if inrange(mdate,`va0',`va1')
    local slope=_b[mc_t]
    local slopep=2*ttail(e(df_r),abs(_b[mc_t]/_se[mc_t]))
    drop mc_t
    post `D' ("joint") ("`o'") (`fall') (`f23') (`f24') (`fmax') ///
        (`tr') (`va') (`va'/`tr') (`vm') (`ag') (`slope') (`slopep')
    * In-time placebo contrasts use the BACKDATED gap (2024 untouched).
    forvalues j=1/4 {
        local ws=ym(2024,`j'+2)
        local we=`ws'+2
        local ref=ym(2024,`j')
        quietly summarize mc_gbd if inrange(mdate,`ws',`we'), meanonly
        local a=r(mean)
        quietly summarize mc_gbd if mdate==`ref', meanonly
        local b=r(mean)
        post `T' ("backdated") ("`o'") ("placebo`j'") (`a'-`b')
    }
    quietly summarize mc_gbd if inrange(mdate,ym(2024,10),ym(2024,12)), meanonly
    post `T' ("backdated") ("`o'") ("helene_octdec2024") (r(mean)-`ag')
    drop mc_gap mc_sq mc_abs mc_gbd mc_sqb
}
postclose `D'
postclose `T'
use "${MC_OUT}/data/03_fit_audit.dta", clear
assert !missing(fit_rmse,bd_train_rmse,bd_holdout_rmse,bd_aug2024_gap)
assert _N==4
isid model outcome
export delimited using "${MC_OUT}/tables/pretrend_train_holdout_ALL.csv", replace
use "${MC_OUT}/data/03_time_placebos.dta", clear
export delimited using "${MC_OUT}/tables/pre2024_time_placebos.csv", replace
use "${MC_OUT}/data/03_fit_audit.dta", clear
list outcome fit_rmse fit_rmse_2022_23 fit_rmse_2024 fit_max_abs_gap, noobs
list outcome bd_train_rmse bd_holdout_rmse bd_rmspe_ratio bd_holdout_mean ///
    bd_holdout_slope_p, noobs
quietly count if bd_rmspe_ratio>2
if r(N)>0 di as error "CAUTION: backdated (2022-23 fitted) weights have >2x out-of-sample/training RMSE for at least one outcome."
