/* New placebo treatment groups of 39 untreated counties.
   All post-event outcomes are excluded from tuning/weight estimation.
   Same four-outcome objective and full pre-eligible donor universe as primary;
   lambda is FROZEN to the 2022-2023-selected primary joint value.
   Revision 14 Nov 2026: placebo groups are fitted exactly like the primary
   synthetic-DiD weights (Jan 2022-Aug 2024, intercept per outcome), and every
   shock is the Oct-Dec 2024 mean gap minus the PRE-EVENT MEAN gap, matching
   the primary estimand. The former holdout column is replaced by the fitted
   Jan-Aug 2024 segment error (pre2024_rmse), which is in-sample.
   Unlike randomized assignment, this is descriptive placebo benchmarking.
   A smaller group-specific eligible pool is inevitable (39 pseudo units removed).
   Report all four outcomes and fit diagnostics, not only favorable ranks. */
version 16.0
local reps=${H_JOINT_PLACEBO_REPS}
local seed=${H_JOINT_PLACEBO_SEED}
local maxit=${H_JOINT_PLACEBO_ITER}
local tol=${H_JOINT_PLACEBO_TOL}
assert `reps'>=19 & `reps'<=999
assert `maxit'>=100
assert `tol'>0 & `tol'<0.01

preserve
    use "${MC_OUT}/data/02_solver_audit.dta", clear
    quietly summarize lambda if model=="joint", meanonly
    local joint_lambda=r(mean)
    assert !missing(`joint_lambda') & `joint_lambda'>0
restore
use "${MC_OUT}/data/02_weighted_county_month.dta", clear
sort mc_uid mdate
isid mc_uid mdate
assert !missing(mc_emp,mc_lf,mc_jobs,mc_ur) if mdate<=ym(2024,12) & mdate!=ym(2024,9)
local firstmonth=ym(2022,1)
local augpos=ym(2024,8)-`firstmonth'+1
local postfirst=ym(2024,10)-`firstmonth'+1
local postlast=ym(2024,12)-`firstmonth'+1
foreach fn in hj_project_simplex hj_fista hj_run {
    capture mata: mata drop `fn'()
}
mata:
real colvector hj_project_simplex(real colvector v)
{
    real colvector u
    real scalar i, acc, rho, theta
    u = -sort(-v,1)
    acc = 0
    rho = 0
    theta = 0
    for (i=1; i<=rows(u); i++) {
        acc = acc+u[i]
        if (u[i]>(acc-1)/i) {
            rho = i
            theta = (acc-1)/i
        }
    }
    if (rho==0) {
        errprintf("Simplex projection failed.\n")
        exit(459)
    }
    return((v:-theta):*(v:>theta))
}

// Accelerated projected gradient on 0.5||A*w-y||^2 + lambda/2||w-u||^2
// subject to w>=0, sum(w)=1. Last 3 output entries: iters, dual gap,
// convergence flag. The gap max_g - min_g is a simplex optimality bound.
real colvector hj_fista(real matrix A, real colvector y,
    real scalar lambda, real scalar maxit, real scalar tol)
{
    real colvector c, v, w, z, wnew, g
    real scalar p, i, it, L, t, tnew, gap, status
    p = cols(A)
    c = quadcross(A,y) :+ lambda/p
    v = J(p,1,1/sqrt(p))
    for (i=1; i<=60; i++) {
        v = quadcross(A,A*v)+lambda*v
        v = v/sqrt(quadcross(v,v))
    }
    L = 1.05*quadcross(v,quadcross(A,A*v)+lambda*v)
    if (L<=1e-12) L=1
    w = J(p,1,1/p)
    z = w
    t = 1
    gap = .
    status = 1
    for (it=1; it<=maxit; it++) {
        g = quadcross(A,A*z)+lambda*z-c
        wnew = hj_project_simplex(z-g/L)
        tnew = (1+sqrt(1+4*t^2))/2
        z = wnew + ((t-1)/tnew)*(wnew-w)
        w = wnew
        t = tnew
        if (mod(it,50)==0 | it==maxit) {
            g = quadcross(A,A*w)+lambda*w-c
            gap = max((0,quadcross(g,w)-min(g)))
            if (gap<=tol) {
                status = 0
                break
            }
        }
    }
    if (it>maxit) it=maxit
    w = w:*(w:>0)
    w = w/sum(w)
    g = quadcross(A,A*w)+lambda*w-c
    gap = max((0,quadcross(g,w)-min(g)))
    status=1
    if (gap<=tol) status=0
    return(w \ it \ gap \ status)
}

void hj_run()
{
    real matrix dat, X, draws, groups, A, rr
    real colvector units, step, donor, trained, y, wobs, idx, dst, grp
    real colvector colsleft, weights, z, fit, diff, scales, index, gout
    real colvector pout, scaled, pre, post, aug, mask, ridx, obsdraw, observed_fit, holdout_fit
    real scalar n, nc, u, o, k, j, r, B, itmax, tol, lam, no, row, valfit
    real scalar ntreated, seed, eff, before, startidx, nbad

    dat=st_data(.,tokens("mc_emp mc_lf mc_jobs mc_ur"))
    units=st_data(.,"mc_uid")
    step=st_data(.,"mc_step")
    donor=st_data(.,"mc_donor")
    trained=st_data(.,"mc_w_joint")
    nc=rows(dat)/57
    if (nc!=floor(nc)) exit(459)
    n=sum((step:==32):&(donor:==1))
    ntreated=sum((step:==32):&(donor:==0))
    if (ntreated!=39 | n<500) exit(459)
    X=J(228,n,.)
    y=J(228,1,0)
    wobs=J(n,1,0)
    j=0
    for(u=1; u<=nc; u++) {
        ridx=(u-1)*57 :+ (1::57)
        if (sum(step[ridx]:==(1::57))!=57) exit(459)
        if (donor[ridx[1]]==1) {
            j++
            wobs[j]=trained[ridx[1]]
        }
        for (o=1; o<=4; o++) {
            idx=(o-1)*57 :+ (1::57)
            if (donor[ridx[1]]==1) X[idx,j]=dat[ridx,o]
            else y[idx]=y[idx]+dat[ridx,o]
        }
    }
    if(j!=n | abs(sum(wobs)-1)>0.0001) exit(459)
    y=y/39
    B=strtoreal(st_local("reps"))
    seed=strtoreal(st_local("seed"))
    itmax=strtoreal(st_local("maxit"))
    tol=strtoreal(st_local("tol"))
    lam=strtoreal(st_local("joint_lambda"))
    rseed(seed)
    draws=J(B,11,.)
    groups=J(B*39,2,.)
    rr=J(4,3,.)
    /* Observed changes from frozen common-weight SCM. */
    for (o=1; o<=4; o++) {
        idx=(o-1)*57 :+ (1::57)
        diff=y[idx]-X[idx,.]*wobs
        diff=diff:-mean(diff[1::32])
        rr[o,1]=mean(diff[strtoreal(st_local("postfirst"))::strtoreal(st_local("postlast"))])
        rr[o,2]=sqrt(mean(diff[1::32]:^2))
        rr[o,3]=sqrt(mean(diff[25::32]:^2))
    }
    st_matrix("HJ_OBS",rr)
    /* Real treated pre-fit on the SAME four-outcome, 2022-scaled loss
       as the placebo fit screen (32 pre-event months, intercept removed).
       Scales are learned from treated 2022. */
    observed_fit=J(128,1,.)
    for(o=1; o<=4; o++) {
        idx=(o-1)*57 :+ (1::12)
        scaled=y[idx]
        if(o==4) before=.15
        else before=.25
        before=max((before,sqrt(mean((scaled:-mean(scaled)):^2))))
        idx=(o-1)*57 :+ (1::32)
        diff=y[idx]-X[idx,.]*wobs
        observed_fit[((o-1)*32+1)::(o*32)]=(diff:-mean(diff))/(before*sqrt(128))
    }
    st_numscalar("HJ_OBS_FIT",sqrt(mean(observed_fit:^2)))
    for(r=1; r<=B; r++) {
        index=order(runiform(n,1),1)
        grp=index[1::39]
        groups[((r-1)*39+1)::(r*39),.] = (J(39,1,r),grp)
        mask=J(n,1,1)
        mask[grp]=J(39,1,0)
        colsleft=selectindex(mask)
        pout=rowsum(X[,grp]):/39
        scales=J(4,1,.)
        A=J(128,rows(colsleft),.)
        z=J(128,1,.)
        for(o=1; o<=4; o++) {
            idx=(o-1)*57 :+ (1::12)
            scaled=pout[idx]
            if(o==4) before=.15
            else before=.25
            scales[o]=max((before,sqrt(mean((scaled:-mean(scaled)):^2))))
            idx=(o-1)*57 :+ (1::32)
            startidx=(o-1)*32
            A[(startidx+1)::(startidx+32),.]=(X[idx,colsleft]:-mean(X[idx,colsleft]))/(scales[o]*sqrt(128))
            z[(startidx+1)::(startidx+32)]=(pout[idx]:-mean(pout[idx]))/(scales[o]*sqrt(128))
        }
        fit=hj_fista(A,z,lam,itmax,tol)
        weights=fit[1::rows(colsleft)]
        diff=pout-X[,colsleft]*weights
        gout=J(4,1,.)
        holdout_fit=J(32,1,.)
        for(o=1; o<=4; o++) {
            idx=(o-1)*57 :+ (1::57)
            scaled=diff[idx]
            scaled=scaled:-mean(scaled[1::32])
            holdout_fit[((o-1)*8+1)::(o*8)]=scaled[25::32]/(scales[o]*sqrt(128))
            gout[o]=mean(scaled[strtoreal(st_local("postfirst"))::strtoreal(st_local("postlast"))])
        }
        valfit=sqrt(mean(holdout_fit:^2))
        draws[r,.]=(r,gout',sqrt(mean((A*weights-z):^2)),valfit,fit[rows(fit)-2],fit[rows(fit)-1],fit[rows(fit)],sum(weights:>1e-8))
    }
    st_matrix("HJ_DRAWS",draws)
    st_matrix("HJ_GROUPS",groups)
}
end

mata: hj_run()
matrix colnames HJ_DRAWS = draw shock_emp shock_lf shock_jobs shock_urate ///
    prefit_rmse pre2024_rmse iterations dual_gap status n_active
matrix colnames HJ_OBS = observed_shock pre_rmse pre2024_rmse
matrix rownames HJ_OBS = emp lf jobs urate
matrix colnames HJ_GROUPS = draw donor_rank

preserve
    clear
    svmat double HJ_DRAWS, names(col)
    assert _N==`reps'
    gen byte converged=(status==0 & dual_gap<=`tol')
    gen double observed_joint_prefit=scalar(HJ_OBS_FIT)
    gen double lambda_locked=`joint_lambda'
    gen long seed=`seed'
    gen byte pseudo_group_size=39
    gen double kkt_tolerance=`tol'
    gen int max_iterations=`maxit'
    save "${MC_OUT}/placebos/joint_placebo_draws.dta", replace
    export delimited using "${MC_OUT}/placebos/joint_placebo_draws.csv", replace
    quietly count if converged==0
    local nbad=r(N)
    quietly summarize prefit_rmse, meanonly
    local obsjoint=scalar(HJ_OBS_FIT)
restore

/* Empirical ranks by outcome and fixed fit screens; NO selection by significance. */
tempname S
postfile `S' str12 outcome str12 screen double observed_shock ///
    int n_placebos n_extreme double left_rank right_rank ///
    using "${MC_OUT}/placebos/joint_placebo_summary.dta", replace
preserve
    clear
    svmat double HJ_DRAWS, names(col)
    forvalues o=1/4 {
        local yn="emp"
        if `o'==2 local yn="lf"
        if `o'==3 local yn="jobs"
        if `o'==4 local yn="urate"
        local obs=el(HJ_OBS,`o',1)
        local base=el(HJ_OBS,`o',2)
        foreach bound in 0 1 1.5 2 {
            local condition "!missing(prefit_rmse)"
            local screen "all"
            /* Fit screens use one shared four-outcome scale; require only
               overall convergence for publication-ready ranks. */
            if `bound'>0 {
                local condition "prefit_rmse <= (`bound'*`obsjoint')"
                local screen "`bound'x_treated"
            }
            quietly count if `condition' & status==0
            local n=r(N)
            quietly count if `condition' & status==0 & shock_`yn'<=`obs'
            local k=r(N)
            local left=.
            if `n'>0 local left=(1+`k')/(1+`n')
            quietly count if `condition' & status==0 & abs(shock_`yn')>=abs(`obs')
            local a=r(N)
            local right=.
            if `n'>0 local right=(1+`a')/(1+`n')
            post `S' ("`yn'") ("`screen'") (`obs') (`n') (`k') (`left') (`right')
        }
    }
restore
postclose `S'
use "${MC_OUT}/placebos/joint_placebo_summary.dta", clear
isid outcome screen
export delimited using "${MC_OUT}/placebos/joint_placebo_summary.csv", replace
list if screen=="all", noobs

/* Map sampled donor ranks to original county FIPS for full traceability. */
tempfile donorid
use "${MC_OUT}/data/02_weighted_county_month.dta", clear
keep if mc_donor==1 & mc_step==32
sort fips
gen double donor_rank=_n
keep donor_rank fips
isid donor_rank
save `donorid'
clear
svmat double HJ_GROUPS, names(col)
* Unsampled eligible donors may exist only in the lookup; sampled ranks must match.
merge m:1 donor_rank using `donorid', assert(match using) keep(match) nogen
assert _N==`reps'*39
isid draw donor_rank
export delimited using "${MC_OUT}/placebos/joint_placebo_counties.csv", replace

if `nbad'>0 {
    di as error "WARNING: `nbad' of `reps' joint placebo fits did not meet `tol' tolerance."
    if ${H_JOINT_PLACEBO_STRICT}==1 exit 459
}
di as result "Joint placebo stage complete. This is benchmarking, NOT random-assignment inference."
