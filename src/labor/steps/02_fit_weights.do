/* DEFINITIVE JOINT-ONLY estimator, synthetic-DiD unit weights.
   No old employment-only SCM is run, read, or reported.
   Revision 14 Nov 2026 (pre-trend redesign):
   - Penalty lambda: UNCHANGED inner protocol (fit 2022, predict 2023).
   - Final weights: fitted on the FULL pre-event window Jan 2022-Aug 2024
     (32 months) with a free intercept per outcome, i.e. each outcome block is
     demeaned over the fit window so the weights match pre-event MOVEMENTS and
     the level difference is left to the DiD step (Arkhangelsky et al. 2021).
   - Backdated diagnostic weights: same objective and lambda fitted on
     2022-2023 only. They are used ONLY for the out-of-sample Jan-Aug 2024
     audit in stage 03, never for estimation.
   Nothing after August 2024 enters any fit. */
version 16.0
di as txt "02: ONE four-outcome joint synthetic-DiD weight vector (Jan 2022-Aug 2024 fit)"
use "${MC_OUT}/data/01_frozen_candidate_panel.dta", clear
sort mc_uid mdate
local maxit=${MC_MAX_ITER}
local tol=${MC_KKT_TOL}
mata:
real colvector mc_project_simplex(real colvector v)
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
real colvector mc_fista(real matrix A, real colvector y,
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
        wnew = mc_project_simplex(z-g/L)
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

// X: vertically stacked 57-month donor paths, 4 outcome blocks.
// Return [scaled X | scaled y]. All scales learned from 2022 ALONE.
// dm=1 removes each outcome block's own mean over months a..b from every
// column (donors and treated): equivalent to a free intercept per outcome.
real matrix mc_stack(real matrix X, real colvector y, real colvector s,
    real colvector which, real scalar a, real scalar b, real scalar dm)
{
    real scalar k, o, tt, n, no, denom
    real colvector src, dest
    real matrix Z
    n = cols(X)
    tt = b-a+1
    no = rows(which)
    Z = J(tt*no,n+1,.)
    for (k=1; k<=no; k++) {
        o = which[k]
        src = (o-1)*57 :+ (a::b)
        dest = ((k-1)*tt+1)::(k*tt)
        denom = s[o]*sqrt(tt*no)
        Z[dest,1..n] = X[src,.]/denom
        Z[dest,n+1] = y[src]/denom
        if (dm) Z[dest,.] = Z[dest,.] :- mean(Z[dest,.])
    }
    return(Z)
}

void mc_fit_joint_only()
{
    real matrix dat, XX, Z, X22, X23, Xfull, cv, fitdiag, paths
    real colvector group, donor, treated, step, y, scale, idx, ridx
    real colvector w, fit, lamgrid, scores, ds, sy, active, ddidx, which, obs
    real colvector wb, fitb, actb, wa, wba, icpt, icptb
    real matrix bddiag
    real scalar N, nn, n, nt, u, o, j, jb, k, q, r, best, chosen, selected
    real scalar mu, sd, floorv, maxit, tol, eff, lambda, jj
    dat = st_data(., tokens("mc_emp mc_lf mc_jobs mc_ur"))
    group = st_data(.,"mc_uid")
    donor = st_data(.,"mc_donor")
    treated = st_data(.,"treated_nc")
    step = st_data(.,"mc_step")
    N = rows(dat)
    if (mod(N,57)!=0) exit(459)
    nn=N/57
    n=sum((step:==32):&(donor:==1))
    nt=sum((step:==32):&(treated:==1))
    if (nt!=39 | n<500) {
        errprintf("Invalid frozen treated or eligible donor count.\n")
        exit(459)
    }
    XX=J(57*4,n,.)
    y=J(57*4,1,0)
    j=0
    for (u=1; u<=nn; u++) {
        ridx=(u-1)*57 :+ (1::57)
        if (sum(group[ridx]:==u)!=57 | sum(step[ridx]:==(1::57))!=57) exit(459)
        if (donor[ridx[1]]==1) j++
        for (o=1; o<=4; o++) {
            idx=(o-1)*57 :+ (1::57)
            if (donor[ridx[1]]==1) XX[idx,j]=dat[ridx,o]
            else if (treated[ridx[1]]==1) y[idx]=y[idx]+dat[ridx,o]
        }
    }
    if (j!=n) exit(459)
    y=y/39
    scale=J(4,1,.)
    for (o=1; o<=4; o++) {
        idx=(o-1)*57 :+ (1::12)
        ds=y[idx]
        mu=mean(ds)
        sd=sqrt(mean((ds:-mu):^2))
        floorv=0.25
        if (o==4) floorv=0.15
        scale[o]=max((floorv,sd))
    }
    st_matrix("MC_SCALES",scale')
    maxit=strtoreal(st_local("maxit"))
    tol=strtoreal(st_local("tol"))
    which=(1\2\3\4)
    Z=mc_stack(XX,y,scale,which,1,12,0)
    X22=Z[,1..n]
    sy=Z[,n+1]
    Z=mc_stack(XX,y,scale,which,13,24,0)
    X23=Z[,1..n]
    ds=Z[,n+1]
    lamgrid=(0.03\.1\.3\1\3\10)
    cv=J(rows(lamgrid),6,.)
    scores=J(rows(lamgrid),1,.)
    for (k=1; k<=rows(lamgrid); k++) {
        fit=mc_fista(X22,sy,lamgrid[k],maxit,tol)
        if (fit[n+3]!=0) {
            errprintf("Joint CV optimizer failed at lambda %g.\n",lamgrid[k])
            exit(459)
        }
        scores[k]=sqrt(mean((X23*fit[1::n]-ds):^2))
        cv[k,.]=(1,lamgrid[k],scores[k],fit[n+1],fit[n+2],fit[n+3])
    }
    best=min(scores)
    selected=1
    for (k=1; k<=rows(lamgrid); k++) {
        if (scores[k]<=1.01*best+1e-12) selected=k
    }
    lambda=lamgrid[selected]
    // Backdated diagnostic fit: 2022-2023 only, intercept, same lambda.
    Z=mc_stack(XX,y,scale,which,1,24,1)
    fitb=mc_fista(Z[,1..n],Z[,n+1],lambda,maxit,tol)
    if (fitb[n+3]!=0) {
        errprintf("Backdated optimizer did not converge; gap %g.\n",fitb[n+2])
        exit(459)
    }
    wb=fitb[1::n]
    bddiag=(fitb[n+1],fitb[n+2],sum(wb:>1e-8),1/sum(wb:^2))
    st_matrix("MC_BD",bddiag)
    // Definitive fit: full pre-event window (months 1-32), intercept.
    Z=mc_stack(XX,y,scale,which,1,32,1)
    Xfull=Z[,1..n]
    sy=Z[,n+1]
    fit=mc_fista(Xfull,sy,lambda,maxit,tol)
    if (fit[n+3]!=0) {
        errprintf("Joint final optimizer did not converge; gap %g.\n",fit[n+2])
        exit(459)
    }
    w=fit[1::n]
    eff=1/sum(w:^2)
    fitdiag=(1,lambda,scores[selected],fit[n+1],fit[n+2],sum(w:>1e-8),eff)
    ddidx=selectindex((donor:==1):&(step:==32))
    if (rows(ddidx)!=n) exit(459)
    j=st_addvar("double","mc_w_joint")
    st_store(ddidx,j,w)
    // Backdated weights are stored too, so the event study can be refitted
    // with a control that never saw Jan-Aug 2024 (stage 04d).
    jb=st_addvar("double","mc_w_joint_bd")
    st_store(ddidx,jb,wb)
    st_matrix("MC_CV",cv)
    st_matrix("MC_FITS",fitdiag)
    // Synthetic paths INCLUDE the fitted intercept: the definitive path is
    // shifted by its Jan 2022-Aug 2024 mean gap, the backdated path by its
    // 2022-2023 mean gap. Gaps therefore average zero over their fit window.
    paths=J(57,13,.)
    active=selectindex(w:>1e-10)
    actb=selectindex(wb:>1e-10)
    wa=w[active]/sum(w[active])
    wba=wb[actb]/sum(wb[actb])
    icpt=J(4,1,.)
    icptb=J(4,1,.)
    for (o=1; o<=4; o++) {
        idx=(o-1)*57 :+ (1::32)
        icpt[o]=mean(y[idx]-XX[idx,active]*wa)
        idx=(o-1)*57 :+ (1::24)
        icptb[o]=mean(y[idx]-XX[idx,actb]*wba)
    }
    st_matrix("MC_ICPT",(icpt',icptb'))
    for (r=1; r<=57; r++) {
        paths[r,1]=strtoreal(st_local("firstmonth"))+r-1
        for (o=1; o<=4; o++) {
            q=(o-1)*57+r
            jj=2+(o-1)*2
            paths[r,jj]=y[q]
            if (sum(XX[q,active]:>=.)==0 & y[q]<.) {
                paths[r,jj+1]=XX[q,active]*wa+icpt[o]
            }
            if (sum(XX[q,actb]:>=.)==0 & y[q]<.) {
                paths[r,9+o]=XX[q,actb]*wba+icptb[o]
            }
        }
    }
    st_matrix("MC_PATHS",paths)
}
end
local firstmonth=ym(2022,1)
mata: mc_fit_joint_only()

preserve
    keep if mc_donor==1 & mc_step==32
    keep fips mc_w_joint mc_w_joint_bd
    isid fips
    quietly summarize mc_w_joint, meanonly
    assert abs(r(sum)-1)<0.00001
    quietly summarize mc_w_joint_bd, meanonly
    assert abs(r(sum)-1)<0.00001
    assert inrange(mc_w_joint,0,1)
    save "${MC_OUT}/data/02_donor_weights.dta", replace
    export delimited using "${MC_OUT}/tables/donor_weights.csv", replace
restore
matrix colnames MC_CV = model_id lambda cv2023_rmse iterations dual_gap status
matrix colnames MC_FITS = model_id lambda cv2023_rmse iterations dual_gap n_active effective_donors
matrix colnames MC_SCALES = scale_emp scale_lf scale_jobs scale_urate
matrix colnames MC_PATHS = mdate mc_y_emp mc_syn_joint_emp ///
    mc_y_lf mc_syn_joint_lf mc_y_jobs mc_syn_joint_jobs ///
    mc_y_urate mc_syn_joint_urate ///
    mc_syn_bd_emp mc_syn_bd_lf mc_syn_bd_jobs mc_syn_bd_urate
matrix colnames MC_BD = bd_iterations bd_dual_gap bd_n_active bd_effective_donors
matrix colnames MC_ICPT = icpt_emp icpt_lf icpt_jobs icpt_urate ///
    bd_icpt_emp bd_icpt_lf bd_icpt_jobs bd_icpt_urate
preserve
    clear
    svmat double MC_CV, names(col)
    gen str16 model="joint"
    assert status==0
    save "${MC_OUT}/data/02_inner_cv.dta", replace
    export delimited using "${MC_OUT}/tables/inner_2023_cv.csv", replace
restore
preserve
    clear
    svmat double MC_FITS, names(col)
    svmat double MC_BD, names(col)
    svmat double MC_ICPT, names(col)
    gen str16 model="joint"
    gen byte fit_first_month=1
    gen byte fit_last_month=32
    assert _N==1
    assert dual_gap<=`tol'
    assert bd_dual_gap<=`tol'
    save "${MC_OUT}/data/02_solver_audit.dta", replace
    export delimited using "${MC_OUT}/tables/solver_audit.csv", replace
restore
preserve
    clear
    svmat double MC_SCALES, names(col)
    export delimited using "${MC_OUT}/tables/training_2022_scales.csv", replace
restore
preserve
    clear
    svmat double MC_PATHS, names(col)
    format mdate %tm
    isid mdate
    save "${MC_OUT}/data/02_paths.dta", replace
    export delimited using "${MC_OUT}/tables/monthly_paths.csv", replace
restore
* One common fitted weight for ALL six matched outcomes, plus the backdated
* vector (2022-23 fit) carried alongside it for the stage-04d robustness refit.
drop mc_w_joint mc_w_joint_bd
merge m:1 fips using "${MC_OUT}/data/02_donor_weights.dta", ///
    keep(master match) nogen
replace mc_w_joint=0 if treated_nc==1
replace mc_w_joint_bd=0 if treated_nc==1
assert !missing(mc_w_joint) if mc_donor==1
assert !missing(mc_w_joint_bd) if mc_donor==1
isid fips mdate
save "${MC_OUT}/data/02_weighted_county_month.dta", replace
di as result "02 complete: one four-outcome synthetic-DiD weight vector (full pre-event fit, intercept)."
