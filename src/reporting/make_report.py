#!/usr/bin/env python3
"""Rebuild readable results and a complete figure atlas from passed stages only."""
from pathlib import Path
import argparse, csv, hashlib, html, json, shutil, tempfile, textwrap
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[2]
NAMES = {'emp':'Resident employment','lf':'Labor force','jobs':'Workplace jobs','urate':'Unemployment rate','diff':'Jobs / resident employment','share':'Employment / labor force','unemp':'Unemployed residents'}
WINDOWS = {'placebo':'Jan–Aug 2024 (pre-event placebo)','immediate':'Oct–Dec 2024','year2025':'Jan–Sep 2025','late':'Nov 2025–Mar 2026','short':'Jan–Sep 2025','medium':'Nov 2025–Sep 2026*','late_extended':'Nov 2025–Sep 2026','aprsep2026':'Apr–Sep 2026'}
BLUE,ORANGE,GRAY='#176B91','#C46822','#64748B'
def key(s): return str(s).replace('mc_y','').replace('mx_','')
def label(s): return NAMES.get(key(s),str(s))
def units(s): return 'Percentage points' if key(s)=='urate' else '100 × log change (≈ %)'
def fmt(x):
    return 'NA' if pd.isna(x) else (f'{float(x):,.0f}' if abs(float(x)) >= 1 and float(x).is_integer() else f'{float(x):.3g}')
def month(s):
    y,m=str(s).split('m'); return pd.Timestamp(int(y),int(m),1)
def main(root=ROOT):
    reports=root/'outputs/reports'; reports.mkdir(parents=True,exist_ok=True)
    status=pd.read_csv(root/'outputs/audit/STAGE_STATUS.csv')
    passed=set(status.loc[status.status=='PASSED','stage'])
    inputs=[]; unavailable=[]
    def read(rel,stage):
        p=root/rel
        if stage not in passed:
            unavailable.append(f'{rel}: stage {stage} not PASSED'); return pd.DataFrame()
        if not p.exists():
            unavailable.append(f'{rel}: missing file'); return pd.DataFrame()
        inputs.append({'path':rel,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'stage':stage})
        return pd.read_csv(p,na_values=['.'],dtype={'fips':str})
    def table(n): return read('outputs/labor/tables/'+n+'.csv','joint_multioutcome')
    effects=table('period_effects_WITH_HOLM'); paths=table('monthly_paths'); fit=table('pretrend_train_holdout_ALL')
    trends=table('twfe_pretrend_tests'); events=table('event_coefficients'); weights=table('donor_weights')
    cv=table('inner_2023_cv'); solver=table('solver_audit'); scales=table('training_2022_scales'); times=table('pre2024_time_placebos')
    ext=table('LAUS_extended_period_effects_EXPLORATORY'); extev=table('LAUS_extended_event_coefficients_EXPLORATORY'); extpre=table('LAUS_extended_pretrends')
    curves=table('gaussian_rebased_curves'); county=table('kernel_county_windows'); slopes=table('kernel_dose_slopes')
    draws=read('outputs/labor/placebos/joint_placebo_draws.csv','joint_placebos')
    ranks=read('outputs/labor/placebos/joint_placebo_summary.csv','joint_placebos')
    migration=read('outputs/migration/census_annual_effects.csv','census_annual_migration')
    groups=read('outputs/migration/census_group_comparison.csv','census_annual_migration')
    coverage=read('outputs/preprocessing/audits/final_panel_coverage_by_month.csv','preinference')
    text=[]
    def add(*s): text.extend(str(x) for x in s)
    def section(s): add('',s.upper(),'='*88)
    def tab(df,cols,renames=None):
        if df.empty: add('Unavailable for this run.'); return
        d=df[cols].rename(columns=renames or {}).copy()
        for c in d:
            if pd.api.types.is_numeric_dtype(d[c]): d[c]=d[c].map(fmt)
        add(d.to_string(index=False))
    add('HURRICANE HELENE — RESULTS AND INTERPRETATION',f'Report generated: {datetime.now().astimezone().isoformat(timespec="seconds")}',
        'Based on the latest recorded stage statuses. Reporting does not refit models.')
    section('Read this first')
    if not effects.empty:
        n=int((effects.holm_p<.05).sum())
        add(f'{n} of {len(effects)} main outcome/window tests are significant at 5% after Holm adjustment.',
            'No significant result does not prove no damage, a small effect, or full recovery.',
            'First full post-event window (October–December 2024), relative to the pre-event mean (Jan 2022–Aug 2024):')
        for _,r in effects[effects.window=='immediate'].iterrows():
            add(f"  {label(r.outcome)}: {r.display_beta:+.2f} {'pp' if key(r.outcome)=='urate' else 'log-%'}; 95% CI [{r.display_lo:.2f}, {r.display_hi:.2f}]; Holm p={r.holm_p:.3f}.")
    if not trends.empty:
        k=len(trends)
        add(f'PRE-EVENT CHECK (weights are fitted on this window, so it is in-sample): linear pre-trend test rejects at 5% for {int((trends.slope_p<.05).sum())}/{k} outcomes; '
            f'half-year-bin test for {int((trends.bin_p<.05).sum())}/{k}; {int(trends.n_sig05.sum())} of {32*k} monthly pre-event coefficients are individually significant.',
            f'Stricter joint tests: last eight pre-event months reject for {int((trends.last8_p<.05).sum())}/{k} outcomes; all 31 monthly coefficients for {int((trends.pretrend_p<.05).sum())}/{k}.',
            'Surface parallelism (no drift, no visible gap) is not proof of parallel trends: the strict joint tests still detect small month-specific deviations.')
    if not fit.empty:
        add(f'OUT-OF-SAMPLE CHECK (weights refitted on 2022–23 only, then predicting Jan–Aug 2024): error is more than twice the training error for {int((fit.bd_rmspe_ratio>2).sum())}/{len(fit)} outcomes.')
    add('Placebo ranks are descriptive comparisons, not randomized causal p-values.',
        'Migration is a separate, unweighted eligible-county comparison, not the labor synthetic control.',
        'Graphs: open HELENE_GRAPH_ATLAS.pdf or graphs/index.html beside this report.')
    section('Study design and how to read the numbers')
    add('39 frozen North Carolina counties; equal county weighting within the treated group.',
        'Synthetic difference-in-differences: one shared donor-weight vector fitted jointly to four labor outcomes',
        'over the full pre-event window (January 2022–August 2024) with a free intercept per outcome.',
        'Penalty chosen by the unchanged inner protocol (fit 2022, predict 2023). Nothing after August 2024 enters the fit.',
        'Effects and event-study coefficients are relative to the pre-event mean gap (uniform time weights), not to August 2024 alone.',
        'Because 2024 is now inside the fit window there is no untouched holdout for the definitive weights;',
        'the out-of-sample evidence is the backdated refit (2022–23 only) reported under Validation.',
        'September 2024 is a transition month and is excluded from the event regressions.',
        'October 2025 LAUS observations are missing; they are not treated as zero.',
        'Main windows: Oct–Dec 2024; Jan–Sep 2025; Nov 2025–Mar 2026.',
        'Log-% means 100 × log change, approximately a percent change for small values.',
        'pp means percentage points. Jobs/employment and employment/labor-force are log ratios.',
        'Negative employment/jobs/labor-force estimates mean a relative decline; positive unemployment-rate estimates mean a relative increase.',
        '95% intervals are individual, county-clustered intervals, not simultaneous Holm-adjusted intervals.',
        'The design was frozen after earlier outcomes had been inspected; this is not preregistration.')
    section('Execution and data coverage')
    tab(status,['stage','status','return_code'])
    add('PASSED means the stage executed successfully, not that its identifying assumptions hold.')
    if not coverage.empty:
        add(f"Calendar: {coverage.mdate.iloc[0]} to {coverage.mdate.iloc[-1]} ({len(coverage)} months). See coverage graph for missing outcome months.")
    section('Main effects — all 18 planned tests')
    if not effects.empty:
        d=effects.copy();d['Outcome']=d.outcome.map(label);d['Period']=d.window.map(WINDOWS)
        tab(d,['Outcome','Period','display_beta','display_lo','display_hi','pvalue','holm_p'],{'display_beta':'Estimate','display_lo':'CI low','display_hi':'CI high','pvalue':'Raw p','holm_p':'Holm p'})
    else: add('Unavailable for this run.')
    section('Validation — pre-event fit and what it does not show')
    if not fit.empty:
        d=fit.copy();d['Outcome']=d.outcome.map(label)
        add('(a) Definitive weights, in-sample pre-event gap (log-%; unemployment rate in pp):')
        tab(d,['Outcome','fit_rmse','fit_rmse_2022_23','fit_rmse_2024','fit_max_abs_gap'],{'fit_rmse':'RMSE Jan22–Aug24','fit_rmse_2022_23':'RMSE 2022–23','fit_rmse_2024':'RMSE Jan–Aug24','fit_max_abs_gap':'Largest monthly gap'})
        add('(b) Backdated weights (fitted on 2022–23 only), out-of-sample Jan–Aug 2024:')
        tab(d,['Outcome','bd_train_rmse','bd_holdout_rmse','bd_rmspe_ratio','bd_aug2024_gap','bd_holdout_slope_p'],{'bd_train_rmse':'Train RMSE','bd_holdout_rmse':'Holdout RMSE','bd_rmspe_ratio':'Ratio','bd_aug2024_gap':'Aug 2024 gap','bd_holdout_slope_p':'Holdout slope p'})
        add('RMSE is average mismatch in outcome units; smaller is better.',
            'Table (a) is flat partly by construction. Table (b) is the honest forecast test: a ratio well above 1 means',
            'donors chosen on 2022–23 did not predict the treated counties in 2024, so the flat fit in (a) relies on 2024 being in the fit window.')
    if not trends.empty:
        d=trends.copy();d['Outcome']=d.outcome.map(label)
        f=np.where(d.outcome.map(key)=='urate',1,100);d['pre_rmse']=d.pre_rmse*f;d['max_abs_pre']=d.max_abs_pre*f
        add('(c) Event-study pre-event tests, definitive weights (county-clustered):')
        tab(d,['Outcome','slope_p','bin_p','n_sig05','max_abs_t','pre_rmse','max_abs_pre','last8_p','pretrend_p'],{'slope_p':'Linear trend p','bin_p':'Half-year bins p','n_sig05':'Months p<.05 (of 32)','max_abs_t':'Largest |t|','pre_rmse':'Pre RMSE','max_abs_pre':'Largest pre coef','last8_p':'Jan–Jul24 vs Aug24 p','pretrend_p':'All 31 months p'})
        add('Linear-trend and bin tests ask whether the pre-event gap drifts; they have limited power with 39 treated counties.',
            'The last two columns are strict joint tests of exact month-by-month equality; small p-values there mean small but',
            'precisely estimated month-specific deviations remain. Compare their size (Pre RMSE, largest coefficient) with the effects.')
    section('Donor model and tuning')
    tab(solver,['lambda','cv2023_rmse','iterations','dual_gap','n_active','effective_donors','bd_n_active','bd_effective_donors'])
    add('cv2023_rmse: penalty selection only (fit 2022, predict 2023). bd_ columns: backdated diagnostic weights.')
    tab(cv,['lambda','cv2023_rmse','dual_gap','status'])
    tab(scales,list(scales.columns))
    if not weights.empty:
        tab(weights.sort_values('mc_w_joint',ascending=False).head(15),['fips','mc_w_joint'],{'fips':'County FIPS','mc_w_joint':'Weight'})
        add('Top 15 shown here; the atlas and donor_weights.csv retain the full distribution.')
    section('Geographic placebo benchmarks')
    if not draws.empty: add(f"{len(draws)} refits; {int(draws.converged.sum())} converged. Same frozen lambda and 39-county pseudo-treated group size.")
    if not ranks.empty:
        d=ranks.copy();d['Outcome']=d.outcome.map(label)
        tab(d,['Outcome','screen','observed_shock','n_placebos','left_rank','right_rank'],{'right_rank':'Two-sided rank','left_rank':'Lower-tail rank'})
    add('Observed shock = Oct–Dec 2024 average gap minus the pre-event (Jan 2022–Aug 2024) mean gap; placebo groups are fitted the same way.',
        'Ranks use a +1 correction. Fit screens compare joint pre-fit errors; they do not establish exchangeable treatment assignment.')
    section('Exploratory labor extension — through September 2026')
    if not ext.empty:
        d=ext.copy(); factor=np.where(d.outcome.map(key)=='urate',1,100)
        for c in ['beta','lower95','upper95']:d[c]=d[c]*factor
        d['Outcome']=d.outcome.map(label);d['Period']=d.window.map(WINDOWS)
        tab(d,['Outcome','Period','beta','lower95','upper95','pvalue'])
    if not extpre.empty:
        d=extpre.copy();d['Outcome']=d.outcome.map(label);tab(d,['Outcome','slope_p','last8_p','pretrend_p'])
    add('Exploratory p-values are unadjusted and outside the 18-test main Holm family.')
    section('Exploratory migration and population')
    tab(migration,['outcome','b2025','lo95','hi95','p2025','pretest_p','N_counties'])
    add('netdom_per1000: net domestic migrants per 1,000 residents. pop_growth_pct: population growth in percent; effect is a percentage-point difference.',
        '2025 is July 2024–July 2025, only partly after Helene. Reference is the annual period ending July 2024.',
        'A negative net-migration contrast can mean fewer arrivals, more departures, or both; it does not count displaced original residents.')
    tab(groups,['year','grp','counties','mean_pop_growth_pct','mean_netdom_per1000'])
    add('Group 1 = treated; group 0 = eligible comparison counties. These means are unweighted across counties.')
    section('Exploratory rainfall and damage curves')
    add('The atlas includes every exported outcome × window × exposure curve, with county points and 95% pointwise bootstrap bands.',
        'These are within-treated-county associations, not causal dose-response estimates.',
        'Immediate = Oct–Dec 2024; short = Jan–Sep 2025; medium = Nov 2025–Sep 2026, except workplace jobs end in March 2026.',
        'HDSI is a composite damage index; rainfall is measured in millimetres. The shared synthetic baseline shifts levels, not exposure slopes.',
        'Each county is measured against its own pre-event mean gap (Jan 2022–Aug 2024), not against August 2024 alone.',
        'Placebo window: Jan–Aug 2024 against the county 2022–23 mean gap. A flat placebo curve is the pre-trend check for these curves.')
    if not slopes.empty:
        d=slopes.copy();d['Outcome']=d.outcome.map(label);d['Window']=d.window.map(WINDOWS)
        add('Linear exposure gradient per one SD of exposure (HC3 standard errors):')
        tab(d,['dose','Window','Outcome','n_counties','mean_change','slope_per_sd','se','pvalue'])
        pl=d[d.window=='placebo']
        add(f'Placebo gradients significant at 5%: {int((pl.pvalue<.05).sum())} of {len(pl)}. Compare their size with the post-event gradients.')
    if not curves.empty:
        d=curves.groupby(['dose','window','outcome'],sort=False).agg(counties=('n_counties','first'),bootstrap_draws=('n_boot','first'),min_effective_n=('effective_n','min'),max_effective_n=('effective_n','max')).reset_index()
        tab(d,list(d.columns));add('Low effective sample sizes signal weak local support; bands are pointwise and not corrected for searching many curves.')
    section('What these results can and cannot support')
    add('The fitted comparisons describe relative post-event changes. The pre-event gap is now flat on the surface, but it is flat in-sample;',
        'the backdated out-of-sample check and the strict joint tests above remain the relevant cautions for a causal reading.',
        'Failure to reject is not evidence of no effect. Small placebo ranks do not override a failed out-of-sample check.',
        'County-clustered standard errors omit some regional dependence and donor-weight estimation uncertainty.',
        'Resident employment (LAUS) and workplace jobs (QCEW) describe different populations.',
        'All source tables remain separate CSV files. No raw county-by-year dump is embedded in this readable report.')
    if unavailable: section('Unavailable outputs'); add(*unavailable)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,'axes.prop_cycle':plt.cycler(color=[BLUE,ORANGE,GRAY]),'savefig.facecolor':'white'})
    staging=Path(tempfile.mkdtemp(prefix='helene-report-',dir=reports)); graph=staging/'graphs';graph.mkdir(); (graph/'.gitkeep').touch()
    catalog=[]
    with PdfPages(staging/'HELENE_GRAPH_ATLAS.pdf') as pdf:
        def figgrid(title,subtitle,rows=2,cols=2):
            fig,axs=plt.subplots(rows,cols,figsize=(14,10),squeeze=False)
            fig.suptitle(title,x=.07,ha='left',fontsize=19,fontweight='bold',y=.975)
            fig.text(.07,.925,subtitle,fontsize=10,color=GRAY)
            fig.subplots_adjust(left=.17,right=.96,bottom=.14,top=.83,hspace=.85,wspace=.65)
            return fig,axs.ravel()
        def save(fig,slug,caption):
            fig.text(.07,.025,textwrap.fill(caption,145),fontsize=9,color=GRAY)
            fig.savefig(graph/(slug+'.png'),dpi=170);pdf.savefig(fig);plt.close(fig)
            catalog.append({'file':slug+'.png','title':fig._suptitle.get_text(),'caption':caption})
        def zero(ax): ax.axhline(0,color=GRAY,lw=.8);ax.grid(axis='y',alpha=.15)
        def dateaxis(ax):
            import matplotlib.dates as md
            ax.xaxis.set_major_locator(md.MonthLocator(interval=12));ax.xaxis.set_major_formatter(md.DateFormatter('%b %Y'))
        if not effects.empty:
            fig,axes=figgrid('Main effects: six outcomes, three windows','Points and individual 95% confidence intervals; Holm correction applies to the 18 tests.',3,2)
            for ax,o in zip(axes,['emp','lf','jobs','urate','diff','share']):
                d=effects[effects.outcome.map(key)==o].set_index('window').reindex(['immediate','year2025','late'])
                y=np.arange(3);b=d.display_beta.to_numpy();ax.errorbar(b,y,xerr=[b-d.display_lo,d.display_hi-b],fmt='o',color=BLUE,capsize=3)
                ax.axvline(0,color=GRAY,lw=.8);ax.set_yticks(y,[WINDOWS[w] for w in d.index],fontsize=8);ax.invert_yaxis();ax.set(title=NAMES[o],xlabel=units(o))
            save(fig,'01_main_effects','Intervals are not multiplicity-adjusted. See the readable report for every raw and Holm-adjusted p-value.')
        if not paths.empty:
            dates=paths.mdate.map(month)
            fig,axes=figgrid('Treated and synthetic-control paths','Four outcomes, centered on 2022 means; synthetic path includes the fitted intercept. Shading: 2024 months, now inside the fit window.')
            for ax,o in zip(axes,['emp','lf','jobs','urate']):
                ax.plot(dates,paths['mc_y_'+o],label='Treated',lw=1.6);ax.plot(dates,paths['mc_syn_joint_'+o],label='Synthetic control',ls='--',lw=1.5)
                pre=dates<pd.Timestamp('2024-09-01');ax.plot(dates[pre],paths.loc[pre,'mc_syn_bd_'+o],label='Backdated (2022–23 fit)',ls=':',lw=1.1,color=GRAY)
                ax.axvspan(pd.Timestamp('2024-01-01'),pd.Timestamp('2024-08-31'),color=GRAY,alpha=.12);ax.axvline(pd.Timestamp('2024-09-01'),color=GRAY,ls=':');ax.set(title=NAMES[o],ylabel=units(o));dateaxis(ax)
            axes[0].legend(frameon=False,fontsize=8)
            save(fig,'02_monthly_paths','Missing observations remain gaps. The vertical dotted line marks September 2024.')
        for ev,title,slug in [(events,'Main monthly event-study estimates','03_main_events'),(extev,'Exploratory monthly labor extension','10_extended_events')]:
            if ev.empty: continue
            outcomes=list(ev.outcome.unique());fig,axes=figgrid(title,'Relative to the pre-event mean (Jan 2022–Aug 2024). Points and individual 95% confidence intervals.',3,2)
            for ax,o in zip(axes,outcomes):
                d=ev[ev.outcome==o].sort_values('event_month');f=1 if key(o)=='urate' else 100
                ax.errorbar(d.event_month,d.beta*f,yerr=[(d.beta-d.lower95)*f,(d.upper95-d.beta)*f],fmt='.',ms=4,elinewidth=.7,color=BLUE,alpha=.85)
                zero(ax);ax.axvline(0,color=GRAY,ls=':');ax.set(title=label(o),ylabel=units(o),xlabel='Months from September 2024')
            for ax in axes[len(outcomes):]:ax.set_visible(False)
            save(fig,slug,'No point is imputed for excluded months. Pre-event coefficients average zero by construction. Exploratory extensions are not Holm-adjusted.')
        if not fit.empty:
            fig,axes=figgrid('Validation before interpreting effects','In-sample pre-event fit, backdated out-of-sample forecast, and formal pre-event tests answer different questions.')
            for ax,subset in zip(axes[:2],[fit[fit.outcome!='urate'],fit[fit.outcome=='urate']]):
                x=np.arange(len(subset));ax.bar(x-.27,subset.fit_rmse_2024,.27,label='Definitive weights, Jan–Aug 2024 (in-sample)');ax.bar(x,subset.bd_train_rmse,.27,label='Backdated weights, 2022–23 (training)',color=GRAY);ax.bar(x+.27,subset.bd_holdout_rmse,.27,label='Backdated weights, Jan–Aug 2024 (out-of-sample)',hatch='//',color=ORANGE);ax.set_xticks(x,subset.outcome.map(label),fontsize=8);ax.set(ylabel='RMSE: '+('pp' if len(subset)==1 else 'log-%'));ax.legend(fontsize=7,frameon=False)
            for ax,cols,ttl in zip(axes[2:],[['slope_p','bin_p'],['last8_p','pretrend_p']],['Drift tests: linear trend (●) and half-year bins (△)','Strict joint tests: Jan–Jul 2024 (●) and all 31 months (△)']):
                if trends.empty:ax.set_visible(False);continue
                y=np.arange(len(trends))
                for col,mk in zip(cols,['o','^']):ax.scatter(trends[col].clip(lower=1e-20),y,color=BLUE,marker=mk,facecolors=BLUE if mk=='o' else 'none')
                ax.set_xscale('log');ax.axvline(.05,color=ORANGE,ls='--',label='5% threshold');ax.set_yticks(y,trends.outcome.map(label),fontsize=8);ax.set(xlabel='p-value (log scale)',title=ttl);ax.title.set_fontsize(10);ax.legend(fontsize=8,frameon=False)
            save(fig,'04_validation','Smaller p-values indicate pre-event differences. The definitive weights are fitted on the whole pre-event window, so their flat fit is not an out-of-sample validation; the hatched bars are.')
        if not extpre.empty:
            fig,axes=figgrid('Exploratory extension: pre-event diagnostics','Formal pre-event tests for the longer LAUS estimation sample.',1,2)
            for ax,col in zip(axes,['pretrend_p','last8_p']):
                ax.scatter(extpre[col],np.arange(len(extpre)),color=BLUE);ax.set_xscale('log');ax.axvline(.05,color=ORANGE,ls='--',label='5% threshold');ax.set_yticks(np.arange(len(extpre)),extpre.outcome.map(label),fontsize=9);ax.set(xlabel='p-value (log scale)',title=col.replace('_',' '));ax.legend(frameon=False,fontsize=8)
            save(fig,'04b_extended_validation','Small p-values indicate pre-event differences. These checks remain exploratory and unadjusted.')
        if not weights.empty:
            fig,axes=figgrid('Donor weights and model tuning','One shared weight vector fitted on Jan 2022–Aug 2024; the penalty is tuned on 2022–2023 only.')
            w=weights.sort_values('mc_w_joint',ascending=False);top=w.head(15).iloc[::-1]
            axes[0].barh(top.fips,100*top.mc_w_joint);axes[0].set(xlabel='Weight (%)',title='Top 15 donors (county FIPS)');axes[0].tick_params(axis='y',labelsize=8)
            axes[1].plot(np.arange(1,len(w)+1),100*w.mc_w_joint.cumsum());axes[1].set(xlabel='Donors ranked by weight',ylabel='Cumulative weight (%)',ylim=(0,102))
            if not cv.empty:axes[2].plot(cv['lambda'],cv.cv2023_rmse,'o-');axes[2].set(xscale='log',xlabel='Regularization parameter',ylabel='2023 validation error')
            axes[3].axis('off')
            if not solver.empty:
                s=solver.iloc[0];axes[3].text(0,.85,f"Selected lambda: {fmt(s['lambda'])}\nActive donors: {fmt(s.n_active)}\nEffective donors: {fmt(s.effective_donors)}\nSolver dual gap: {fmt(s.dual_gap)}",va='top',fontsize=14,linespacing=1.8)
            save(fig,'05_donors','County identifiers are FIPS codes. The effective donor count reflects weight concentration.')
        if not draws.empty and not ranks.empty:
            fig,axes=figgrid('Geographic placebo distributions','Same four-outcome fitting procedure; 39 counties in each pseudo-treated group.')
            for ax,o in zip(axes,['emp','lf','jobs','urate']):
                ax.hist(draws.loc[draws.converged==1,'shock_'+o],bins=24,color=BLUE,alpha=.75)
                r=ranks[(ranks.outcome==o)&(ranks.screen=='all')].iloc[0];ax.axvline(r.observed_shock,color=ORANGE,lw=2,label=f'Observed {r.observed_shock:+.2f}');ax.set(title=NAMES[o],xlabel=units(o),ylabel='Placebo draws');ax.legend(fontsize=8,frameon=False)
            save(fig,'06_placebo_distributions','Ranks are descriptive benchmarks, not randomized causal p-values. All converged draws are shown.')
            fig,axes=figgrid('Placebo fit and sensitivity to fit screens','Convergence and pre-event fit are reported separately from the post-event benchmark.')
            axes[0].scatter(draws.prefit_rmse,draws.pre2024_rmse,s=12,alpha=.6);axes[0].axvline(draws.observed_joint_prefit.iloc[0],ls='--',color=ORANGE);axes[0].set(xlabel='Scaled pre-event RMSE (Jan 2022–Aug 2024)',ylabel='Scaled RMSE, Jan–Aug 2024 segment')
            axes[1].scatter(draws.iterations,draws.dual_gap,s=12);axes[1].axhline(draws.kkt_tolerance.iloc[0],color=ORANGE,ls='--');axes[1].set(xlabel='Iterations',ylabel='Solver dual gap')
            for ax,col in zip(axes[2:],['left_rank','right_rank']):
                for o in ranks.outcome.unique():
                    d=ranks[ranks.outcome==o];ax.plot(np.arange(len(d)),d[col],marker='o',label=NAMES[o])
                ax.set_xticks(np.arange(len(d)),d.screen,rotation=15,fontsize=8);ax.set(ylabel='Lower-tail rank' if col=='left_rank' else 'Two-sided rank',ylim=(0,1.05));ax.legend(fontsize=7,frameon=False)
            save(fig,'07_placebo_diagnostics','Dashed lines: observed joint pre-event fit error (left) and frozen convergence tolerance (right). Both axes of the first panel are in-sample.')
        if not times.empty:
            fig,axes=figgrid('Pre-event time-placebo contrasts (backdated weights)','Weights fitted on 2022–23 only, so every 2024 contrast is out-of-sample. Compare earlier gap changes with the Oct–Dec 2024 shock.')
            for ax,o in zip(axes,['emp','lf','jobs','urate']):
                d=times[times.outcome==o];ax.bar(np.arange(len(d)),d.reference_change,color=[BLUE]*4+[ORANGE]);zero(ax);ax.set_xticks(np.arange(len(d)),['Mar–May vs Jan','Apr–Jun vs Feb','May–Jul vs Mar','Jun–Aug vs Apr','Helene']);ax.tick_params(axis='x',labelrotation=20,labelsize=8);ax.set(title=NAMES[o],ylabel=units(o))
            save(fig,'08_time_placebos','Earlier checks are all in 2024: three-month average gap minus the indicated reference month; Helene bar is Oct–Dec 2024 minus August 2024 for the same backdated weights. No intervals are provided.')
        if not ext.empty:
            fig,axes=figgrid('Exploratory period effects through September 2026','Unadjusted individual 95% intervals; separate from the 18 main tests.',3,2)
            outcomes=list(ext.outcome.unique())
            for ax,o in zip(axes,outcomes):
                d=ext[ext.outcome==o];f=1 if key(o)=='urate' else 100;b=d.beta*f;y=np.arange(len(d));ax.errorbar(b,y,xerr=[b-d.lower95*f,d.upper95*f-b],fmt='o',capsize=3);ax.axvline(0,color=GRAY,lw=.8);ax.set_yticks(y,d.window.map(WINDOWS),fontsize=8);ax.set(title=label(o),xlabel=units(o))
            for ax in axes[len(outcomes):]:ax.set_visible(False)
            save(fig,'09_extended_periods','Do not interpret failure to reject as proof of full recovery.')
        if not groups.empty and not migration.empty:
            fig,axes=figgrid('Exploratory annual migration and population','Unweighted eligible-county comparison; 2025 spans July 2024–July 2025.')
            for ax,col,name in zip(axes[:2],['mean_netdom_per1000','mean_pop_growth_pct'],['Net domestic migration per 1,000','Population growth (%)']):
                for g,sty in [(1,'-'),(0,'--')]:
                    d=groups[(groups.grp==g)&(groups.year>=2021)];ax.plot(d.year,d[col],sty,marker='o',label='Treated' if g else 'Eligible controls')
                ax.set(ylabel=name,xticks=range(2021,2026));ax.legend(fontsize=8,frameon=False)
            for ax,(_,r) in zip(axes[2:],migration.iterrows()):
                ax.errorbar([r.b2025],[0],xerr=[[r.b2025-r.lo95],[r.hi95-r.b2025]],fmt='o',capsize=5);ax.axvline(0,color=GRAY);ax.set(yticks=[],title='2025 contrast: '+('net migration' if r.outcome=='netdom_per1000' else 'population growth'),xlabel='Migrants per 1,000' if r.outcome=='netdom_per1000' else 'Percentage points')
            save(fig,'11_migration','Net migration cannot distinguish fewer arrivals from more departures. Individual 95% intervals; exploratory, unadjusted.')
        if not curves.empty:
            for (dose,window),d in curves.groupby(['dose','window'],sort=False):
                fig,axes=figgrid('Exploratory exposure curves: '+('damage index' if dose=='hdsi' else 'rainfall'),WINDOWS[window]+' — within-treated-county associations; pointwise bootstrap bands.',3,2)
                for ax,o in zip(axes,['emp','lf','jobs','urate','unemp']):
                    sub=d[d.outcome==o].sort_values('intensity');dots=county[(county.outcome==o)&(county.window==window)&(county.complete==1)] if not county.empty else pd.DataFrame()
                    if not dots.empty:ax.scatter(dots[dose],dots.change,s=12,color=GRAY,alpha=.4)
                    ax.plot(sub.intensity,sub.fitted,color=BLUE);ax.fill_between(sub.intensity,sub.lower,sub.upper,color=BLUE,alpha=.18);zero(ax);ax.set(title=NAMES[o],xlabel='HDSI (index)' if dose=='hdsi' else 'Rainfall (mm)',ylabel=units(o))
                axes[5].axis('off');axes[5].text(0,.8,'Dots: complete-window county changes\nrelative to the county pre-event mean gap\nLine: local-linear Gaussian fit\nBand: 95% pointwise bootstrap interval\n\nExploratory association, not a causal curve.',fontsize=10,va='top',linespacing=1.8)
                save(fig,'12_kernel_'+dose+'_'+window,'* In the medium window, workplace jobs end in March 2026; other outcomes end in September 2026. Low local support increases uncertainty.')
        if not coverage.empty:
            fig,axes=figgrid('Source coverage by month','Observed county counts; missingness remains explicit.',1,1);ax=axes[0]
            for c,name,sty in [('n_laus','LAUS','-'),('n_qcew','QCEW','--'),('n_pool','Eligible donor pool',':')]:ax.plot(coverage.mdate.map(month),coverage[c],sty,label=name)
            ax.set(ylabel='Counties',ylim=(0,max(coverage.n_laus.max(),coverage.n_qcew.max())*1.08));ax.legend(frameon=False);dateaxis(ax)
            save(fig,'13_coverage','Zero coverage marks missing source observations, not zero employment. QCEW ends in March 2026.')
        if not catalog:
            fig,axes=figgrid('No verified results available','No numerical figures were generated.');[ax.axis('off') for ax in axes];save(fig,'00_unavailable','See STAGE_STATUS.csv and the readable report.')
    add('',f'Figure atlas: {len(catalog)} pages. Every graph is also available as a PNG in graphs/.',
        'Input hashes and figure inventory: graphs/manifest.json.')
    (staging/'HELENE_COMPLETE_RESULTS.txt').write_text('\n'.join(text)+'\n')
    body=''.join(f'<section id="g{i}"><h2>{html.escape(c["title"])}</h2><img src="{c["file"]}" alt="{html.escape(c["title"]+". "+c["caption"])}"><p>{html.escape(c["caption"])}</p></section>' for i,c in enumerate(catalog))
    nav=''.join(f'<li><a href="#g{i}">{html.escape(c["title"])}</a></li>' for i,c in enumerate(catalog))
    (graph/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Helene results atlas</title><style>body{font:17px/1.6 system-ui;color:#203040;max-width:1200px;margin:40px auto;padding:0 24px;background:#f5f7fa}h1{font-size:40px}section{background:white;padding:24px;margin:30px 0;border-radius:12px}img{width:100%;height:auto}a{color:#176B91}p{max-width:85ch}</style><h1>Hurricane Helene: results atlas</h1><p>Complete current-run figures. The pre-event fit is in-sample; these charts are not a certificate of causality.</p><p><a href="../HELENE_COMPLETE_RESULTS.txt">Readable results</a> · <a href="../HELENE_GRAPH_ATLAS.pdf">Download the full PDF</a></p><ol>'+nav+'</ol>'+body+'</html>')
    (graph/'manifest.json').write_text(json.dumps({'generated':datetime.now().astimezone().isoformat(),'figures':catalog,'inputs':inputs,'unavailable':unavailable},indent=2))
    # Publish only after all figures and the report were built successfully.
    dest=reports/'graphs'
    if dest.exists(): shutil.rmtree(dest)
    shutil.move(str(graph),str(dest))
    for n in ['HELENE_COMPLETE_RESULTS.txt','HELENE_GRAPH_ATLAS.pdf']: (staging/n).replace(reports/n)
    staging.rmdir()
    (root/'outputs/audit/REPORT_BUILD_OK.txt').write_text(f'SUCCESS: {len(catalog)} figures; {len(inputs)} input tables.\n')
    print(f'Built readable report and {len(catalog)} figure pages in {reports}')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=ROOT);args=parser.parse_args();main(args.root)
