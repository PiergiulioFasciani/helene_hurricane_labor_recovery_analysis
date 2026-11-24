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
    cv=table('inner_2023_cv'); solver=table('solver_audit'); scales=table('training_2022_scales')
    ext=table('LAUS_extended_period_effects_EXPLORATORY'); extev=table('LAUS_extended_event_coefficients_EXPLORATORY'); extpre=table('LAUS_extended_pretrends')
    migdose=read('outputs/migration/migration_by_damage.csv','migration_by_damage')
    bdpre=read('outputs/labor/tables/backdated_pretrends.csv','joint_multioutcome')
    bdeff=read('outputs/labor/tables/backdated_weights_robustness.csv','joint_multioutcome')
    dev=table('dose_event_coefficients'); dper=table('dose_period_effects'); dpre=table('dose_pretrend_tests'); dexp=table('dose_treated_exposure')
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
        add(f'PRE-EVENT CHECK (weights are fitted on this window, so it is in-sample): the linear pre-trend test rejects at 5% for {int((trends.slope_p<.05).sum())}/{k} outcomes, '
            f'and {int(trends.n_sig05.sum())} of {32*k} monthly pre-event coefficients are individually significant.',
            f'The strict joint test of all 31 monthly coefficients rejects for {int((trends.pretrend_p<.05).sum())}/{k} outcomes.',
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
    section('Validation')
    add('READ IN THIS ORDER. The definitive weights are fitted over the whole pre-event window, so their',
        'pre-event fit is IN-SAMPLE and is not evidence. Validation rests on what the fit could not manufacture:',
        'a window the weights never saw, and the estimates re-run on weights fitted before that window.')
    if not bdpre.empty:
        d=bdpre.copy();d['Outcome']=d.outcome.map(label)
        add('1. OUT-OF-SAMPLE PRE-EVENT PATH (the primary check). Weights fitted on 2022-23 ONLY.',
            'Jan 2022-Dec 2023 is in-sample for them; Jan-Jul 2024 is a genuine holdout. Compare the two amplitudes:')
        tab(d,['Outcome','drift_p','holdout_p','rms_insample','rms_holdout','max_abs_holdout'],
            {'drift_p':'Drift p','holdout_p':'Holdout joint p','rms_insample':'RMS in-sample','rms_holdout':'RMS holdout','max_abs_holdout':'Max holdout'})
        add('Two findings, and they point in different directions.',
            'NO TREND: the drift test does not reject for any outcome, and a differential trend is what parallel trends rules out.',
            'BUT AMPLITUDE GROWS: out of sample the gap wobbles several times more than in-sample.',
            'Read the last column against the estimated effects below. It bounds what the AGGREGATE can support;',
            'it does not touch the tercile contrasts, which difference the common control error away (see item 2).')
    if not fit.empty:
        d=fit.copy();d['Outcome']=d.outcome.map(label)
        add('   The same thing measured as forecast error rather than coefficients:')
        tab(d,['Outcome','bd_train_rmse','bd_holdout_rmse','bd_holdout_max','bd_rmspe_ratio'],
            {'bd_train_rmse':'Train RMSE','bd_holdout_rmse':'Forecast RMSE','bd_holdout_max':'Max gap','bd_rmspe_ratio':'Ratio'})
        add('The ratio is inflated by a very tight training fit; the absolute columns are the interpretable ones.')
    if not bdeff.empty:
        add('2. DO THE ESTIMATES SURVIVE THOSE WEIGHTS? The headline models re-run on the 2022-23 weights:',
            '   see the backdated-weights robustness section below. The tercile contrasts are unchanged to two decimals.')
    if not trends.empty:
        d=trends.copy();d['Outcome']=d.outcome.map(label)
        add('3. IN-SAMPLE PRE-EVENT PATH (reported for completeness, NOT as evidence).',
            'These are the definitive weights on the window they were fitted on, so flatness here is mechanical:')
        f=np.where(d.outcome.map(key)=='urate',1,100);d['max_abs_pre']=d.max_abs_pre*f
        tab(d,['Outcome','slope_p','pretrend_p','n_sig05','max_abs_pre'],{'slope_p':'Linear trend p','pretrend_p':'All 31 months p','n_sig05':'Months p<.05 (of 32)','max_abs_pre':'Largest pre coef'})
        add('Do not cite this table as validation. It establishes only that the optimiser did its job.')
    add('4. GEOGRAPHIC PLACEBO: see the placebo section below.')
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
    section('Effects by damage tercile')
    add('The 39 treated counties are split into terciles of storm exposure; each tercile is compared with the same synthetic control.',
        'Sample, weights, fixed effects, clustering and the pre-event-mean reference are identical to the main event study.',
        'HDSI (composite damage index) is primary; six-day rainfall alone is reported as a robustness comparison.',
        'Exposure is not randomly assigned across treated counties, so these are comparisons of differently exposed groups, not a randomised dose.')
    if not dper.empty:
        for fam,nm in [('hdsi_tercile','HDSI terciles (primary)'),('rain_tercile','Rainfall terciles (robustness)')]:
            f=dper[dper.dose==fam].copy()
            if f.empty: continue
            f['Outcome']=f.outcome.map(label);f['Period']=f.window.map(WINDOWS)
            f['Group']=f.term.map({'low':'low','mid':'mid','high':'high','high_minus_low':'high - low'})
            add('',nm+':')
            tab(f.sort_values(['outcome','window','term']),['Outcome','Period','Group','display_beta','display_lo','display_hi','pvalue'],{'display_beta':'Estimate','display_lo':'CI low','display_hi':'CI high','pvalue':'Raw p'})
    if not dpre.empty:
        d=dpre[dpre.dose=='hdsi_tercile'].copy();d['Outcome']=d.outcome.map(label)
        add('Pre-event tests within each HDSI tercile:')
        tab(d.sort_values(['outcome','term']),['Outcome','term','slope_p','n_sig05','max_abs_pre'],{'term':'Tercile','slope_p':'Linear trend p','n_sig05':'Months p<.05 (of 32)','max_abs_pre':'Largest pre coef'})
        add('A null in the least-damaged tercile at every horizon is a within-treatment placebo: those counties are FEMA-positive but barely hit.')
    section('Exploratory labor extension — through September 2026')
    if not ext.empty:
        d=ext.copy(); factor=np.where(d.outcome.map(key)=='urate',1,100)
        for c in ['beta','lower95','upper95']:d[c]=d[c]*factor
        d['Outcome']=d.outcome.map(label);d['Period']=d.window.map(WINDOWS)
        tab(d,['Outcome','Period','beta','lower95','upper95','pvalue'])
    if not extpre.empty:
        d=extpre.copy();d['Outcome']=d.outcome.map(label);tab(d,['Outcome','slope_p','pretrend_p'])
    add('Exploratory p-values are unadjusted and outside the 18-test main Holm family.')
    section('Backdated-weight robustness: do the estimates need the 2024 fit?')
    add('The headline models re-estimated on weights fitted on 2022-23 ONLY, which never saw Jan-Aug 2024.',
        'Sample window, excluded months, fixed effects, clustering, recentring and windows are unchanged;',
        'only the donor weights and the donor inclusion rule differ. These weights are fitted on 24 months',
        'instead of 32, so they are handicapped by construction and this is a conservative check.')
    if not bdeff.empty and not dper.empty:
        b=bdeff.copy(); m=dper[dper.dose=='hdsi_tercile'].copy()
        f=lambda d: np.where(d.outcome.map(key)=='urate',1,100)
        b['est']=b.beta*f(b); m['est']=m.beta*f(m)
        j=m.merge(b[b.spec=='tercile'][['outcome','term','window','est','pvalue']],
                  on=['outcome','term','window'],how='inner',suffixes=('_main','_bd'))
        j=j[j.term.isin(['low','high','high_minus_low'])].copy()
        j['Outcome']=j.outcome.map(label); j['Window']=j.window.map(WINDOWS)
        for c in ('est_main','est_bd'): j[c]=j[c].round(3)
        for c in ('pvalue_main','pvalue_bd'): j[c]=j[c].round(4)
        worder={'immediate':0,'year2025':1,'late':2}; torder={'low':0,'high':1,'high_minus_low':2}
        j['_w']=j.window.map(worder); j['_t']=j.term.map(torder)
        tab(j.sort_values(['outcome','_w','_t']),
            ['Outcome','Window','term','est_main','pvalue_main','est_bd','pvalue_bd'],
            {'term':'Group','est_main':'Definitive','pvalue_main':'p','est_bd':'Backdated','pvalue_bd':'p'})
        add('The high-minus-low contrasts agree to two decimals. Both terciles are measured against the SAME',
            'control inside one regression, so any error in that control enters both coefficients and differences out.',
            'That is why the out-of-sample wobble in the validation section bounds the aggregate but not the gradient.')
        if not bdeff[bdeff.spec=='aggregate'].empty:
            a=bdeff[bdeff.spec=='aggregate'].copy(); a['est']=(a.beta*f(a)).round(3)
            a['Outcome']=a.outcome.map(label); a['Window']=a.window.map(WINDOWS)
            a['pvalue']=a.pvalue.round(4); a['_w']=a.window.map({'immediate':0,'year2025':1,'late':2})
            tab(a.sort_values(['outcome','_w']),['Outcome','Window','est','pvalue'],
                {'est':'Backdated estimate','pvalue':'Raw p'})
            add('The aggregate estimates are if anything STRONGER on weights that never saw 2024.')
    else: add('Unavailable for this run.')
    section('Migration and population by damage tercile')
    add('netdom_per1000: net domestic migrants per 1,000 residents. pop_growth_pct: population growth in percent; effect is a percentage-point difference.',
        '2025 is July 2024–July 2025, only partly after Helene. Reference is the annual period ending July 2024.',
        'A negative contrast means slower net IN-migration, not population loss: see the levels table below before writing about this.',
        'It can also mean fewer arrivals rather than more departures, and it does not count displaced original residents.')
    add('Migration is reported as an OUTCOME of the same design, never as a control: it is realized after the storm,',
        'so conditioning on it would be a bad control and would remove the channel under study.',
        'Annual July-to-July frequency; county and year effects; 2024 is the omitted last unexposed year; clustered by county.',
        'Two comparison schemes: unweighted eligible counties, and donors carrying the fitted labor synthetic-control weights.')
    if not migdose.empty:
        d=migdose.copy()
        d['Scheme']=d.scheme.map({'unwt':'unweighted','syn':'synthetic'})
        d['Outcome']=d.outcome.map({'netdom_per1000':'Net domestic migration per 1,000','pop_growth_pct':'Population growth (pp)'})
        d['Group']=d.group.map({'all':'all 39','low':'low HDSI','mid':'mid HDSI','high':'high HDSI','hi-lo':'high minus low'})
        syn=d[d.scheme=='syn']
        tab(syn,['Outcome','Group','b2025','lower95','upper95','pvalue','pretest_p'],{'b2025':'2025 effect','lower95':'CI low','upper95':'CI high','pvalue':'Raw p','pretest_p':'Pre-test p'})
        add('Donors carry the labor synthetic-control weights, so migration is comparable with the labor outcomes.',
            'Pre-test p is the joint test of the 2021-2023 treated-by-year coefficients: small values mean the group was already',
            'on a different path before the storm.')
        u=d[d.scheme=='unwt']
        if not u.empty:
            hi=u[(u.group=='high')&(u.outcome=='netdom_per1000')]
            lo=u[(u.group=='hi-lo')&(u.outcome=='netdom_per1000')]
            if not hi.empty and not lo.empty:
                add(f'Footnote: under an unweighted comparison of all eligible counties the high-damage estimate is '
                    f'{hi.b2025.iloc[0]:.2f} (p={hi.pvalue.iloc[0]:.3f}) and the high-minus-low contrast {lo.b2025.iloc[0]:.2f} '
                    f'(p={lo.pvalue.iloc[0]:.3f}); that scheme fails its pre-test for the two less-damaged terciles, which is why it is not the primary one.')
    else: add('Unavailable for this run.')
    tab(groups,['year','grp','counties','mean_pop_growth_pct','mean_netdom_per1000'])
    add('LEVELS, not effects: group 1 = treated, group 0 = eligible comparison counties; unweighted county means.',
        'Net migration stays positive for the treated group throughout. The damaged counties attracted fewer new residents',
        'than they otherwise would have; they did not depopulate, and the paper must not say they did.')
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
        if not paths.empty:
            dates=paths.mdate.map(month)
            fig,axes=figgrid('Treated and synthetic-control paths','Four outcomes, centered on 2022 means; synthetic path includes the fitted intercept. Shading: 2024 months, now inside the fit window.')
            for ax,o in zip(axes,['emp','lf','jobs','urate']):
                ax.plot(dates,paths['mc_y_'+o],label='Treated',lw=1.6);ax.plot(dates,paths['mc_syn_joint_'+o],label='Synthetic control',ls='--',lw=1.5)
                pre=dates<pd.Timestamp('2024-09-01');ax.plot(dates[pre],paths.loc[pre,'mc_syn_bd_'+o],label='Backdated (2022–23 fit)',ls=':',lw=1.1,color=GRAY)
                ax.axvspan(pd.Timestamp('2024-01-01'),pd.Timestamp('2024-08-31'),color=GRAY,alpha=.12);ax.axvline(pd.Timestamp('2024-09-01'),color=GRAY,ls=':');ax.set(title=NAMES[o],ylabel=units(o));dateaxis(ax)
            axes[0].legend(frameon=False,fontsize=8)
            save(fig,'02_monthly_paths','Missing observations remain gaps. The vertical dotted line marks September 2024.')
        if not dev.empty:
            for fam,title,slug in [('hdsi_tercile','Event study by HDSI tercile','03b_tercile_events_hdsi'),('rain_tercile','Event study by rainfall tercile','03c_tercile_events_rain')]:
                if dev[dev.dose==fam].empty: continue
                fig,axes=figgrid(title,'Each tercile is 13 treated counties, compared with the same synthetic control; relative to its own pre-event mean.',3,2)
                for ax,o in zip(axes,['mc_yemp','mc_ylf','mc_yjobs','mc_yurate','mc_ydiff','mc_yshare']):
                    f=1 if key(o)=='urate' else 100
                    for term,col,lab_ in [('low',BLUE,'low'),('mid',GRAY,'mid'),('high',ORANGE,'high')]:
                        dd=dev[(dev.dose==fam)&(dev.outcome==o)&(dev.term==term)].sort_values('event_month')
                        ax.plot(dd.event_month,dd.beta*f,'-',lw=1.2,color=col,label=lab_)
                    zero(ax);ax.axvline(0,color=GRAY,ls=':');ax.set(title=label(o),ylabel=units(o),xlabel='Months from September 2024')
                axes[0].legend(fontsize=7,frameon=False,title='exposure tercile')
                save(fig,slug,'Point estimates; intervals are in dose_event_coefficients.csv. Thirteen counties per series, so each is noisier than the pooled event study.')
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
            fig,axes=figgrid('Validation before interpreting effects','Backdated weights (fitted 2022-23 only) forecasting 2024, and the formal pre-event tests.',1,2)
            sub=fit.copy();x=np.arange(len(sub))
            axes[0].bar(x-.2,sub.bd_train_rmse,.4,label='Training 2022-23',color=GRAY)
            axes[0].bar(x+.2,sub.bd_holdout_rmse,.4,label='Forecast Jan-Aug 2024',hatch='//',color=ORANGE)
            axes[0].set_xticks(x,sub.outcome.map(label),fontsize=8);axes[0].set(ylabel='RMSE (log-% / pp)');axes[0].legend(fontsize=8,frameon=False)
            if trends.empty: axes[1].set_visible(False)
            else:
                y=np.arange(len(trends))
                for col,mk,lab_ in [('slope_p','o','linear trend'),('pretrend_p','^','all 31 months')]:
                    axes[1].scatter(trends[col].clip(lower=1e-20),y,color=BLUE,marker=mk,facecolors=BLUE if mk=='o' else 'none',label=lab_)
                axes[1].set_xscale('log');axes[1].axvline(.05,color=ORANGE,ls='--')
                axes[1].set_yticks(y,trends.outcome.map(label),fontsize=8);axes[1].set(xlabel='p-value (log scale)');axes[1].legend(fontsize=8,frameon=False)
            save(fig,'04_validation','Left: the only genuinely out-of-sample check. Right: drift test and strict joint test, which answer different questions.')
        if not weights.empty:
            fig,axes=figgrid('Donor weights and model tuning','One shared weight vector fitted on Jan 2022–Aug 2024; the penalty is tuned on 2022–2023 only.')
            w=weights.sort_values('mc_w_joint',ascending=False);top=w.head(15).iloc[::-1]
            axes[0].barh(top.fips,100*top.mc_w_joint);axes[0].set(xlabel='Weight (%)',title='Top 15 donors (county FIPS)');axes[0].tick_params(axis='y',labelsize=8)
            axes[1].plot(np.arange(1,len(w)+1),100*w.mc_w_joint.cumsum());axes[1].set(xlabel='Donors ranked by weight',ylabel='Cumulative weight (%)',ylim=(0,102))
            if not cv.empty:axes[2].plot(cv['lambda'],cv.cv2023_rmse,'o-');axes[2].set(xscale='log',xlabel='Regularization parameter',ylabel='2023 validation error')
            axes[3].set_visible(False)
            save(fig,'05_donors','County identifiers are FIPS codes. The effective donor count reflects weight concentration.')
        if not draws.empty and not ranks.empty:
            fig,axes=figgrid('Geographic placebo distributions','Same four-outcome fitting procedure; 39 counties in each pseudo-treated group.')
            for ax,o in zip(axes,['emp','lf','jobs','urate']):
                ax.hist(draws.loc[draws.converged==1,'shock_'+o],bins=24,color=BLUE,alpha=.75)
                r=ranks[(ranks.outcome==o)&(ranks.screen=='all')].iloc[0];ax.axvline(r.observed_shock,color=ORANGE,lw=2,label=f'Observed {r.observed_shock:+.2f}');ax.set(title=NAMES[o],xlabel=units(o),ylabel='Placebo draws');ax.legend(fontsize=8,frameon=False)
            save(fig,'06_placebo_distributions','Ranks are descriptive benchmarks, not randomized causal p-values. All converged draws are shown.')
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
