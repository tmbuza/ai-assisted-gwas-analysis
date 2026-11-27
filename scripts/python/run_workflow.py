"""Inspect quantitative-trait GWAS summary statistics; no genotype model fitting."""
import argparse,json,platform
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from common import load_inputs,sha256

sns.set_theme(style='whitegrid',context='notebook',font_scale=.95)
BLUE='#1976a3';TEAL='#009b8e';ORANGE='#d56a33';GREY='#c8ced4'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'figure.facecolor':'white','axes.titleweight':'bold','grid.alpha':.2,'grid.linewidth':.6})

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--config',type=Path,default=Path('config/default.json'))
    parser.add_argument('--check-only',action='store_true')
    args=parser.parse_args()
    try:
        t,meta=load_inputs(args.data_dir)
        cfg=json.loads(args.config.read_text())
        if not 0<=cfg['minimum_maf']<.5 or not 0<=cfg['minimum_info']<=1 or not 0<cfg['association_threshold']<1:
            raise ValueError('Invalid MAF, INFO or association thresholds.')
        if not isinstance(cfg['top_n'],int) or cfg['top_n']<1:raise ValueError('top_n must be a positive integer.')
        t['maf']=np.minimum(t.effect_allele_frequency,1-t.effect_allele_frequency)
        t['fails_maf']=t.maf<cfg['minimum_maf'];t['fails_info']=t.imputation_info<cfg['minimum_info']
        t['retained']=~(t.fails_maf|t.fails_info)
        t['palindromic']=t.apply(lambda r:set([r.effect_allele,r.other_allele]) in [set('AT'),set('CG')],axis=1)
        retained=t[t.retained].copy().sort_values(['chromosome','position'])
        if retained.empty:raise ValueError('No variants remain after filtering.')
    except (ValueError,OSError,KeyError,TypeError) as error:parser.exit(1,f'Input check failed: {error}\n')
    print(f'Input checks passed: {len(t):,} variant rows; {len(retained):,} retained; {t.chromosome.nunique()} autosomes represented.')
    if args.check_only:return
    figures=args.output_dir/'figures';tables=args.output_dir/'tables'
    figures.mkdir(parents=True,exist_ok=True);tables.mkdir(parents=True,exist_ok=True)
    # Remove only legacy exports of this workflow's named figures.
    for stem in ['manhattan','qq','effect-forest','frequency-info','effect-frequency','sample-size',
                 'filter-audit','precision-sample-size','p-consistency','regional-position']:
        for ext in ['pdf','svg']:
            legacy=figures/f'{stem}.{ext}'
            if legacy.exists():legacy.unlink()
    label='SYNTHETIC • fictional variants' if meta['synthetic'] else meta['study_id']
    def save(fig,name):
        fig.suptitle(label,fontsize=10,color='#52606d',y=.995)
        fig.tight_layout(rect=(0,0,1,.94))
        fig.savefig(figures/f'{name}.png',dpi=200)
        plt.close(fig)
    t.to_csv(tables/'variant-qc.tsv',sep='\t',index=False)
    retained['z_from_beta_se']=retained.beta/retained.standard_error
    retained['wald_p_from_beta_se']=2*stats.norm.sf(abs(retained.z_from_beta_se))
    retained['minus_log10_p']=-np.log10(retained.p_value)
    retained['ci_lower']=retained.beta-stats.norm.ppf(.975)*retained.standard_error
    retained['ci_upper']=retained.beta+stats.norm.ppf(.975)*retained.standard_error
    retained['passes_threshold']=retained.p_value<=cfg['association_threshold']
    retained['wald_log10_p_difference']=retained.minus_log10_p - (-stats.norm.logsf(abs(retained.z_from_beta_se))/np.log(10)-np.log10(2))
    retained.to_csv(tables/'variant-results.tsv',sep='\t',index=False)
    top=retained.nsmallest(cfg['top_n'],'p_value')
    top.to_csv(tables/'top-variants.tsv',sep='\t',index=False)
    threshold=-np.log10(cfg['association_threshold'])
    # Manhattan coordinates preserve chromosome order and physical spacing within chromosomes.
    offset=0; ticks=[];names=[];parts=[]
    for chromosome,frame in retained.groupby('chromosome',sort=True):
        frame=frame.copy();frame['plot_position']=frame.position+offset
        ticks.append((frame.plot_position.min()+frame.plot_position.max())/2);names.append(str(int(chromosome)))
        offset=frame.plot_position.max()+2000000;parts.append(frame)
    plotted=pd.concat(parts);plotted[['variant_id','chromosome','position','plot_position','minus_log10_p']].to_csv(tables/'manhattan-display.tsv',sep='\t',index=False)
    fig,ax=plt.subplots(figsize=(14,5))
    for c,frame in plotted.groupby('chromosome'):
        ax.scatter(frame.plot_position,frame.minus_log10_p,s=6,c='#244a73' if int(c)%2 else '#d68a35',alpha=.9,linewidths=0,rasterized=True)
    ax.axhline(threshold,c='#555555',ls='--',lw=1.1,label=f"Reference p={cfg['association_threshold']:g}")
    ax.set(xticks=ticks,xticklabels=names,xlabel='Autosome',ylabel='−log10(source p-value)',title='Association landscape • alternating autosome colours');ax.legend(frameon=False)
    ax.tick_params(axis='x',labelsize=8,labelrotation=45)
    ax.margins(y=.14)
    ax.grid(axis='x',visible=False)
    for _, frame in plotted[plotted.passes_threshold].groupby('chromosome'):
        row=frame.loc[frame.p_value.idxmin()]
        ax.annotate(row.variant_id,(row.plot_position,row.minus_log10_p),xytext=(4,7),textcoords='offset points',fontsize=7)
    save(fig,'manhattan')
    # QQ only for a declared genome-wide file, never a hits-only selection.
    diagnostic={}
    if meta['summary_scope']=='genome_wide':
        pvals=np.sort(retained.p_value.to_numpy());m=len(pvals);ranks=np.arange(1,m+1)
        expected=-np.log10((ranks-.5)/m);observed=-np.log10(pvals)
        lower=-np.log10(stats.beta.ppf(.975,ranks,m+1-ranks));upper=-np.log10(stats.beta.ppf(.025,ranks,m+1-ranks))
        qq=pd.DataFrame(dict(expected_minus_log10_p=expected,observed_minus_log10_p=observed,pointwise_lower=lower,pointwise_upper=upper))
        qq.to_csv(tables/'qq-display.tsv',sep='\t',index=False)
        inflation=float(np.median(stats.chi2.isf(pvals,1))/stats.chi2.ppf(.5,1))
        diagnostic={'lambda_gc':inflation,'lambda_definition':'median(chi2.isf(source p,1))/chi2.ppf(0.5,1); 1-df diagnostic only',
                    'qq_envelope':'95% pointwise under independent Uniform(0,1) p-values; not simultaneous and not LD-adjusted'}
        fig,axes=plt.subplots(1,2,figsize=(11,4.8));order=np.argsort(expected)
        for ax in axes:
            ax.fill_between(expected[order],lower[order],upper[order],color=GREY,alpha=.6,label='95% pointwise independent-null envelope')
            ax.scatter(expected,observed,s=6,color=BLUE,alpha=.6,rasterized=True)
            limit=expected.max();ax.plot([0,limit],[0,limit],c=ORANGE,ls='--',lw=1)
            ax.set(xlabel='Expected −log10(p)',ylabel='Observed −log10(p)')
        axes[0].set(xlim=(0,limit*1.08),ylim=(0,max(observed.max(),limit)*1.09),title='Full association tail')
        axes[0].text(.04,.94,f'λGC = {inflation:.3f}',transform=axes[0].transAxes,va='top')
        axes[1].set(xlim=(0,2),ylim=(0,2),title='Central distribution • same variants')
        axes[1].set_aspect('equal',adjustable='box')
        axes[0].legend(fontsize=6.5,loc='lower right',frameon=False)
        save(fig,'qq')
    else:
        # Avoid stale diagnostics when a folder previously held a full-file run.
        for old in [figures/f'qq.{ext}' for ext in ['png','pdf','svg']]+[tables/'qq-display.tsv']:
            if old.exists():old.unlink()
        diagnostic={'qq_status':'Not generated for selected variants; do not infer global calibration.'}
    fig,ax=plt.subplots(figsize=(12,6));y=np.arange(len(top))
    ax.errorbar(top.beta,y,xerr=stats.norm.ppf(.975)*top.standard_error,fmt='o',color=BLUE,ecolor=TEAL,capsize=3)
    ax.axvline(0,c=GREY,ls='--');ax.set(yticks=y,yticklabels=[f'{r.variant_id} | chr{int(r.chromosome)}:{int(r.position):,} | {r.effect_allele}/{r.other_allele}' for r in top.itertuples()],xlabel=f"Beta ({meta['effect_units']})",title='Top variants • position, effect/other alleles and nominal 95% intervals');ax.invert_yaxis()
    ax.tick_params(axis='y',labelsize=8)
    save(fig,'effect-forest')
    fig,axes=plt.subplots(1,2,figsize=(11,4.4))
    maf_bins=np.linspace(0,.5,31);info_bins=np.linspace(0,1,31)
    axes[0].hist(t.maf,bins=maf_bins,histtype='step',linewidth=1.6,color='#7b858f',label='All input rows')
    sns.histplot(x=retained.maf,bins=maf_bins,color=TEAL,alpha=.65,ax=axes[0],label='Retained')
    axes[0].axvline(cfg['minimum_maf'],color=ORANGE,ls='--',label=f"MAF ≥ {cfg['minimum_maf']:g}")
    axes[0].set(xlabel='Minor allele frequency',ylabel='Variant count',title='Frequency mixture • shared histogram bins');axes[0].legend(frameon=False,fontsize=7)
    sns.histplot(x=t.imputation_info,bins=info_bins,color=BLUE,ax=axes[1])
    axes[1].axvline(cfg['minimum_info'],color=ORANGE,ls='--',label=f"INFO ≥ {cfg['minimum_info']:g}")
    axes[1].set(xlabel='Source INFO score',ylabel='Variant count',title='Imputation quality before filtering');axes[1].legend(frameon=False,fontsize=7)
    pd.DataFrame({'bin_lower':maf_bins[:-1],'bin_upper':maf_bins[1:],
        'input_variants':np.histogram(t.maf,bins=maf_bins)[0],
        'retained_variants':np.histogram(retained.maf,bins=maf_bins)[0]}).to_csv(tables/'maf-histogram-display.tsv',sep='\t',index=False)
    save(fig,'frequency-info')
    fig,ax=plt.subplots(figsize=(8,4.5));hit=retained.passes_threshold
    ax.scatter(retained.loc[~hit,'maf'],retained.loc[~hit,'beta'],s=6,c=GREY,alpha=.5,rasterized=True)
    ax.scatter(retained.loc[hit,'maf'],retained.loc[hit,'beta'],s=13,c=ORANGE,label=f"p ≤ {cfg['association_threshold']:g}",alpha=.8)
    ax.axhline(0,c=BLUE,lw=.8);ax.set(xlabel='Minor allele frequency',ylabel='Beta (source units)',title='Effect estimates across allele frequencies');ax.legend(frameon=False,fontsize=8);save(fig,'effect-frequency')
    fig,ax=plt.subplots(figsize=(8,4.5));sns.histplot(x=retained.sample_size,bins=30,color=BLUE,ax=ax)
    ax.axvline(retained.sample_size.median(),color=ORANGE,ls='--',label=f'Median N={retained.sample_size.median():,.0f}');ax.legend(frameon=False,fontsize=8)
    ax.set(xlabel='Per-variant sample size',ylabel='Variant count',title='Sample size can vary across variants');save(fig,'sample-size')
    stages=pd.DataFrame({'stage':['Input','After MAF','After MAF + INFO'], 'variants':[len(t),int((~t.fails_maf).sum()),len(retained)]})
    stages['excluded_at_step']=np.r_[0,-np.diff(stages.variants)]
    stages['excluded_percent_of_previous']=np.r_[0,100*(-np.diff(stages.variants))/stages.variants.to_numpy()[:-1]]
    stages['retained_percent_of_input']=100*stages.variants/len(t)
    stages.to_csv(tables/'filter-summary.tsv',sep='\t',index=False)
    fig,ax=plt.subplots(figsize=(10,4.5));ax.barh(stages.stage,stages.variants,color=[GREY,BLUE,TEAL]);ax.invert_yaxis()
    for i,row in stages.iterrows():
        note='Input reference' if i==0 else f"removed {int(row.excluded_at_step):,} ({row.excluded_percent_of_previous:.1f}% of prior step)"
        ax.text(row.variants,i,f' {row.variants:,} • {note}',va='center',fontsize=8)
    ax.set(xlabel='Variant rows (not participants)',title='Filtering audit • sequential counts',xlim=(0,len(t)*1.85));ax.grid(axis='y',visible=False);save(fig,'filter-audit')
    # Precision is described alongside MAF; this is not a fitted association model.
    precision=retained[['variant_id','sample_size','standard_error','maf','imputation_info']].copy()
    precision.to_csv(tables/'precision-display.tsv',sep='\t',index=False)
    fig,ax=plt.subplots(figsize=(8,4.7))
    points=ax.scatter(retained.sample_size,retained.standard_error,c=retained.maf,cmap='viridis',s=7,alpha=.5,rasterized=True,vmin=0,vmax=.5)
    ax.set(yscale='log',xlabel='Per-variant N',ylabel='Standard error (log scale)',title='Precision depends on coverage and allele frequency')
    fig.colorbar(points,ax=ax,label='Minor allele frequency');save(fig,'precision-sample-size')
    # Compute log p stably even if the derived probability underflows to zero.
    retained['minus_log10_wald_p']=-stats.norm.logsf(abs(retained.z_from_beta_se))/np.log(10)-np.log10(2)
    retained[['variant_id','minus_log10_p','minus_log10_wald_p','wald_log10_p_difference']].to_csv(tables/'p-consistency-display.tsv',sep='\t',index=False)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5));limit=max(retained.minus_log10_p.max(),retained.minus_log10_wald_p.max())
    axes[0].scatter(retained.minus_log10_p,retained.minus_log10_wald_p,s=7,color=TEAL,alpha=.45,rasterized=True)
    axes[0].plot([0,limit],[0,limit],ls='--',color=ORANGE)
    axes[0].set(xlabel='−log10(source p)',ylabel='−log10(derived Wald p)',title='Agreement depends on source test')
    axes[1].scatter(retained.minus_log10_p,retained.wald_log10_p_difference,s=6,color=BLUE,alpha=.4,rasterized=True)
    axes[1].axhline(0,color=ORANGE,ls='--')
    axes[1].ticklabel_format(axis='y',style='sci',scilimits=(-3,3))
    axes[1].set(xlabel='−log10(source p)',ylabel='Source −log10(p) minus Wald −log10(p)',title='Difference panel • inspect magnitude')
    max_difference=retained.wald_log10_p_difference.abs().max()
    if max_difference<1e-8:
        axes[1].set_ylim(-1e-6,1e-6)
        axes[1].set_title('Agreement at numerical rounding scale')
    axes[1].text(.03,.94,f'Max |Δ| = {max_difference:.2e}',transform=axes[1].transAxes,va='top',fontsize=8)
    save(fig,'p-consistency')
    # No LD colours or gene labels are invented for a position-only zoom.
    strongest=retained.loc[retained.p_value.idxmin()];region=retained[(retained.chromosome==strongest.chromosome)&(abs(retained.position-strongest.position)<=250000)]
    region.to_csv(tables/'region-display.tsv',sep='\t',index=False)
    fig,ax=plt.subplots(figsize=(9,4.5));ax.scatter(region.position/1e6,region.minus_log10_p,s=16,color=TEAL,alpha=.8)
    ax.scatter([strongest.position/1e6],[strongest.minus_log10_p],marker='*',s=160,color=ORANGE,label=strongest.variant_id)
    ax.axhline(threshold,c=GREY,ls='--');ax.set(xlabel=f'Chromosome {int(strongest.chromosome)} position (Mb; {meta["genome_build"]})',ylabel='−log10(source p-value)',title='Position-only zoom • ±250 kb around smallest p-value');ax.legend(frameon=False);save(fig,'regional-position')
    summary=dict(input_variants=len(t),retained_variants=len(retained),excluded_variants=len(t)-len(retained),
        at_or_below_p_threshold_variants=int(retained.passes_threshold.sum()),palindromic_retained=int(retained.palindromic.sum()),
        max_absolute_wald_log10_p_difference=float(retained.wald_log10_p_difference.abs().max()),**diagnostic)
    (tables/'analysis-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    paths=[args.data_dir/'summary-statistics.tsv',args.data_dir/'study-metadata.json',args.config]
    record={'inputs_sha256':{str(f):sha256(f) for f in paths},'code_sha256':{f.name:sha256(f) for f in Path(__file__).parent.glob('*.py')},'config':cfg,'study_metadata':meta,'summary':summary,
            'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'seaborn':sns.__version__},
            'limitations':['No genotype association model, LD clumping, fine mapping, causal inference or ancestry correction is performed.','Top variants are correlated candidates, not independent loci.','95% intervals are unadjusted post-selection descriptive intervals.']}
    (tables/'run-record.json').write_text(json.dumps(record,indent=2)+'\n')
    print(f'Results saved to {args.output_dir}')

if __name__=='__main__':main()
