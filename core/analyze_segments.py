from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.stats import wilcoxon
B=Path(__file__).resolve().parent.parent/'analysis_inputs';E=B
C=['PPARGC1B','ESRRA','HNF4A','SDHAF1','MGME1'];defs=json.loads((B/'program_definitions.json').read_text())
def bh(v):
 v=np.asarray(v,float);out=np.full(len(v),np.nan);ix=np.flatnonzero(np.isfinite(v));o=ix[np.argsort(v[ix])];n=len(o)
 if n:out[o]=np.minimum(1,np.minimum.accumulate((v[o]*n/np.arange(1,n+1))[::-1])[::-1])
 return out
x=pd.read_csv(E/'segment_counts.tsv',sep='\t',index_col=0).groupby(level=0).sum();m=pd.read_csv(E/'segment_samples.tsv',sep='\t',index_col=0);assert list(x.columns)==list(m.index);assert np.allclose(x.sum(),m.total_umi);l=np.log2(x.div(m.total_umi,axis=1)*1e6+1);summ=[];pairs=[];allscore=[]
for norm,a in [('CPM',l),('Rank',l.rank(axis=0,pct=True))]:
 s=a.loc[C].T.copy()
 for p,d in defs.items():s[p]=a.loc[[g for g in d['genes_excluding_candidates'] if g in a.index]].mean(axis=0)
 s=s.join(m);s['normalization']=norm;allscore.append(s.reset_index())
 for ref in ['PT-S1','PT-S2','PT-S3']:
  for th in [10,20,30]:
   for stratum in ['All','CKD']:
    a1=s[(s.state=='Adaptive PT')&(s.n_cells>=th)].set_index('donor');a0=s[(s.state==ref)&(s.n_cells>=th)].set_index('donor')
    if stratum=='CKD':a1=a1[a1.condition=='CKD'];a0=a0[a0.condition=='CKD']
    ds=a1.index.intersection(a0.index)
    for f in C+list(defs):
     d=a1.loc[ds,f]-a0.loc[ds,f];p=wilcoxon(d,method='auto').pvalue if len(d)>=3 and (d!=0).any() else (1. if len(d)>=3 else np.nan)
     summ.append(dict(normalization=norm,reference=ref,threshold=th,stratum=stratum,feature=f,family='candidate' if f in C else 'program',n_pairs=len(d),median_difference=d.median(),n_decreased=(d<0).sum(),n_increased=(d>0).sum(),p=p))
     if th==20:pairs += [dict(normalization=norm,reference=ref,stratum=stratum,feature=f,donor=dn,condition=a1.loc[dn,'condition'],difference=v) for dn,v in d.items()]
r=pd.DataFrame(summ);r['family_q']=np.nan;r['ppargc1b_three_reference_q']=np.nan
for _,g in r.groupby(['normalization','reference','threshold','stratum','family']):r.loc[g.index,'family_q']=bh(g.p)
for _,g in r[r.feature=='PPARGC1B'].groupby(['normalization','threshold','stratum']):r.loc[g.index,'ppargc1b_three_reference_q']=bh(g.p)
r.to_csv(B/'segment_contrasts.csv',index=False);pd.DataFrame(pairs).to_csv(B/'segment_differences.csv',index=False);pd.concat(allscore).to_csv(B/'segment_scores.csv',index=False)
print(r[(r.normalization=='CPM')&(r.feature=='PPARGC1B')&(r.threshold==20)].to_string(index=False))
# Core cross-cohort effect summaries; each cohort retains its state definition.
p=pd.read_csv(B/'program_candidate_paired_differences.csv');p=p[(p.normalization=='CPM')&(((p.cohort=='Discovery')&(p.reference=='PT Injured'))|((p.cohort=='Validation')&(p.reference=='PT')))]
p=p[['cohort','feature','donor','condition','difference']];e=pd.read_csv(B/'external_differences.csv');e=e[(e.normalization=='CPM')&(e.comparison=='Adaptive PT')][['feature','donor','condition','difference']];e['cohort']='External';p=pd.concat([p,e],ignore_index=True)
rng=np.random.default_rng(61013);out=[]
for (cohort,f),g in p.groupby(['cohort','feature'],sort=True):
 v=g.sort_values('donor').difference.to_numpy();boot=np.median(v[rng.integers(0,len(v),size=(10000,len(v)))],axis=1);lo,hi=np.quantile(boot,[.025,.975]);out.append(dict(cohort=cohort,feature=f,n_pairs=len(v),median_difference=np.median(v),bootstrap95_low=lo,bootstrap95_high=hi,n_decreased=(v<0).sum(),n_increased=(v>0).sum(),bootstrap_resamples=10000,seed=61013))
p.to_csv(B/'cross_cohort_donor_differences.csv',index=False);pd.DataFrame(out).to_csv(B/'cross_cohort_effects.csv',index=False);print('Cross cohort effects computed with donor resampling')
