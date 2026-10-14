from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.stats import wilcoxon
B=Path(__file__).resolve().parent.parent/'analysis_inputs'
C=['PPARGC1B','ESRRA','HNF4A','SDHAF1','MGME1'];defs=json.loads((B/'program_definitions.json').read_text())
def bh(v):
 v=np.asarray(v);i=np.argsort(v);r=np.empty(len(v));r[i]=np.minimum(1,np.minimum.accumulate((v[i]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]);return r
x=pd.read_csv(B/'external_counts.tsv',sep='\t',index_col=0);m=pd.read_csv(B/'external_samples.tsv',sep='\t',index_col=0);assert list(x.columns)==list(m.index);assert np.allclose(x.sum(),m.total_umi)
x=x.groupby(level=0).sum();l=np.log2(x.div(m.total_umi,axis=1)*1e6+1);rank=l.rank(axis=0,pct=True);cov=[];summ=[];pairs=[];scores=[]
for norm,a in [('CPM',l),('Rank',rank)]:
 s=a.loc[C].T.copy()
 for p,d in defs.items():
  gs=[g for g in d['genes_excluding_candidates'] if g in a.index];assert len(gs)/len(d['genes_excluding_candidates'])>=.8
  cov.append(dict(program=p,defined=len(d['genes_excluding_candidates']),matched=len(gs),missing=';'.join(set(d['genes_excluding_candidates'])-set(gs))))
  s[p]=a.loc[gs].mean(axis=0)
 s=s.join(m);s['normalization']=norm;scores.append(s.reset_index())
 for case in ['Adaptive PT','Degenerative PT']:
  for th in [10,20,30]:
   for stratum in ['All','CKD','AKI','Ref']:
    a1=s[(s.state==case)&(s.n_cells>=th)].set_index('donor');a0=s[(s.state=='Mature PT')&(s.n_cells>=th)].set_index('donor')
    if stratum!='All':a1=a1[a1.condition==stratum];a0=a0[a0.condition==stratum]
    ds=a1.index.intersection(a0.index)
    for f in C+list(defs):
     d=a1.loc[ds,f]-a0.loc[ds,f];p=wilcoxon(d,method='auto').pvalue if len(d)>=3 and (d!=0).any() else (1. if len(d)>=3 else np.nan)
     summ.append(dict(normalization=norm,comparison=case,reference='Mature PT',threshold=th,stratum=stratum,feature=f,family='candidate' if f in C else 'program',n_pairs=len(d),median_difference=d.median(),n_decreased=(d<0).sum(),n_increased=(d>0).sum(),p=p))
     if th==20 and stratum=='All':pairs += [dict(normalization=norm,comparison=case,feature=f,donor=dn,condition=a1.loc[dn,'condition'],difference=v) for dn,v in d.items()]
r=pd.DataFrame(summ);r['q']=np.nan
for _,g in r.groupby(['normalization','comparison','threshold','stratum','family']):
 if g.p.notna().all():r.loc[g.index,'q']=bh(g.p)
r.to_csv(B/'external_contrasts.csv',index=False);pd.DataFrame(pairs).to_csv(B/'external_differences.csv',index=False);pd.concat(scores).to_csv(B/'external_scores.csv',index=False);pd.DataFrame(cov).drop_duplicates().to_csv(B/'external_program_coverage.csv',index=False)
main=r[(r.normalization=='CPM')&(r.comparison=='Adaptive PT')&(r.threshold==20)&(r.stratum.isin(['All','CKD']))];print(main.to_string(index=False))
