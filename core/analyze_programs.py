from pathlib import Path
import sys,json
import numpy as np,pandas as pd
from scipy.stats import wilcoxon,spearmanr,rankdata
sys.stdout.reconfigure(encoding='utf8',errors='replace')
B=Path(__file__).resolve().parent.parent;V=B/'analysis_inputs'
defs=json.loads((V/'program_definitions.json').read_text(encoding='utf8'))
genes=['PPARGC1B','ESRRA','HNF4A','SDHAF1','MGME1']
def bh(p):
 p=np.asarray(p,float);out=np.full(len(p),np.nan);ix=np.flatnonzero(np.isfinite(p));order=ix[np.argsort(p[ix])];n=len(order)
 if n:out[order]=np.minimum(1,np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1])
 return out
def load_filtered(path,states):return pd.concat([x[x.state.isin(states)] for x in pd.read_csv(path,chunksize=250000)],ignore_index=True)
summaries=[];pairrows=[];scoreall=[];coverage=[];correlations=[]
for cohort,path,states in [('Discovery',V/'discovery_all_gene_pt_pseudobulk.csv.gz',['PT Inflammatory','PT Injured','PT S1','PT S2','PT S3']),('Validation',V/'all_gene_pseudobulk.csv.gz',['PTVCAM1','PT'])]:
 x=load_filtered(path,states);x['sample']=x.donor+'|'+x.state
 meta=x.groupby('sample').agg(donor=('donor','first'),state=('state','first'),condition=('condition','first'),n_cells=('n_cells','first'),total_umi=('total_umi','first'))
 common=x.groupby('gene')['sample'].nunique();measured=set(common[common==len(meta)].index)
 mat=x.pivot(index='gene',columns='sample',values='log2cpm')
 ranks=mat.loc[sorted(measured)].rank(axis=0,method='average',pct=True)
 for normalization,matrix in {'CPM':mat,'Rank':ranks}.items():
  scores=matrix.loc[[g for g in genes if g in matrix.index]].T.copy()
  for program,d in defs.items():
   gs=[g for g in d['genes_excluding_candidates'] if g in measured]
   coverage.append(dict(cohort=cohort,normalization=normalization,program=program,defined_genes=len(d['genes_excluding_candidates']),matched_genes=len(gs),excluded_unmatched=';'.join(g for g in d['genes_excluding_candidates'] if g not in measured)))
   assert len(gs)/len(d['genes_excluding_candidates'])>=.8
   scores[program]=matrix.loc[gs].mean(axis=0)
  scores=scores.join(meta);scores['cohort']=cohort;scores['normalization']=normalization;scoreall.append(scores.reset_index())
  refs=['PT Injured','PT S1','PT S2','PT S3'] if cohort=='Discovery' else ['PT']
  for ref in refs:
   for threshold in [10,20,30]:
    for stratum in (['UUO'] if cohort=='Discovery' else ['All','DKD','Control']):
     ss=scores[scores.n_cells>=threshold]
     if stratum!='All':ss=ss[ss.condition==stratum]
     aa=ss[ss.state==states[0]].set_index('donor');bb=ss[ss.state==ref].set_index('donor');ds=aa.index.intersection(bb.index)
     for feature in genes+list(defs):
      diff=(aa.loc[ds,feature]-bb.loc[ds,feature]).dropna();n=len(diff)
      p=float(wilcoxon(diff,method='auto').pvalue) if n>=3 and (diff!=0).any() else np.nan
      loo=[float(diff.drop(d).median()) for d in diff.index] if n>1 else []
      summaries.append(dict(cohort=cohort,normalization=normalization,reference=ref,threshold=threshold,stratum=stratum,feature=feature,family='candidate' if feature in genes else 'program',n_pairs=n,median_difference=float(diff.median()) if n else np.nan,n_decreased=int((diff<0).sum()),n_increased=int((diff>0).sum()),wilcoxon_p=p,loo_min=min(loo) if loo else np.nan,loo_max=max(loo) if loo else np.nan))
      if threshold==20 and stratum in ['All','UUO']:pairrows += [dict(cohort=cohort,normalization=normalization,reference=ref,feature=feature,donor=d,condition=str(aa.loc[d,'condition']),difference=float(v)) for d,v in diff.items()]
  for state in states:
   ss=scores[(scores.state==state)&(scores.n_cells>=20)]
   if len(ss)<5:continue
   for gene in genes:
    for program in ['PT_TRANSPORT_SEED','HALLMARK_OXIDATIVE_PHOSPHORYLATION']:
     rho,p=spearmanr(ss[gene],ss[program]);correlations.append(dict(cohort=cohort,normalization=normalization,state=state,gene=gene,program=program,n=len(ss),rho=float(rho),p=float(p),type='within_state'))
     if ss.condition.nunique()==2:
      a=rankdata(ss[gene]);b=rankdata(ss[program]);z=np.column_stack([np.ones(len(ss)),(ss.condition==ss.condition.unique()[0]).astype(float)])
      ra=a-z@np.linalg.lstsq(z,a,rcond=None)[0];rb=b-z@np.linalg.lstsq(z,b,rcond=None)[0]
      correlations.append(dict(cohort=cohort,normalization=normalization,state=state,gene=gene,program=program,n=len(ss),rho=float(np.corrcoef(ra,rb)[0,1]),p=np.nan,type='condition_adjusted_rank_residual_descriptive'))
 print(cohort,'summarized',len(meta),'donor/states',flush=True)
res=pd.DataFrame(summaries);res['family_bh_q']=np.nan;res['all_10_features_bh_q']=np.nan
for _,g in res.groupby(['cohort','normalization','reference','threshold','stratum','family']):res.loc[g.index,'family_bh_q']=bh(g.wilcoxon_p)
for _,g in res.groupby(['cohort','normalization','reference','threshold','stratum']):res.loc[g.index,'all_10_features_bh_q']=bh(g.wilcoxon_p)
res.to_csv(V/'program_candidate_contrasts.csv',index=False);pd.DataFrame(pairrows).to_csv(V/'program_candidate_paired_differences.csv',index=False)
pd.concat(scoreall).to_csv(V/'donor_program_scores.csv',index=False);pd.DataFrame(coverage).to_csv(V/'program_coverage.csv',index=False)
c=pd.DataFrame(correlations);c['bh_q']=np.nan
for _,g in c[c.type=='within_state'].groupby(['cohort','normalization']):c.loc[g.index,'bh_q']=bh(g.p)
c.to_csv(V/'within_state_correlations.csv',index=False)
print(res[(res.cohort=='Validation')&(res.normalization=='CPM')&(res.threshold==20)&(res.stratum=='All')].to_string(index=False))
