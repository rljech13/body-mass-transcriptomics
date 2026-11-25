from pathlib import Path
import os,json,hashlib
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy.stats import pearsonr
R=Path(__file__).resolve().parents[1];O=R/'results/analysis';seeds=[20261023,20261024,20261025];meta=pd.read_csv(R/'inputs/samples.tsv',sep='\t').set_index('sample');traits=pd.read_csv(R/'inputs/species_traits.tsv',sep='\t').set_index('sample');sp=traits.index.tolist();x=np.log10(traits.mass_g.to_numpy())
def bh(a):
 a=np.asarray(a);ix=np.argsort(a);q=np.minimum.accumulate((a[ix]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];r=np.empty_like(q);r[ix]=np.minimum(q,1);return r
sets={};rows=[]
for seed in seeds:
 e=pd.read_csv(R/'results'/str(seed)/'regulon_edges.tsv.gz',sep='\t');sets[seed]={r:set(t.target) for r,t in e.groupby('regulon')};auc=pd.read_csv(R/'results'/str(seed)/'AUCell.tsv',sep='\t',index_col=0);a=auc.groupby(meta.loc[auc.index,'species']).mean().loc[sp];temp=[]
 for reg in a:
  y=a[reg].to_numpy();r,p=pearsonr(x,y) if np.ptp(y)>1e-12 else (np.nan,1.);temp.append(dict(seed=seed,regulon=reg,mass_r=r,mass_p=p))
 t=pd.DataFrame(temp);t['mass_q']=bh(t.mass_p);rows.append(t)
A=pd.concat(rows);A.to_csv(O/'mass_association_all_seeds.tsv',sep='\t',index=False)
base=sets[seeds[0]];j=[]
for reg,ss in base.items():
 gg=[sets[s].get(reg,set()) for s in seeds];r=dict(regulon=reg,TF=reg.split('(')[0],recovered=sum(bool(z) for z in gg),primary_targets=len(ss),targets_all3=len(set.intersection(*gg)))
 for i,s in enumerate(seeds[1:],2):r[f'Jaccard_seed{i}']=len(ss&gg[i-1])/len(ss|gg[i-1]) if ss|gg[i-1] else 0.
 j.append(r)
J=pd.DataFrame(j);J.to_csv(O/'regulon_stability.tsv',sep='\t',index=False)
T=pd.read_csv(O/'module_regulon_results.tsv',sep='\t');T['q_positive_runs']=(T.q.lt(.05).astype(int)+T.overlap_q_seed2.lt(.05).astype(int)+T.overlap_q_seed3.lt(.05).astype(int));T.to_csv(O/'module_regulon_results.tsv',sep='\t',index=False)
s={'primary_regulons':len(base),'recovered_in_all3':int(J.recovered.eq(3).sum()),'median_target_Jaccard_seed2':float(J.Jaccard_seed2.median()),'median_target_Jaccard_seed3':float(J.Jaccard_seed3.median()),'primary_significant_module_pairs':T[T.q.lt(.05)].replace({np.nan:None}).to_dict('records'),'module_pairs_q_positive_all3':int(T.q_positive_runs.eq(3).sum()),'mass_ERG_FOXP4_all_seeds':A[A.regulon.isin(['ERG(+)','FOXP4(+)'])].replace({np.nan:None}).to_dict('records'),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(O/'stability_summary.json').write_text(json.dumps(s,indent=2,allow_nan=False));print(json.dumps(s,indent=2))
