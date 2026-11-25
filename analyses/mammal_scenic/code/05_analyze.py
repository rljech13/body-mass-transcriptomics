from pathlib import Path
import os,json,pickle,hashlib
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd
np.object=object
from scipy.stats import hypergeom,pearsonr,t as tdist
from ctxcore.rnkdb import FeatherRankingDatabase
from pyscenic.aucell import aucell
R=Path(__file__).resolve().parents[1];P=R.parent;O=R/'results/analysis';O.mkdir(exist_ok=True);mods=['orange','violet','cyan','skyblue3'];seeds=[20261023,20261024,20261025]
def bh(x):
 x=np.asarray(x,float);assert np.isfinite(x).all();ix=np.argsort(x,kind='stable');q=np.minimum.accumulate((x[ix]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1];out=np.empty_like(q);out[ix]=np.minimum(q,1);return out
cfg=json.loads((R/'config/existing_reference.json').read_text());db=FeatherRankingDatabase(cfg['database'],name='existing');E=pd.read_pickle(R/'inputs/expression.pkl');meta=pd.read_csv(R/'inputs/samples.tsv',sep='\t').set_index('sample');M=pd.read_csv(R/'inputs/frozen_WGCNA_members.tsv',sep='\t');all_selected=set(M.loc[M.module.isin(mods),'feature']);groups={mod:set(M.loc[M.module.eq(mod),'feature']) for mod in mods};U=set(M.feature)&set(E.columns)&set(db.genes);(O/'background_genes.txt').write_text('\n'.join(sorted(U))+'\n')
traits=pd.read_csv(P/'growth_regulators_2026-10-21_v1/inputs/ma2016_full.samples.tsv',sep='\t').set_index('sample');species=traits.index.tolist();mass=np.log10(traits.loc[species,'mass_g'].to_numpy());assert set(species)==set(meta.species)
eig=pd.read_csv(P/'growth_regulators_2026-10-21_v1/results/wgcna/ma2016_full__raw/eigengenes.tsv',sep='\t').pivot(index='sample',columns='module',values='eigengene').loc[species]
# All statistical families specified before motif pruning, separate from seed diagnostics.
edges={};sets={};tests={}
for seed in seeds:
 a=pd.read_csv(R/'results'/str(seed)/'regulon_edges.tsv.gz',sep='\t');edges[seed]=a;sets[seed]={reg:set(t.target) for reg,t in a.groupby('regulon')};tfmap=a.groupby('regulon').TF.first().to_dict();rows=[]
 for reg,gg in sets[seed].items():
  tf=tfmap[reg];u=U-{tf};tg=gg&u
  for mod in mods:
   group=groups[mod]&u;hits=group&tg;N=len(u);K=len(tg);n=len(group);k=len(hits);p=float(hypergeom.sf(k-1,N,K,n))
   rows.append(dict(seed=seed,module=mod,regulon=reg,TF=tf,N_background=N,K_targets_background=K,n_module_covered=n,k_shared=k,n_module_total=len(groups[mod]),p=p,fold_enrichment=k/(n*K/N) if n*K else np.nan,genes=';'.join(sorted(hits))))
 t=pd.DataFrame(rows);t['q']=bh(t.p);t.to_csv(O/f'overlap_{seed}.tsv',sep='\t',index=False);tests[seed]=t
T=tests[seeds[0]].copy();eligible=T.K_targets_background.ge(5);T['q_ge5_diagnostic']=np.nan;T.loc[eligible,'q_ge5_diagnostic']=bh(T.loc[eligible,'p']);st=[]
for row in T.itertuples():
 base=sets[seeds[0]][row.regulon];ss=[sets[s].get(row.regulon,set()) for s in seeds];core=set.intersection(*ss);target_j=[len(base&ss[i])/len(base|ss[i]) if base|ss[i] else np.nan for i in [1,2]]
 st.append(dict(regulon=row.regulon,module=row.module,TF=row.TF,seeds_recovered=sum(bool(g) for g in ss),Jaccard_seed2=target_j[0],Jaccard_seed3=target_j[1],targets_all3=len(core),module_genes_all3=len(core&groups[row.module]),module_genes_all3_names=';'.join(sorted(core&groups[row.module])),overlap_q_seed2=float(tests[seeds[1]].set_index(['regulon','module']).q.get((row.regulon,row.module),1.)),overlap_q_seed3=float(tests[seeds[2]].set_index(['regulon','module']).q.get((row.regulon,row.module),1.))))
S=pd.DataFrame(st);S.to_csv(O/'seed_stability.tsv',sep='\t',index=False);T=T.merge(S,on=['regulon','module','TF']);T=T.sort_values(['q','p','module','TF']);T.to_csv(O/'module_regulon_results.tsv',sep='\t',index=False)
a=edges[seeds[0]].copy();a=a[a.target.isin(all_selected)].merge(M[['feature','module','kME','weight']].rename(columns={'feature':'target'}),on='target');a=a[a.module.isin(mods)];a['recovered_seed2']=[g in sets[seeds[1]].get(r,set()) for r,g in zip(a.regulon,a.target)];a['recovered_seed3']=[g in sets[seeds[2]].get(r,set()) for r,g in zip(a.regulon,a.target)];a['binding_measured']=False;a.to_csv(O/'module_regulon_genes.tsv',sep='\t',index=False)
# Standard SCENIC AUC activity is a rank-based signature score, not absolute TF activity.
auc=pd.read_csv(R/'results'/str(seeds[0])/'AUCell.tsv',sep='\t',index_col=0);auc_species=auc.groupby(meta.loc[auc.index,'species']).mean().loc[species];auc_species.to_csv(O/'AUCell_species.tsv',sep='\t');exprsp=E.groupby(meta.loc[E.index,'species']).mean().loc[species];activity=[]
for reg in auc_species.columns:
 y=auc_species[reg].to_numpy();r,p=pearsonr(mass,y) if np.ptp(y)>1e-12 else (np.nan,1.)
 tf=edges[seeds[0]].loc[edges[seeds[0]].regulon.eq(reg),'TF'].iloc[0]
 activity.append(dict(regulon=reg,TF=tf,n_species=len(species),mass_r=r,mass_p=p,own_RNA_mass_r=pearsonr(mass,exprsp[tf])[0] if tf in exprsp.columns and np.ptp(exprsp[tf])>1e-12 else np.nan))
A=pd.DataFrame(activity);A['mass_q']=bh(A.mass_p);A.to_csv(O/'regulon_mass_associations.tsv',sep='\t',index=False)
# Exclude all293 selected-module genes and TF itself from an auxiliary activity signature.
with (R/'results'/str(seeds[0])/'regulons.pkl').open('rb') as f:regs=pickle.load(f)
proxy=[];proxyrows=[]
for reg in regs:
 keep=set(reg.gene2weight)-all_selected-{reg.transcription_factor}
 if len(keep)>=5:proxy.append(reg.copy(gene2weight={g:reg.gene2weight[g] for g in keep}))
if proxy:
 pa=aucell(E,proxy,auc_threshold=.05,num_workers=8,seed=20261023);ps=pa.groupby(meta.loc[pa.index,'species']).mean().loc[species];ps.to_csv(O/'disjoint_AUCell_species.tsv',sep='\t')
 for row in T.itertuples():
  if row.regulon not in ps.columns:continue
  y=ps[row.regulon].to_numpy();z=eig[row.module].to_numpy();r,p=pearsonr(y,z) if np.ptp(y)>1e-12 else (np.nan,1.)
  design=np.column_stack([np.ones(len(mass)),mass]);yr=y-design@np.linalg.lstsq(design,y,rcond=None)[0];zr=z-design@np.linalg.lstsq(design,z,rcond=None)[0];pr=float(pearsonr(yr,zr)[0]) if np.std(yr)>1e-12 and np.std(zr)>1e-12 else np.nan;pp=float(2*tdist.sf(abs(pr)*np.sqrt((len(mass)-3)/max(1e-16,1-pr*pr)),len(mass)-3)) if np.isfinite(pr) else 1.
  proxyrows.append(dict(regulon=row.regulon,TF=row.TF,module=row.module,n_nonmodule_targets=len(sets[seeds[0]][row.regulon]-all_selected-{row.TF}),module_r=r,module_p=p,partial_mass_r=pr,partial_mass_p=pp))
 Q=pd.DataFrame(proxyrows);Q['module_q']=bh(Q.module_p);Q['partial_mass_q']=bh(Q.partial_mass_p);Q.to_csv(O/'disjoint_activity_support.tsv',sep='\t',index=False)
# Normalization diagnostic: project raw-normalized species means using frozen WGCNA PC weights (network not refit).
project=[]
for mod in mods:
 members=M[M.module.eq(mod)].set_index('feature');gg=sorted(set(members.index)&set(exprsp.columns));x=exprsp[gg];Z=(x-x.mean())/x.std(ddof=1);score=Z.to_numpy()@members.loc[gg,'weight'].to_numpy();project.append(dict(module=mod,n_genes=len(gg),r_with_frozen_eigengene=float(pearsonr(score,eig[mod])[0]),r_mass=float(pearsonr(score,mass)[0]),scope='Diagnostic projection on new normalization, frozen genes/weights; WGCNA not refit'))
pd.DataFrame(project).to_csv(O/'normalization_projection_diagnostic.tsv',sep='\t',index=False)
summary=dict(background=len(U),module_coverage={m:len(groups[m]&U) for m in mods},seeds={s:json.loads((R/'results'/str(s)/'SCENIC_summary.json').read_text()) for s in seeds},primary_overlap_tests=len(T),primary_overlap_hits=int(T.q.lt(.05).sum()),TFs_with_primary_module_hits=sorted(T.loc[T.q.lt(.05),'TF'].unique()),mass_association_tests=len(A),mass_association_hits=int(A.mass_q.lt(.05).sum()),top_module_results=T.head(15).replace({np.nan:None}).to_dict('records'),normalization_diagnostic=project,E2F_results=T[T.TF.str.startswith('E2F')].replace({np.nan:None}).to_dict('records'),causal_regulator_established=False,method='Full pySCENIC on26 Ma2016 bulk donor-labelled groups; preexisting human v10 short-promoter ranking; exploratory within-study results',postprocessing_revision='v2: evaluate all primary regulons x4 as declared; >=5-background-target family retained only as sensitivity; constant AUCell profiles assigned p=1 before BH after independent R check exposed NaNs in v1',code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(O/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False));print(json.dumps({k:v for k,v in summary.items() if k not in ['top_module_results','E2F_results']},indent=2),flush=True)
