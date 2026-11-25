from pathlib import Path
import json,hashlib,sys
import numpy as np,pandas as pd
from scipy import stats,linalg
O=Path(__file__).resolve().parents[1];R=O.parent;H=Path('/path/to/user-home/chicken/kleiber_birds_gemcat')
def bh(p):
 a=np.asarray(p,float);out=np.full(len(a),np.nan);ok=np.flatnonzero(np.isfinite(a));ii=ok[np.argsort(a[ok])];out[ii]=np.minimum(1,np.minimum.accumulate((a[ii]*len(ii)/np.arange(1,len(ii)+1))[::-1])[::-1]);return out
def fit(x,y,C,extra=None,loo=False):
 x=np.asarray(x,float);y=np.asarray(y,float)
 if not np.isfinite(x).all() or not np.isfinite(y).all() or np.std(y)<1e-12 or len(x)<5:return dict(beta=np.nan,p=np.nan,r=np.nan,n=len(x))
 X=np.column_stack([np.ones(len(x)),x] + ([] if extra is None else [np.asarray(extra)]));L=np.linalg.cholesky(C);W=linalg.solve_triangular(L,X,lower=True);v=linalg.solve_triangular(L,y,lower=True);b=np.linalg.lstsq(W,v,rcond=None)[0];df=len(x)-np.linalg.matrix_rank(W);err=v-W@b;se=np.sqrt((err@err/df)*np.linalg.pinv(W.T@W)[1,1]);t=b[1]/max(se,1e-30)
 z=dict(beta=float(b[1]),p=float(2*stats.t.sf(abs(t),df)),r=float(t/np.sqrt(t*t+df)),n=len(x),df=int(df))
 if loo:
  zz=[]
  for i in range(len(x)):
   ix=np.arange(len(x))!=i;zz.append(fit(x[ix],y[ix],C[np.ix_(ix,ix)],None if extra is None else np.asarray(extra)[ix]))
  z.update(loo_sign_fraction=float(np.mean([np.sign(t['beta'])==np.sign(b[1]) for t in zz])),loo_max_p=float(max(t['p'] for t in zz)),loo_min_abs_r=float(min(abs(t['r']) for t in zz)))
 return z
def load(panel):
 M=pd.read_csv(O/'inputs'/f'{panel}.log.tsv',sep='\t',index_col=0);t=pd.read_csv(O/'inputs'/f'{panel}.samples.tsv',sep='\t');t=t.rename(columns={t.columns[0]:'sample'})
 f=pd.read_csv(O/'inputs'/f'{panel}.features.tsv',sep='\t');f=f.rename(columns={f.columns[0]:'feature'}).set_index('feature')
 sp=list(dict.fromkeys(t.species));X=M.T.join(t.set_index('sample').species).groupby('species',sort=False).mean().loc[sp].T
 cov=np.load(O/'inputs'/f'{panel}.cov.npz');assert list(cov['species'])==sp
 trait=t.drop_duplicates('species').set_index('species').loc[sp];C=cov['C'];C=C/np.mean(np.diag(C))
 # strict single-symbol gene layer; family orthogroups retained separately in WGCNA.
 fs=f[f.symbols.notna() & ~f.symbols.str.contains(';',na=False)];fs=fs[~fs.symbols.duplicated(keep=False)]
 G=X.loc[X.index.intersection(fs.index)].rename(index=fs.symbols.to_dict());Z=G.sub(G.mean(axis=1),axis=0).div(G.std(axis=1),axis=0);Z=Z[np.isfinite(Z).all(axis=1)]
 return dict(panel=panel,M=M,meta=t,feature=f,X=X,G=G,Z=Z,C=C,sp=sp,mass=np.log10(trait.mass_g.to_numpy(float)),background=Z.loc[~Z.index.isin(set(sum(json.loads((O/'config/fixed_sets.json').read_text())['bird' if panel.startswith('bird') else 'mammal'].values(),[]))|set(json.loads((O/'config/analysis_plan.initial.json').read_text())['candidate_TFs']))].mean(axis=0).to_numpy())
def save(name,rr):
 d=pd.DataFrame(rr)
 if 'p' in d:d['q_all_tests']=bh(d.p)
 d.to_csv(O/'results'/('full_'+name),sep='\t',index=False);return d
plan=json.loads((O/'config/analysis_plan.initial.json').read_text());TF=plan['candidate_TFs'];sets=json.loads((O/'config/fixed_sets.json').read_text())
inv=pd.read_csv(O/'config/full_inventory.tsv',sep='\t');P={p:load(p) for p in inv.panel}
priorpath=H/'results/lab_regulator_panel_2026-09-29/reference/CollecTRI_regulons.csv';prior=pd.read_csv(priorpath).groupby(['source','target'],as_index=False).weight.mean();prior=prior[prior.source.isin(TF)&prior.weight.ne(0)].copy()
prior.to_csv(O/'reference/candidate_CollecTRI_edges.tsv',sep='\t',index=False)
def groups(panel):return sets['bird' if panel.startswith('bird') else 'mammal']
def group_features(p,gg):
 return p['feature'].index[p['feature'].symbols.isin(gg)].intersection(p['X'].index).tolist()
def syms(p,features):
 return set(p['feature'].reindex(features).symbols.dropna().loc[lambda s:~s.str.contains(';')])
module_mass=[];module_groups=[];module_tf=[];projection=[];tfmass=[];tfgroup=[];coverage=[];edges=[];scores=[];learnedmass=[];learnedmembership=[]
for panel,p in P.items():
 Z=p['Z'];story=set(g for v in groups(panel).values() for g in v);cx=p['C'];mass=p['mass'];print('ANALYSE',panel,len(Z),flush=True)
 # Prior-defined RNA footprints with complete exclusion of every fixed story set and TF self.
 for tf,a in prior.groupby('source'):
  aa=a[a.target.isin(Z.index)&~a.target.isin(story|{tf})];coverage.append(dict(panel=panel,TF=tf,n_prior=len(a),n_measured_no_story=len(aa),TF_measured=tf in Z.index))
  if len(aa)<10:continue
  w=aa.set_index('target').weight;raw=(Z.loc[w.index].T@w.to_numpy())/np.abs(w).sum();sc=(Z.loc[w.index].sub(p['background'],axis=1).T@w.to_numpy())/np.abs(w).sum()
  for typ,sc0 in [('raw_signed',raw),('background_centered_signed',sc)]:
   y=np.asarray(sc0);tfmass.append(dict(panel=panel,TF=tf,score=typ,n_targets=len(aa),**fit(mass,y,cx,loo=True),ols_p=fit(mass,y,np.eye(len(mass)))['p']))
   for i,sp in enumerate(p['sp']):scores.append(dict(panel=panel,TF=tf,score=typ,species=sp,activity=y[i]))
  for group,gg in groups(panel).items():
   gg=sorted(set(gg)&set(Z.index))
   if len(gg)<2:continue
   y=Z.loc[gg].mean(axis=0).to_numpy()-p['background'];tfgroup.append(dict(panel=panel,TF=tf,group=group,n_group=len(gg),n_targets=len(aa),**fit(np.asarray(sc),y,cx,extra=np.column_stack([mass,p['background']]),loo=True)))
 # Evidence for edges does not treat 22 libraries as 22 independent species.
 for row in prior.itertuples():
  if row.source not in Z.index or row.target not in Z.index or row.target in story|{row.source}:continue
  f=fit(Z.loc[row.source].to_numpy(),Z.loc[row.target].to_numpy(),cx,extra=np.column_stack([mass,p['background']]),loo=True)
  edges.append(dict(panel=panel,TF=row.source,target=row.target,prior_sign=float(np.sign(row.weight)),**f))
 for mode in ['raw','sample_centered']:
  wd=O/'results/wgcna'/f'{panel}__{mode}';members=pd.read_csv(wd/'members.tsv',sep='\t');e=pd.read_csv(wd/'eigengenes.tsv',sep='\t').merge(p['meta'][['sample','species']],on='sample').groupby(['species','module']).eigengene.mean().unstack().loc[p['sp']]
  for module,mm in members[members.module.ne('grey')].groupby('module'):
   fs=mm.feature.tolist();y=e[module].to_numpy();module_mass.append(dict(panel=panel,mode=mode,module=module,n_features=len(mm),variance_PC1=mm.variance_PC1.iloc[0],**fit(mass,y,cx,loo=True),ols_p=fit(mass,y,np.eye(len(mass)))['p']))
   modgenes=syms(p,fs);universe=syms(p,members.feature)
   for tf,aa in prior.groupby('source'):
    target=set(aa.target)&universe;k=len(modgenes&target)
    if len(target)<10:continue
    table=[[k,len(modgenes)-k],[len(target)-k,len(universe-modgenes-target)]]
    odds,pv=stats.fisher_exact(table,alternative='greater')
    module_tf.append(dict(panel=panel,mode=mode,module=module,TF=tf,n_module_symbols=len(modgenes),n_targets=len(target),overlap=k,odds=odds,p=pv,overlap_genes=';'.join(sorted(modgenes&target))))
   # Remove ALL own-story genes from the module score before testing any story association.
   keep=mm[~mm.feature.isin(group_features(p,story))]
   if len(keep)>=20:
    x=p['X'].loc[keep.feature]
    if mode=='sample_centered':x=x-p['X'].mean(axis=0)
    z=x.sub(x.mean(axis=1),axis=0).div(x.std(axis=1),axis=0);sc=keep.weight.to_numpy()@z.to_numpy()
    for group,gg in groups(panel).items():
     ff=group_features(p,gg)
     if len(ff)<2:continue
     gx=p['X'].loc[ff];gz=gx.sub(gx.mean(axis=1),axis=0).div(gx.std(axis=1),axis=0).mean(axis=0).to_numpy()
     module_groups.append(dict(panel=panel,mode=mode,module=module,group=group,n_group=len(ff),n_module_no_story=len(keep),**fit(sc,gz,cx,extra=mass,loo=True)))
  if panel in ['ma2016_full']:
   validations=['gse105160_full','wagner_full']
   for vp in validations:
    v=P[vp]
    for module,mm in members[members.module.ne('grey')].groupby('module'):
     m0=mm[mm.feature.isin(v['X'].index)];cover=len(m0)/len(mm)
     if len(m0)<20 or cover<.8:continue
     x=v['X'].loc[m0.feature]
     if mode=='sample_centered':x=x-v['X'].mean(axis=0)
     z=x.sub(x.mean(axis=1),axis=0).div(x.std(axis=1),axis=0)
     y=m0.weight.to_numpy()@z.to_numpy();cors=[np.corrcoef(zz,y)[0,1] for zz in z.to_numpy()]
     projection.append(dict(training=panel,validation=vp,mode=mode,module=module,coverage=cover,n_features=len(m0),median_gene_score_cor=float(np.median(cors)),**fit(v['mass'],y,v['C'],loo=True),ols_p=fit(v['mass'],y,np.eye(len(y)))['p']))
E=save('regulon_edges.tsv.gz',edges)
# Learn only in discovery panels. q correction includes ALL candidate edges in all six panels.
for train,valids in [('ma2016_full',['ma2016_full','gse105160_full','wagner_full'])]:
 for criterion in ['q10','nominal05']:
  ok=E.q_all_tests.le(.1) if criterion=='q10' else E.p.le(.05)
  es=E[E.panel.eq(train)&ok&E.loo_sign_fraction.ge(.8)]
  for tf,aa in es.groupby('TF'):
   if len(aa)<10:continue
   frozen=aa[['target','beta','prior_sign']].copy();frozen['learned_sign']=np.sign(frozen.beta)
   for rr in frozen.to_dict('records'):learnedmembership.append(dict(training=train,criterion=criterion,TF=tf,**rr))
   for vp in valids:
    v=P[vp];a=frozen[frozen.target.isin(v['Z'].index)];cover=len(a)/len(frozen)
    if len(a)<10 or cover<.8:continue
    w=a.learned_sign.to_numpy();y=w@v['Z'].loc[a.target].sub(v['background'],axis=1).to_numpy()/len(a)
    learnedmass.append(dict(training=train,validation=vp,criterion=criterion,TF=tf,n_targets=len(a),coverage=cover,prior_sign_agreement=float(np.mean(a.learned_sign==a.prior_sign)),**fit(v['mass'],y,v['C'],loo=True)))
for name,rows in [('module_mass.tsv',module_mass),('module_fixed_groups_no_overlap.tsv',module_groups),('module_curated_TF_enrichment.tsv',module_tf),('module_projection_mass.tsv',projection),('prior_TF_mass.tsv',tfmass),('prior_TF_fixed_groups.tsv',tfgroup),('learned_regulon_mass.tsv',learnedmass)]:
 d=save(name,rows);print(name,len(d),'q05',int(d.q_all_tests.lt(.05).sum()) if 'q_all_tests' in d else 0,flush=True)
pd.DataFrame(coverage).to_csv(O/'results/full_regulon_coverage.tsv',sep='\t',index=False);pd.DataFrame(scores).to_csv(O/'results/full_prior_TF_activity.tsv',sep='\t',index=False);pd.DataFrame(learnedmembership).to_csv(O/'results/full_learned_regulon_membership.tsv',sep='\t',index=False)
(O/'audit/statistics_full_execution.json').write_text(json.dumps(dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),prior_sha256=hashlib.sha256(priorpath.read_bytes()).hexdigest(),units='species, phylogenetic covariance fixed lambda1, lambda0 sensitivity, leave-one-species-out; no random-effect donor inference',WGCNA='traitblind networks; below15 independent species exploratory; transfer validates frozen membership/weights with within-panel gene standardization',edges='partial TF-target association given mass and global standardized-expression mean; neither direction nor binding is inferred',human_review='not_completed',training_enabled=False),indent=2))
