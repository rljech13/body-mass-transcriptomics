"""Whole measured-map linear sums, no trait-driven membership selection."""
from pathlib import Path
import argparse, json, hashlib, sys, platform, itertools
import numpy as np, pandas as pd, scipy
from scipy import stats, sparse
ap=argparse.ArgumentParser();ap.add_argument('--results-root',type=Path,required=True);ap.add_argument('--genesets-root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
R,GS,O=a.results_root,a.genesets_root,a.out
for d in ['config','results','audit']:(O/d).mkdir(parents=True,exist_ok=True)
manifest=[]
def read(p,root,label,**kw):
 b=p.read_bytes();manifest.append(dict(origin=label,reference=str(p.relative_to(root)),sha256=hashlib.sha256(b).hexdigest(),bytes=len(b)))
 return pd.read_csv(p,sep='\t',**kw)
def bh(p):
 p=np.asarray(p,float);q=np.ones(len(p));ok=np.flatnonzero(np.isfinite(p));q[ok]=stats.false_discovery_control(p[ok]);return q
plan=dict(primary='log10(sum of all measured fixed-map genes in linear normalized RNA), PGLS lambda=1',
 catalogs='all KEGG and Reactome stored maps; all GO BP/MF/CC in stored MSigDB2026.1 C5; no Hallmark/HPO',
 admission='at least 5 measured genes and positive map sum in every admitted species; measured finite nonnegative genes, true zeros retained; no upper map-size limit',
 coverage='canonical gene count and measured fraction reported, no zero imputation; coverage>=50% secondary stratum',
 composition_modes=['panel_available','all3_common'],
 multiple_testing='BH within panel and across all panel-map tests separately per composition mode; same-sign max(two-sided p) conjunction tests across all3 or any panel pair, BH over all eligible terms/pairs; no independence assumption',
 sensitivity='ordinary unweighted Pearson/OLS on log10 scales; fixed diagonal covariance lambda0; leave-one-species for globally BH significant or conjunction nominal terms',
 null='1000 Brownian covariance simulated traits per panel; common simulated trait shared by every map within panel; descriptive hit-count calibration only; no cross-panel independence assumption',
 seed=2026102307,B=1000,trait='corrected species admission; Microtus pennsylvanicus 49g',human_review='not_completed',training_enabled=False)
(O/'config/plan.json').write_text(json.dumps(plan,indent=2))
cats=[]
for source,file in [('KEGG','kegg_hsa_all_pathways.tsv'),('Reactome','mammal_catalogs/reactome.tsv')]:
 t=read(GS/file,GS,'gene_sets')
 for row in t.itertuples():
  name=(str(row.pathway)+' '+row.name) if source=='KEGG' else row.gs_name
  cats.append(dict(source=source,term=name,genes=sorted(set(str(row.genes).split(';')))))
t=read(O/'config/go_all.tsv',O,'current_run')
for row in t.itertuples():cats.append(dict(source=row.collection.replace(':','_'),term=row.term,genes=sorted(set(str(row.genes).split(';')))))
catalog=pd.DataFrame(cats);catalog['set_id']=catalog.source+':'+catalog.term;assert catalog.set_id.is_unique
catalog['n_canonical']=catalog.genes.map(len)
catalog.assign(genes=catalog.genes.map(';'.join)).to_csv(O/'config/catalog.tsv.gz',sep='\t',index=False)
traits=read(R/'nonbird_common_state_2026-10-21_v1/config/species_trait_admission.corrected.tsv',R,'results');data={};panels=['ma2016','gse105160','wagner'];speciesrows=[]
for panel in panels:
 p=R/('nonbird_regulatory_followup_2026-10-21_v1/results/ma2016_admitted_M.tsv.gz' if panel=='ma2016' else f'nonbird_common_state_2026-10-21_v1/results/{panel}_admission_before_gene_filter_M.tsv.gz')
 M=read(p,R,'results',index_col=0)
 p=R/f'growth_regulators_2026-10-21_v1/inputs/{panel}_full.cov.npz';b=p.read_bytes();manifest.append(dict(origin='results',reference=str(p.relative_to(R)),sha256=hashlib.sha256(b).hexdigest(),bytes=len(b)))
 cv=np.load(p,allow_pickle=True);sp=list(map(str,cv['species']));tr=traits[traits.panel.eq(panel)&traits.status.eq('admitted')].set_index('label');keep=[s for s in sp if s in tr.index];ix=[sp.index(s) for s in keep];C=cv['C'][np.ix_(ix,ix)];C/=np.diag(C).mean();M=M[keep]
 good=np.isfinite(M).all(1)&M.ge(0).all(1);M=M.loc[good];assert M.index.is_unique
 x=np.log10(tr.loc[keep,'accepted_mass_g'].to_numpy(float));data[panel]=dict(M=M,C=C,x=x,species=keep)
 np.savez_compressed(O/'results'/f'{panel}_design.npz',C=C,mass_log10=x,species=keep)
 for i,s in enumerate(keep):speciesrows.append(dict(panel=panel,species=s,scientific_name=tr.loc[s,'scientific_name'],mass_g=10**x[i],sum_available=float(M[s].sum()),n_measured=len(M)))
pd.DataFrame(speciesrows).to_csv(O/'results/species.tsv',sep='\t',index=False)
common=set.intersection(*(set(d['M'].index) for d in data.values()))
def fit(Y,x,C):
 Q=np.linalg.inv(C);one=np.ones(len(x));P=Q-np.outer(Q@one,Q@one)/(one@Q@one);xx=x@P@x;xy=Y@P@x;yy=np.einsum('ij,ij->i',Y@P,Y);r=np.clip(xy/np.sqrt(np.maximum(xx*yy,1e-300)),-1+1e-14,1-1e-14);p=2*stats.t.sf(np.abs(r)*np.sqrt((len(x)-2)/(1-r*r)),len(x)-2)
 return xy/xx,p,r,P
allrows=[];scoreblocks={};cover=[];nullrows=[];rng=np.random.default_rng(plan['seed'])
for mode in plan['composition_modes']:
 for panel,d in data.items():
  M=d['M'] if mode=='panel_available' else d['M'].loc[sorted(common)];gindex={g:i for i,g in enumerate(M.index)};members=[sorted(set(g)&set(gindex)) for g in catalog.genes];nn=np.array(list(map(len,members)));ids=np.flatnonzero(nn>=5)
  rr=[];cc=[]
  for j in ids:
   rr.extend([j]*nn[j]);cc.extend(gindex[g] for g in members[j])
  A=sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(catalog),len(M)));sums=np.asarray(A@M.to_numpy());positive=np.isfinite(sums).all(1)&(sums>0).all(1);ids=ids[positive[ids]];Y=np.log10(sums[ids]);base=catalog.iloc[ids].drop(columns='genes').copy();base['mode']=mode;base['panel']=panel;base['n_species']=len(d['x']);base['n_measured']=nn[ids];base['coverage']=nn[ids]/base.n_canonical;base['genes']=[';'.join(members[i]) for i in ids]
  beta,p,r,P=fit(Y,d['x'],d['C']);base['beta']=beta;base['p']=p;base['r_phylo']=r;base['q_panel']=bh(p);b0,p0,r0,_=fit(Y,d['x'],np.eye(len(d['x'])));base['beta_OLS']=b0;base['p_Pearson']=p0;base['r_Pearson']=r0;base['q_Pearson_panel']=bh(p0);_,pdg,_,_=fit(Y,d['x'],np.diag(np.diag(d['C'])));base['p_diagonal']=pdg;allrows.append(base)
  testedset=set(ids)
  cover.extend(dict(mode=mode,panel=panel,set_id=catalog.set_id[i],n_canonical=int(catalog.n_canonical[i]),n_measured=int(nn[i]),status='tested' if i in testedset else 'fewer_than_5_measured_or_zero_sum') for i in range(len(catalog)))
  scoreblocks[mode,panel]=(ids,Y,members);np.savez_compressed(O/'results'/f'{mode}_{panel}_scores.npz',catalog_index=ids,log10_sum=Y)
  if mode=='panel_available':
   # Each null trait sees every map: dependencies from overlapping memberships are preserved.
   xs=rng.multivariate_normal(np.zeros(len(d['x'])),d['C'],size=plan['B']);xx=np.einsum('ij,ij->i',xs@P,xs);yy=np.einsum('ij,ij->i',Y@P,Y);nr=np.clip((Y@P@xs.T)/np.sqrt(yy[:,None]*xx),-1+1e-14,1-1e-14);ps=2*stats.t.sf(np.abs(nr)*np.sqrt((len(d['x'])-2)/(1-nr*nr)),len(d['x'])-2);nom=(ps<.05).sum(0);nbh=np.array([(bh(ps[:,i])<.05).sum() for i in range(plan['B'])]);np.savez_compressed(O/'results'/f'{panel}_null.npz',nominal_hits=nom,BH_hits=nbh)
   nullrows.append(dict(panel=panel,observed_nominal=int((p<.05).sum()),null_nominal_median=float(np.median(nom)),nominal_count_p=float((1+(nom>=(p<.05).sum()).sum())/(1+plan['B'])),observed_BH=int((bh(p)<.05).sum()),null_BH_median=float(np.median(nbh)),BH_count_p=float((1+(nbh>=(bh(p)<.05).sum()).sum())/(1+plan['B']))) )
  print(mode,panel,'genes',len(M),'maps',len(ids),'PGLS_BH',int((base.q_panel<.05).sum()),'Pearson_BH',int((base.q_Pearson_panel<.05).sum()),flush=True)
out=pd.concat(allrows,ignore_index=True);out['q_all_tests']=1.
for mode,z in out.groupby('mode'):out.loc[z.index,'q_all_tests']=bh(z.p)
out.to_csv(O/'results/map_tests.tsv.gz',sep='\t',index=False);pd.DataFrame(cover).to_csv(O/'results/coverage.tsv.gz',sep='\t',index=False);pd.DataFrame(nullrows).to_csv(O/'results/null_summary.tsv',sep='\t',index=False)
reps=[]
for mode,z in out.groupby('mode'):
 for pp in [panels]+list(map(list,itertools.combinations(panels,2))):
  pvt=z[z.panel.isin(pp)].pivot(index='set_id',columns='panel',values='p').dropna();bet=z[z.panel.isin(pp)].pivot(index='set_id',columns='panel',values='beta').loc[pvt.index,pp];same=(bet.gt(0).all(1)|bet.lt(0).all(1));pmax=pvt.max(1).where(same,1.)
  r=pd.DataFrame(dict(set_id=pvt.index,mode=mode,panels='+'.join(pp),same_sign=same.to_numpy(),p_rep=pmax.to_numpy()));r['q_rep']=bh(r.p_rep)
  for panel in pp:r[panel+'_p']=pvt[panel].to_numpy();r[panel+'_beta']=bet[panel].to_numpy()
  reps.append(r)
rep=pd.concat(reps,ignore_index=True);rep['q_all_pairs']=np.nan
for mode,z in rep[rep.panels.str.count('\\+').eq(1)].groupby('mode'):rep.loc[z.index,'q_all_pairs']=bh(z.p_rep)
rep.to_csv(O/'results/replication.tsv.gz',sep='\t',index=False)
# Fixed shortlist rule for diagnostics; no changes to memberships or primary tests.
diagids=set(out.loc[out.q_all_tests.lt(.05),'set_id'])|set(rep.loc[rep.p_rep.lt(.05),'set_id']);diagnostics=[];driver=[]
for mode in plan['composition_modes']:
 for panel,d in data.items():
  ids,Y,members=scoreblocks[mode,panel];lookup={catalog.set_id[i]:j for j,i in enumerate(ids)}
  for sid in sorted(diagids&set(lookup)):
   j=lookup[sid];i=ids[j];vals=[];ps=[]
   for k,s in enumerate(d['species']):
    mask=np.arange(len(d['x']))!=k;b,p,_,_=fit(Y[j:j+1,mask],d['x'][mask],d['C'][np.ix_(mask,mask)]);vals.append(float(b[0]));ps.append(float(p[0]))
   original=out[(out['mode']==mode)&out.panel.eq(panel)&out.set_id.eq(sid)].iloc[0]
   diagnostics.append(dict(mode=mode,panel=panel,set_id=sid,loo_sign_fraction=float(np.mean(np.sign(vals)==np.sign(original.beta))),loo_max_p=max(ps),loo_min_beta=min(vals),loo_max_beta=max(vals)))
   gs=members[i];gm=d['M'].loc[gs];share=gm.div(gm.sum(0),axis=1).mean(1).sort_values(ascending=False);driver.append(dict(mode=mode,panel=panel,set_id=sid,top_gene=share.index[0],top_gene_mean_share=float(share.iloc[0]),top3_mean_share=float(share.iloc[:3].sum())))
pd.DataFrame(diagnostics).to_csv(O/'results/shortlist_leave_one_species.tsv.gz',sep='\t',index=False);pd.DataFrame(driver).to_csv(O/'results/shortlist_gene_dominance.tsv.gz',sep='\t',index=False)
summary=dict(catalog_counts=catalog.source.value_counts().to_dict(),catalog_total=len(catalog),common_genes=len(common),panels=[],replication=[],limits=['Relative normalized RNA sums, not absolute RNA per cell or flux.','Human-reference catalog symbols use inherited mapping; not a new native orthology audit.','Partial measured map composition differs between panels in panel_available; all3_common fixes composition.','Small species counts, different source contexts and pathway overlap remain.','Post hoc whole-catalog screen, not independent validation of previously selected own5.'],human_review='not_completed',training_enabled=False)
for (mode,panel),z in out.groupby(['mode','panel']):summary['panels'].append(dict(mode=mode,panel=panel,n_tests=len(z),nominal=int(z.p.lt(.05).sum()),BH_panel=int(z.q_panel.lt(.05).sum()),BH_all_tests=int(z.q_all_tests.lt(.05).sum()),q_min=float(z.q_panel.min()),Pearson_BH=int(z.q_Pearson_panel.lt(.05).sum()),coverage_median=float(z.coverage.median()),maps_ge50pct=int(z.coverage.ge(.5).sum())))
for (mode,pp),z in rep.groupby(['mode','panels']):summary['replication'].append(dict(mode=mode,panels=pp,n_terms=len(z),nominal=int(z.p_rep.lt(.05).sum()),BH=int(z.q_rep.lt(.05).sum()),q_min=float(z.q_rep.min())))
(O/'results/summary.json').write_text(json.dumps(summary,indent=2));(O/'audit/input_manifest.json').write_text(json.dumps(manifest,indent=2));(O/'audit/environment.json').write_text(json.dumps(dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,platform=platform.platform()),indent=2));print(json.dumps(summary,indent=2),flush=True)
