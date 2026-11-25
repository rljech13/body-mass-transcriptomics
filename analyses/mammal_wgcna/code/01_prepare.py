from pathlib import Path
import os,sys,json,hashlib,re,collections
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd,scipy
O=Path(__file__).resolve().parents[1];R=O.parent;H=Path('/path/to/user-home/chicken/kleiber_birds_gemcat');B=R/'sizes6_three_gaps_2026-10-19_v1';S=R/'sizes4_execution_2026-10-16_v3';V=R/'nonbird_regulatory_followup_2026-10-21_v1';A=R/'nonbird_common_state_2026-10-21_v1'
def use(p):
 p=Path(p);inputs.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size));return p
def read(p,**kw):return pd.read_csv(use(p),sep='\t',**kw)
inputs=[];manifest=[];cov={}
birdsets=json.loads(use(S/'05_regulation/program_memberships.json').read_text());birdsets={k:v for k,v in birdsets.items() if not k.startswith('EXPLORATORY')}
mammalsets=json.loads(use(R/'nonbird_decision_2026-10-21_v1/config/memberships.json').read_text())
(O/'config/fixed_sets.json').write_text(json.dumps({'bird':birdsets,'mammal':mammalsets},indent=2))
def emit(panel,M,meta,features,C):
 M=M.astype(float);assert M.index.is_unique and M.columns.is_unique
 assert np.isfinite(M.to_numpy()).all() and (M.to_numpy()>=0).all()
 M=M.loc[(M.ge(1).sum(1)>=max(2,int(np.ceil(M.shape[1]*.5)))) & (np.log2(M+1).std(1)>1e-8)]
 # No trait-derived filtering or ranking.
 M.to_csv(O/'inputs'/f'{panel}.linear.tsv.gz',sep='\t');np.log2(M+1).to_csv(O/'inputs'/f'{panel}.log.tsv',sep='\t')
 meta=meta.set_index('sample').loc[M.columns].reset_index();meta.to_csv(O/'inputs'/f'{panel}.samples.tsv',sep='\t',index=False)
 features=features.set_index('feature').loc[M.index].reset_index();features.to_csv(O/'inputs'/f'{panel}.features.tsv',sep='\t',index=False)
 sp=list(dict.fromkeys(meta.species));cov[panel]=C
 assert C.shape==(len(sp),len(sp)) and np.linalg.eigvalsh(C).min()>0
 np.savez_compressed(O/'inputs'/f'{panel}.cov.npz',C=C,species=np.array(sp))
 manifest.append(dict(panel=panel,n_features=len(M),n_libraries=M.shape[1],n_species=len(sp),network_unit='source_library' if panel.startswith('bird') else 'species_mean',trait_unit='species',physical_donor='unresolved_or_species_mean',log_transform='log2(1+normalized_linear)',features_filtered_before_trait_test=True))
tr=read(A/'config/species_trait_admission.corrected.tsv');tr=tr[tr.status.eq('admitted')]
tz=np.load(use(R/'sizes5_execution_2026-10-16_v1/01_mammals/existing_phylogenetic/dated_covariances.npz'));ts=list(tz['species'])
ms={p:read(V/'results'/f'{p}_admitted_M.tsv.gz',index_col=0) if p=='ma2016' else read(A/'results'/f'{p}_admission_before_gene_filter_M.tsv.gz',index_col=0) for p in ['ma2016','gse105160','wagner']}
common=sorted(set.intersection(*(set(m.index) for m in ms.values())))
for panel,M in ms.items():
 t=tr[tr.panel.eq(panel)];sp=t.label.tolist();M=M.loc[common,sp];meta=pd.DataFrame(dict(sample=sp,species=sp,scientific_name=t.scientific_name.tolist(),mass_g=t.accepted_mass_g.to_numpy(),context='adult_fibroblast',sex='unresolved',code=sp))
 feat=pd.DataFrame(dict(feature=M.index,symbols=M.index,primary_symbol=M.index));ix=[ts.index(x) for x in t.scientific_name];C=tz['covariances'][0][np.ix_(ix,ix)];emit(panel,M,meta,feat,C)
reg=read(B/'02_measurements/sample_registry.tsv');M=read(B/'02_measurements/native_complete_family7_sample_expression.tsv',index_col=0);ann=read(B/'02_measurements/native_complete_family7_mapping.tsv',dtype=str)
ann=ann.rename(columns={'orthogroup':'feature','chicken_symbols':'symbols'});ann['primary_symbol']=ann.symbols.where(~ann.symbols.str.contains(';'),ann.feature)
bt=read(S/'02_panels_traits/traits_long.tsv');bt=bt[bt.trait_id.eq('adult_mass_g')];bm=bt.set_index('species').value.to_dict();sci=bt.set_index('species').species_id.to_dict()
zz=np.load(use(S/'02_panels_traits/dated_covariances.npz'));tree=list(zz['species'])
for panel,context in [('bird_HG','fib_high_glucose'),('bird_F12','fib_F12'),('bird_E_label','epi_source_label')]:
 t=reg[reg.context.eq(context)].copy();t['sample']=t.species+'|'+t.column;t['mass_g']=t.species.map(bm);t['scientific_name']=t.species.map(sci);t['sex']='unresolved';sp=list(dict.fromkeys(t.species));ix=[tree.index(sci[x]) for x in sp];C=zz['covariances'][0][np.ix_(ix,ix)]
 emit(panel,M[t['sample']],t[['sample','species','scientific_name','mass_g','context','sex','code']],ann[['feature','symbols','primary_symbol']],C)
pd.DataFrame(manifest).to_csv(O/'config/panel_inventory.tsv',sep='\t',index=False)
(O/'audit/preparation.json').write_text(json.dumps(dict(inputs=inputs,panels=manifest,mammal_common_before_filter=len(common),bird_native_scope='complete-family7 including reciprocal within-species annotation bridge; single-copy motif scope separately',excluded='PGC C37/C38 already excluded by accepted source sample registry',E_label='anatomy and physical donors unresolved; not adult comparison',python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2))
print(pd.DataFrame(manifest).to_string(index=False),flush=True)
