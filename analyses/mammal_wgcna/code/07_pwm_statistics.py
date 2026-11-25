from pathlib import Path
import sys,json,hashlib,itertools
import numpy as np,pandas as pd
from scipy import stats,linalg
from common import O,R,H,bh,fit,load
mode=sys.argv[1];meta_path=O/'reference/bird_native_promoters.tsv' if mode=='bird' else O/'reference/human/promoters.tsv';fimo_path=O/'reference/bird_candidate_fimo.tsv' if mode=='bird' else O/'reference/human/fimo.tsv'
if mode=='human_full':meta_path=O/'reference/human_plus_own5/promoters.tsv';fimo_path=O/'reference/human_plus_own5/fimo.tsv'
T=pd.read_csv(meta_path,sep='\t');F=pd.read_csv(fimo_path,sep='\t',comment='#');motifs=pd.read_csv(O/'reference/candidate_motifs.tsv',sep='\t');sets=json.loads((O/'config/fixed_sets.json').read_text());TF=json.loads((O/'config/analysis_plan.initial.json').read_text())['candidate_TFs']
B=R/'sizes6_three_gaps_2026-10-19_v1';ann=pd.read_csv(B/'02_measurements/native_ortholog7_mapping.tsv',sep='\t',dtype=str);og2sym=ann.set_index('orthogroup').symbol.to_dict()
panels=['bird_HG','bird_F12','bird_E_label'] if mode=='bird' else ['ma2016','gse105160','wagner'];panels=['ma2016_full','gse105160_full','wagner_full'] if mode=='human_full' else panels;contexts={'bird_HG':'fib_high_glucose','bird_F12':'fib_F12','bird_E_label':'epi_source_label'}
rows=[];edge_rows=[];matrices=[];coverage=[]
def exact_tail(labels,hits,cases):
 dist=np.array([1.]);expected=0.
 for st in np.unique(labels):
  ix=np.flatnonzero(labels==st);N=len(ix);K=int(hits[ix].sum());n=int(cases[ix].sum());expected+=n*K/N
  if n:
   k=np.arange(n+1);dist=np.convolve(dist,stats.hypergeom.pmf(k,N,K,n))
 obs=int(hits[cases].sum());return obs,expected,float(dist[obs:].sum())
def residual_edges(x,Y,C,mass,bg):
 Q=np.column_stack([np.ones(len(x)),mass,bg]);L=np.linalg.cholesky(C);Q=linalg.solve_triangular(L,Q,lower=True);x=linalg.solve_triangular(L,x,lower=True);Y=linalg.solve_triangular(L,Y,lower=True)
 x=x-Q@np.linalg.lstsq(Q,x,rcond=None)[0];Y=Y-Q@np.linalg.lstsq(Q,Y,rcond=None)[0];r=(x@Y)/(np.linalg.norm(x)*np.linalg.norm(Y,axis=0));r=np.clip(r,-1+1e-14,1-1e-14);df=len(x)-np.linalg.matrix_rank(Q)-1;tt=r*np.sqrt(df/(1-r*r));return r,2*stats.t.sf(abs(tt),df)
for panel in panels:
 print('PWM',panel,flush=True);p=load(panel);groups=sets['bird' if mode=='bird' else 'mammal'];t=T[T.context.eq(contexts[panel])].copy() if mode=='bird' else T.copy()
 if mode=='bird':
  nsp=p['meta'].species.nunique();t=t[t.species.isin(p['sp'])];genes=sorted(set(t.groupby('gene').species.nunique().loc[lambda z:z.eq(nsp)].index)&set(p['G'].index))
 else:nsp=1;genes=sorted(set(t.gene)&set(p['G'].index))
 t=t[t.gene.isin(genes)];cov=t.groupby('gene')[['GC','CpG']].median().reindex(genes);cov['expression']=p['G'].reindex(genes).mean(axis=1)
 labels=np.zeros(len(genes),int)
 for j,col in enumerate(['GC','CpG','expression']):labels+=(cov[col].to_numpy()>cov[col].median()).astype(int)*(2**j)
 zcov=(cov-cov.mean())/cov.std();g2i={g:i for i,g in enumerate(genes)};seq2gene=t.set_index('sequence_name').gene.to_dict();names=set(seq2gene)
 feature2symbol={f:og2sym[f] for f in p['X'].index if f in og2sym} if mode=='bird' else {f:f for f in p['X'].index}
 tests=[]
 for group,gg in groups.items():
  seen=sorted(set(gg)&set(genes));tests.append(dict(kind='fixed_group',group=group,mode='fixed',requested=len(gg),genes=seen,complete=set(gg)<=set(genes)))
 for network_mode in ['raw','sample_centered']:
  members=pd.read_csv(O/'results/wgcna'/f'{panel}__{network_mode}'/'members.tsv',sep='\t')
  for module,mm in members[members.module.ne('grey')].groupby('module'):
   gg=sorted(set(feature2symbol.get(f) for f in mm.feature)-{None});seen=sorted(set(gg)&set(genes))
   tests.append(dict(kind='WGCNA_module',group=module,mode=network_mode,requested=len(mm),genes=seen,complete=len(seen)>=.8*len(mm)))
 for threshold in [1e-4,1e-5]:
  hit=np.zeros((len(genes),len(motifs)),dtype=np.int16);mi={m:i for i,m in enumerate(motifs.motif_id)}
  ff=F[F.sequence_name.isin(names)&F['p-value'].le(threshold)].drop_duplicates(['sequence_name','motif_id'])
  for rr in ff.itertuples(index=False):hit[g2i[seq2gene[rr.sequence_name]],mi[rr.motif_id]]+=1
  binary=hit>=int(np.ceil(nsp/2));np.savez_compressed(O/'results'/f'{panel}.PWM.{threshold}.npz',genes=np.array(genes),motif_ids=motifs.motif_id.to_numpy(dtype=str),counts=hit,H=binary,labels=labels)
  for test in tests:
   case=np.array([g in set(test['genes']) for g in genes]);n=int(case.sum());minimum=2 if test['kind']=='fixed_group' else 20
   if n<minimum:
    coverage.append(dict(panel=panel,kind=test['kind'],group=test['group'],mode=test['mode'],n_requested=test['requested'],n_covered=n,complete=test['complete'],status='insufficient_sequence_scope'));continue
   # Residual matching imbalance against the exact stratified-null expectation.
   expect=np.zeros(3)
   for st in np.unique(labels):
    ix=labels==st;expect+=case[ix].sum()*zcov.to_numpy()[ix].mean(axis=0)/n
   balance=float(np.max(np.abs(zcov.to_numpy()[case].mean(axis=0)-expect)))
   for j,mo in enumerate(motifs.itertuples()):
    obs,ex,pv=exact_tail(labels,binary[:,j],case);rows.append(dict(panel=panel,kind=test['kind'],group=test['group'],mode=test['mode'],threshold=threshold,motif_id=mo.motif_id,TF=mo.TF,n_universe=len(genes),n_requested=test['requested'],n_covered=n,complete=test['complete'],n_species=nsp,min_species_with_site=int(np.ceil(nsp/2)),observed=obs,expected=ex,p=pv,matching_max_abs_SMD=balance,hit_genes=';'.join(np.array(genes)[case&binary[:,j]])))
  if threshold!=1e-4:continue
  # Genome-wide motif-pruned putative edges, adjusted for mass and global expression.
  story=set(g for vv in groups.values() for g in vv)
  for tf in TF:
   if tf not in p['Z'].index:continue
   jj=[j for j,tt in enumerate(motifs.TF) if tf in tt.split('::')]
   if not jj:continue
   target=np.array(genes)[binary[:,jj].any(axis=1)];target=[g for g in target if g in p['Z'].index and g!=tf]
   if not target:continue
   x=p['Z'].loc[tf].to_numpy();Y=p['Z'].loc[target].to_numpy().T;rr,pv=residual_edges(x,Y,p['C'],p['mass'],p['background'])
   loo=[]
   for k in range(len(x)):
    ix=np.arange(len(x))!=k;r0,p0=residual_edges(x[ix],Y[ix],p['C'][np.ix_(ix,ix)],p['mass'][ix],p['background'][ix]);loo.append(r0)
   stability=(np.sign(np.array(loo))==np.sign(rr)).mean(axis=0)
   for j,g in enumerate(target):edge_rows.append(dict(panel=panel,TF=tf,target=g,partial_r=float(rr[j]),p=float(pv[j]),loo_sign_fraction=float(stability[j]),target_in_fixed_story=g in story,n_species=len(x),motif_support='native_CAGE_cluster_majority_species' if mode=='bird' else 'human_reference_promoter',binding_demonstrated=False))
D=pd.DataFrame(rows);D['q_mode_all_tests']=bh(D.p);D.to_csv(O/'results'/f'PWM_{mode}_enrichment.tsv.gz',sep='\t',index=False)
E=pd.DataFrame(edge_rows)
if len(E):E['q_mode_all_edges']=bh(E.p)
E.to_csv(O/'results'/f'motif_regulon_{mode}_edges.tsv.gz',sep='\t',index=False);pd.DataFrame(coverage).drop_duplicates().to_csv(O/'results'/f'PWM_{mode}_coverage_gaps.tsv',sep='\t',index=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
(O/'audit'/f'PWM_{mode}_execution.json').write_text(json.dumps(dict(script_sha256=sha(Path(__file__)),metadata_sha256=sha(meta_path),fimo_sha256=sha(fimo_path),n_tests=len(D),n_edges=len(E),matching='conditional hypergeometric tail convolved over 8 GC/CpG/expression median strata; cases retained in population; coarse matching assumption, residual SMD exported',null_scope='gene-set enrichment; motif occurrence aggregated across species, not independent species tests',motif_specificity='shared monomer/heterodimer PWMs do not identify unique TF or demonstrate binding',human_review='not_completed'),indent=2))
print('DONE',mode,'tests',len(D),'edges',len(E),'minq',D.q_mode_all_tests.min(),flush=True)
