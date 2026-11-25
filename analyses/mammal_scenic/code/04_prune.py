from pathlib import Path
import os,json,pickle,time,hashlib,argparse
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
np.object=object
import pandas as pd
from ctxcore.rnkdb import FeatherRankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df,df2regulons
from pyscenic.aucell import aucell
from threadpoolctl import threadpool_limits
threadpool_limits(1)
R=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,required=True);args=ap.parse_args();seed=args.seed;O=R/'results'/str(seed);expr=pd.read_pickle(R/'inputs/expression.pkl');adj=pd.read_csv(O/'adjacencies.tsv.gz',sep='\t');t=time.time()
 modules=list(modules_from_adjacencies(adj,expr,rho_mask_dropouts=False));pd.DataFrame([dict(module=i,TF=m.transcription_factor,genes=len(m),context=';'.join(sorted(m.context))) for i,m in enumerate(modules)]).to_csv(O/'coexpression_modules.tsv',sep='\t',index=False)
 with (O/'coexpression_modules.pkl').open('wb') as f:pickle.dump(modules,f)
 print('MODULES',seed,len(modules),flush=True);dfs=[];windows={};summ=[]
 for tag in ['500up100down','10kb_each_side']:
  p=R/'reference'/f'human_{tag}.regions_vs_motifs.rankings.feather';db=FeatherRankingDatabase(str(p),name='human_GRCh38_JASPAR2024_'+tag)
  info=dict(seed=seed,window=tag,ranking_genes=db.total_genes,ranking_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),rank_threshold=1500,auc_threshold=.05,NES_threshold=3.,num_workers=16,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
  (O/f'prune_{tag}_before.json').write_text(json.dumps(info,indent=2));df=prune2df([db],modules,str(R/'reference/motif_annotations.tsv'),rank_threshold=1500,auc_threshold=.05,nes_threshold=3.,num_workers=16,client_or_address='custom_multiprocessing',module_chunksize=100,filter_for_annotation=True)
  df.to_pickle(O/f'motif_enrichment_{tag}.pkl');df.to_csv(O/f'motif_enrichment_{tag}.csv');regs=df2regulons(df) if len(df) else []
  with (O/f'regulons_{tag}.pkl').open('wb') as f:pickle.dump(regs,f)
  windows[tag]={r.name:set(r.gene2weight)-{r.transcription_factor} for r in regs};dfs.append(df);summ.append(dict(window=tag,regulons=len(regs),motif_enrichment_rows=len(df)));print('PRUNED',seed,tag,len(regs),flush=True)
 df=pd.concat(dfs);df.to_pickle(O/'motif_enrichment_union.pkl');regs=df2regulons(df) if len(df) else []
 with (O/'regulons.pkl').open('wb') as f:pickle.dump(regs,f)
 edges=[];inv=[]
 for r in regs:
  tf=r.transcription_factor;targets=set(r.gene2weight)-{tf};inv.append(dict(regulon=r.name,TF=tf,targets_nonself=len(targets),windows=sum(r.name in v for v in windows.values()),targets_both_windows=len(set.intersection(*[v.get(r.name,set()) for v in windows.values()]))))
  for g in sorted(targets):edges.append(dict(regulon=r.name,TF=tf,target=g,importance=r.gene2weight[g],window_count=sum(g in v.get(r.name,set()) for v in windows.values())))
 pd.DataFrame(inv,columns=['regulon','TF','targets_nonself','windows','targets_both_windows']).to_csv(O/'regulon_inventory.tsv',sep='\t',index=False);pd.DataFrame(edges,columns=['regulon','TF','target','importance','window_count']).to_csv(O/'regulon_edges.tsv.gz',sep='\t',index=False)
 if regs:
  auc=aucell(expr,regs,auc_threshold=.05,num_workers=8,seed=seed);auc.to_csv(O/'AUCell.tsv',sep='\t')
 (O/'SCENIC_summary.json').write_text(json.dumps(dict(seed=seed,coexpression_modules=len(modules),windows=summ,union_regulons=len(regs),nonself_edges=len(edges),seconds=time.time()-t,expression='26 bulk donor-labelled groups; not single cells',reference='custom human GRCh38+JASPAR2024; not native species occupancy'),indent=2));print('SCENIC_DONE',seed,len(regs),len(edges),time.time()-t,flush=True)
if __name__=='__main__':main()
