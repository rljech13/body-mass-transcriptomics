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
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<22),b''):h.update(b)
 return h.hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,required=True);args=ap.parse_args();seed=args.seed;O=R/'results'/str(seed);expr=pd.read_pickle(R/'inputs/expression.pkl');adj=pd.read_csv(O/'adjacencies.tsv.gz',sep='\t');t=time.time();cfg=json.loads((R/'config/existing_reference.json').read_text())
 modules=list(modules_from_adjacencies(adj,expr,rho_mask_dropouts=False));pd.DataFrame([dict(module=i,TF=m.transcription_factor,genes=len(m),context=';'.join(sorted(m.context))) for i,m in enumerate(modules)]).to_csv(O/'coexpression_modules.tsv',sep='\t',index=False)
 with (O/'coexpression_modules.pkl').open('wb') as f:pickle.dump(modules,f)
 print('MODULES',seed,len(modules),flush=True)
 p=Path(cfg['database']);ann=Path(cfg['annotation']);assert sha(p)==cfg['database_sha256'];assert sha(ann)==cfg['annotation_sha256'];db=FeatherRankingDatabase(str(p),name='human_hg38_v10_clustered_500up100down_existing')
 info=dict(seed=seed,ranking_genes=db.total_genes,ranking_sha256=cfg['database_sha256'],annotation_sha256=cfg['annotation_sha256'],rank_threshold=1500,auc_threshold=.05,NES_threshold=3.,num_workers=16,code_sha256=sha(Path(__file__)),reference='Existing official hg38 full-tx 500bpUp100down v10-clustered; all5876motifs match HGNC snapshot m1e-5')
 (O/'prune_before.json').write_text(json.dumps(info,indent=2));df=prune2df([db],modules,str(ann),rank_threshold=1500,auc_threshold=.05,nes_threshold=3.,num_workers=16,client_or_address='custom_multiprocessing',module_chunksize=100,filter_for_annotation=True)
 df.to_pickle(O/'motif_enrichment.pkl');df.to_csv(O/'motif_enrichment.csv');regs=df2regulons(df) if len(df) else []
 with (O/'regulons.pkl').open('wb') as f:pickle.dump(regs,f)
 edges=[];inv=[]
 for r in regs:
  tf=r.transcription_factor;targets=set(r.gene2weight)-{tf};inv.append(dict(regulon=r.name,TF=tf,targets_nonself=len(targets)))
  for g in sorted(targets):edges.append(dict(regulon=r.name,TF=tf,target=g,importance=r.gene2weight[g]))
 pd.DataFrame(inv,columns=['regulon','TF','targets_nonself']).to_csv(O/'regulon_inventory.tsv',sep='\t',index=False);pd.DataFrame(edges,columns=['regulon','TF','target','importance']).to_csv(O/'regulon_edges.tsv.gz',sep='\t',index=False)
 if regs:
  auc=aucell(expr,regs,auc_threshold=.05,num_workers=8,seed=seed);auc.to_csv(O/'AUCell.tsv',sep='\t')
 (O/'SCENIC_summary.json').write_text(json.dumps(dict(seed=seed,coexpression_modules=len(modules),motif_enrichment_rows=len(df),union_regulons=len(regs),nonself_edges=len(edges),seconds=time.time()-t,expression='26 bulk donor-labelled groups; not single cells',reference='existing human hg38 v10-clustered short-promoter DB; not native species occupancy',two_window_sensitivity=False),indent=2));print('SCENIC_DONE',seed,len(regs),len(edges),time.time()-t,flush=True)
if __name__=='__main__':main()
