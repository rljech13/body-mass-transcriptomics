import pathlib,json,time,hashlib,inspect,multiprocessing as mp,os,argparse
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
np.object=object
import pandas as pd
from arboreto.core import infer_partial_network,SGBM_KWARGS,EARLY_STOP_WINDOW_LENGTH,to_tf_matrix
from arboreto.algo import _prepare_input
from threadpoolctl import threadpool_limits
threadpool_limits(1)
R=pathlib.Path(__file__).resolve().parents[1]
def infer(i):
 return infer_partial_network(regressor_type='GBM',regressor_kwargs=SGBM_KWARGS,tf_matrix=TFX,tf_matrix_gene_names=TFN,target_gene_name=GN[i],target_gene_expression=X[:,i],include_meta=False,early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,seed=SEED)
def main():
 global X,TFX,TFN,GN,SEED
 ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,required=True);args=ap.parse_args();SEED=args.seed;out=R/'results'/str(SEED);out.mkdir(exist_ok=True)
 start=time.time();expr=pd.read_pickle(R/'inputs/expression.pkl');tfs=(R/'inputs/tf_names.txt').read_text().splitlines();X,GN,tfs=_prepare_input(expr,None,tfs);TFX,TFN=to_tf_matrix(X,GN,tfs)
 info=dict(seed=SEED,workers=24,n_profiles=X.shape[0],genes=X.shape[1],TFs=len(TFN),regressor_kwargs=SGBM_KWARGS,early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,code_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),input_sha256=hashlib.sha256((R/'inputs/expression.pkl').read_bytes()).hexdigest(),official_core='arboreto.core.infer_partial_network; GRNBoost2 defaults; fork-shared dispatcher')
 (out/'GRN_before.json').write_text(json.dumps(info,indent=2));(R/'audit/arboreto_infer_partial_network.py').write_text(inspect.getsource(infer_partial_network));print(json.dumps(info),flush=True)
 parts=[]
 with mp.get_context('fork').Pool(24) as pool:
  for i,a in enumerate(pool.imap_unordered(infer,range(len(GN)),chunksize=2),1):
   parts.append(a)
   if i%1000==0 or i==len(GN):print(json.dumps(dict(finished_targets=i,total=len(GN),elapsed_seconds=time.time()-start)),flush=True)
 adj=pd.concat(parts).sort_values(['importance','TF','target'],ascending=[False,True,True],kind='stable');assert not adj.duplicated(['TF','target']).any();adj.to_csv(out/'adjacencies.tsv.gz',sep='\t',index=False)
 (out/'GRN_summary.json').write_text(json.dumps(dict(positive_edges=len(adj),seconds=time.time()-start,seed=SEED,official_core_GRNBoost2=True),indent=2));print('GRN_DONE',SEED,len(adj),time.time()-start,flush=True)
if __name__=='__main__':main()
