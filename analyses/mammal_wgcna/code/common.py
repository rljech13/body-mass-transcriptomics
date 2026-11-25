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
 d.to_csv(O/'results'/name,sep='\t',index=False);return d
