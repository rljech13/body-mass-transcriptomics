from common import *
treefiles={'mammal':R/'sizes5_execution_2026-10-16_v1/01_mammals/existing_phylogenetic/dated_covariances.npz','bird':R/'sizes4_execution_2026-10-16_v3/02_panels_traits/dated_covariances.npz'}
def manyfit(x,Y,C):
 X=np.column_stack([np.ones(len(x)),x]);L=np.linalg.cholesky(C);X=linalg.solve_triangular(L,X,lower=True);Y=linalg.solve_triangular(L,Y,lower=True);B=np.linalg.lstsq(X,Y,rcond=None)[0];err=Y-X@B;df=len(x)-2;se=np.sqrt(np.sum(err*err,axis=0)/df*np.linalg.inv(X.T@X)[1,1]);tt=B[1]/np.maximum(se,1e-30);return B[1],2*stats.t.sf(abs(tt),df)
out=[]
for scope,invname,actname in [('common','panel_inventory.tsv','prior_TF_activity.tsv'),('full','full_inventory.tsv','full_prior_TF_activity.tsv')]:
 inv=pd.read_csv(O/'config'/invname,sep='\t');activity=pd.read_csv(O/'results'/actname,sep='\t')
 for panel in inv.panel:
  p=load(panel);clade='bird' if panel.startswith('bird') else 'mammal';t=np.load(treefiles[clade]);tree_sp=list(t['species']);meta=p['meta'].drop_duplicates('species').set_index('species').loc[p['sp']];ix=[tree_sp.index(x) for x in meta.scientific_name];Y=[];labels=[]
  for mode in ['raw','sample_centered']:
   e=pd.read_csv(O/'results/wgcna'/f'{panel}__{mode}'/'eigengenes.tsv',sep='\t').merge(p['meta'][['sample','species']],on='sample').groupby(['species','module']).eigengene.mean().unstack().loc[p['sp']]
   for module in e.columns:
    if module=='grey':continue
    Y.append(e[module].to_numpy());labels.append(dict(kind='module',mode=mode,name=module))
  a=activity[activity.panel.eq(panel)]
  for (tf,score),aa in a.groupby(['TF','score']):
   y=aa.set_index('species').loc[p['sp']].activity.to_numpy();Y.append(y);labels.append(dict(kind='prior_TF',mode=score,name=tf))
  Y=np.array(Y).T
  for k,C in enumerate(t['covariances'][:,ix][:,:,ix]):
   b,pv=manyfit(p['mass'],Y,C)
   for j,label in enumerate(labels):out.append(dict(scope=scope,panel=panel,tree=k,beta=b[j],p=pv[j],**label))
d=pd.DataFrame(out);d['q_tree_all_tests']=np.nan
for _,ix in d.groupby(['scope','kind','tree']).groups.items():d.loc[ix,'q_tree_all_tests']=bh(d.loc[ix,'p'])
d.to_csv(O/'results/tree_sensitivity_long.tsv.gz',sep='\t',index=False)
rr=[]
for key,a in d.groupby(['scope','panel','kind','mode','name']):
 b=a.loc[a.tree.eq(0)].beta.iloc[0];rr.append(dict(zip(['scope','panel','kind','mode','name'],key),n_trees=len(a),p_median=a.p.median(),p_min=a.p.min(),p_max=a.p.max(),q_median=a.q_tree_all_tests.median(),q_max=a.q_tree_all_tests.max(),same_sign_fraction=np.mean(np.sign(a.beta)==np.sign(b)),fraction_q05=np.mean(a.q_tree_all_tests<.05)))
pd.DataFrame(rr).to_csv(O/'results/tree_sensitivity_summary.tsv',sep='\t',index=False)
(O/'audit/tree_sensitivity_execution.json').write_text(json.dumps(dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),n_trees=100,scope='same frozen module and prior-regulon scores under all available dated covariance trees; uncertainty of species sampling or stage is not removed',input_sha256={k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in treefiles.items()}),indent=2));print('TREE DONE',len(d),flush=True)
