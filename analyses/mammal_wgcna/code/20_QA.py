from common import *
import itertools,subprocess,datetime
checks=[]
def check(name,condition,detail=None):
 checks.append(dict(check=name,pass_=bool(condition),detail=detail));assert condition,(name,detail)
a=pd.read_csv(O/'audit/GLS_python.tsv',sep='\t').merge(pd.read_csv(O/'audit/GLS_baseR.tsv',sep='\t'),on='panel',suffixes=('_py','_R'))
check('independent_baseR_GLS',np.allclose(a.beta_py,a.beta_R,rtol=1e-10)&np.allclose(a.p_py,a.p_R,rtol=1e-10))
v=np.random.default_rng(20261021).uniform(size=1000);check('BH_matches_SciPy',np.allclose(bh(v),stats.false_discovery_control(v)))
# Exhaustive tiny conditional null, independent of sequence enrichment data.
labels=np.array([0,0,0,1,1,1]);hits=np.array([1,0,1,0,1,0]);cases=np.array([1,1,0,1,0,0],bool);observed=hits[cases].sum();vals=[]
for i in itertools.combinations(range(3),2):
 for j in range(3,6):vals.append(hits[list(i)+[j]].sum())
dist=np.convolve(stats.hypergeom.pmf(np.arange(3),3,2,2),stats.hypergeom.pmf(np.arange(2),3,1,1))
check('stratified_hypergeometric_equals_exhaustive',np.isclose(np.mean(np.array(vals)>=observed),dist[observed:].sum()))
def fasta(p):
 d={};key=None
 for l in p.read_text().splitlines():
  if l.startswith('>'):key=l[1:].split()[0];d[key]=''
  elif key:d[key]+=l.strip()
 return d
rng=np.random.default_rng(20261021);sitechecks=[];export=[]
for clade,fp,mp,fa in [('bird',O/'reference/bird_candidate_fimo.tsv',O/'reference/bird_native_promoters.tsv',O/'reference/bird_native_promoters.fa'),('human',O/'reference/human/fimo.tsv',O/'reference/human/promoters.tsv',O/'reference/human/promoters.fa')]:
 f=pd.read_csv(fp,sep='\t',comment='#');t=pd.read_csv(mp,sep='\t').set_index('sequence_name');seq=fasta(fa);ix=rng.choice(len(f),min(500,len(f)),replace=False)
 for row in f.iloc[ix].itertuples(index=False):
  s=seq[row.sequence_name][int(row.start)-1:int(row.stop)]
  if row.strand=='-':s=s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
  assert s.upper()==str(row.matched_sequence).upper()
 check(clade+'_FIMO_site_coordinates_and_strands',True,{'random_sites_checked':len(ix)})
 check(clade+'_FASTA_metadata_hashes',all(hashlib.sha256(seq[g].encode()).hexdigest()==t.loc[g,'sequence_sha256'] for g in t.index))
 # Evidence export for fixed groups; site sequences remain linked to exact genomic window metadata.
 story=set(sum(json.loads((O/'config/fixed_sets.json').read_text())['bird' if clade=='bird' else 'mammal'].values(),[]));selected=set(t.index[t.gene.isin(story)])
 z=f[f.sequence_name.isin(selected)].merge(t.reset_index(),on='sequence_name',how='left',suffixes=('_site','_promoter'));z['clade']=clade;export.append(z)
pd.concat(export,ignore_index=True).to_csv(O/'results/fixed_group_PWM_sites.tsv.gz',sep='\t',index=False)
check('common_activity_not_overwritten_by_full',set(pd.read_csv(O/'results/prior_TF_activity.tsv',sep='\t').panel)==set(pd.read_csv(O/'config/panel_inventory.tsv',sep='\t').panel))
check('full_activity_names_distinct',set(pd.read_csv(O/'results/full_prior_TF_activity.tsv',sep='\t').panel)==set(pd.read_csv(O/'config/full_inventory.tsv',sep='\t').panel))
# Recompute PC1 score from exported weights: verifies transform, gene/sample alignment and selected full scope.
for panel,module in [('ma2016','red'),('ma2016_full','brown')]:
 m=pd.read_csv(O/'inputs'/f'{panel}.log.tsv',sep='\t',index_col=0);a=pd.read_csv(O/'results/wgcna'/f'{panel}__raw'/'members.tsv',sep='\t');a=a[a.module.eq(module)];x=m.loc[a.feature].to_numpy();z=(x-a.mean_log.to_numpy()[:,None])/a.sd_log.to_numpy()[:,None];score=a.weight.to_numpy()@z;e=pd.read_csv(O/'results/wgcna'/f'{panel}__raw'/'eigengenes.tsv',sep='\t');ref=e[e.module.eq(module)].set_index('sample').loc[m.columns].eigengene.to_numpy();check(panel+'_PC1_export_reconstruction',np.allclose(score,ref,rtol=1e-9,atol=1e-9))
env={'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__}
for command,label in [([sys.executable,'-m','pip','freeze'],'pip_freeze'),(['/path/to/user-home/miniconda3/envs/meme/bin/fimo','--version'],'fimo'),(['git','-C',str(H),'rev-parse','HEAD'],'git_commit'),(['git','-C',str(H),'status','--porcelain'],'git_dirty')]:
 r=subprocess.run(command,capture_output=True,text=True);env[label]={'returncode':r.returncode,'output':r.stdout,'stderr':r.stderr}
(O/'audit/environment.json').write_text(json.dumps(env,indent=2));(O/'audit/final_checks.json').write_text(json.dumps(dict(checks=checks,all_passed=True,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2))
print('QA',len(checks),'passed',flush=True)
