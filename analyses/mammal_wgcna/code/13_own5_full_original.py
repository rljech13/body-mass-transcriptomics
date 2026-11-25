from common import *
sets=json.loads((O/'config/fixed_sets.json').read_text())['mammal'];V=R/'nonbird_regulatory_followup_2026-10-21_v1';A=R/'nonbird_common_state_2026-10-21_v1'
massrows=[];linkrows=[];modrows=[]
for base in ['ma2016','gse105160','wagner']:
 p=load(base+'_full');src=V/'results/ma2016_admitted_M.tsv.gz' if base=='ma2016' else A/'results'/f'{base}_admission_before_gene_filter_M.tsv.gz';M=pd.read_csv(src,sep='\t',index_col=0).loc[:,p['sp']];Z=np.log2(M+1);Z=Z.sub(Z.mean(axis=1),axis=0).div(Z.std(axis=1),axis=0)
 activity=pd.read_csv(O/'results/full_prior_TF_activity.tsv',sep='\t');activity=activity[activity.panel.eq(base+'_full')&activity.score.eq('background_centered_signed')].pivot(index='TF',columns='species',values='activity').reindex(columns=p['sp'])
 for group,gg in sets.items():
  seen=sorted(set(gg)&set(M.index));complete=set(gg)<=set(M.index)
  if len(seen)<2:continue
  ysum=np.log2(1+M.loc[seen].sum(axis=0)).to_numpy();yz=Z.loc[seen].mean(axis=0,skipna=False).to_numpy()
  for typ,y in [('log_sum_linear',ysum),('mean_gene_Z',yz)]:
   massrows.append(dict(panel=base,group=group,score=typ,n_requested=len(gg),n_measured=len(seen),complete=complete,missing=';'.join(sorted(set(gg)-set(seen))),**fit(p['mass'],y,p['C'],loo=True)))
  for tf,row in activity.iterrows():
   for typ,y in [('log_sum_linear',ysum),('mean_gene_Z',yz)]:
    linkrows.append(dict(panel=base,TF=tf,group=group,group_score=typ,n_requested=len(gg),n_measured=len(seen),complete=complete,**fit(row.to_numpy(),y,p['C'],extra=np.column_stack([p['mass'],p['background']]),loo=True)))
  for mode in ['raw','sample_centered']:
   mm=pd.read_csv(O/'results/wgcna'/f'{base}_full__{mode}'/'members.tsv',sep='\t');story=set(sum(sets.values(),[]))
   for module,a in mm[mm.module.ne('grey')].groupby('module'):
    a=a[~a.feature.isin(story)]
    if len(a)<20:continue
    x=p['X'].loc[a.feature]
    if mode=='sample_centered':x=x-p['X'].mean(axis=0)
    zz=x.sub(x.mean(axis=1),axis=0).div(x.std(axis=1),axis=0);sc=a.weight.to_numpy()@zz.to_numpy()
    modrows.append(dict(panel=base,mode=mode,module=module,group=group,n_group=len(seen),complete=complete,n_module_no_story=len(a),**fit(sc,yz,p['C'],extra=np.column_stack([p['mass'],p['background']]),loo=True)))
for name,rows in [('own5_full_original_mass.tsv',massrows),('own5_full_original_TF_links.tsv',linkrows),('own5_full_original_module_links.tsv',modrows)]:
 d=pd.DataFrame(rows);d['q_all_tests']=bh(d.p);d.to_csv(O/'results'/name,sep='\t',index=False);print(name,len(d),'q05',int(d.q_all_tests.lt(.05).sum()),flush=True)
(O/'audit/own5_full_original_execution.json').write_text(json.dumps(dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Own5 scores from original full species-admitted matrices before common-intersection or expression filtering. Complete flags retained; absent Wagner genes not set to zero.',regulator_links='heldout full-panel footprints; PGLS given mass AND heldout global expression background; no story-gene overlap in footprint or module score'),indent=2))
