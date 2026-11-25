from common import *
hallpath=H/'data/genesets/mammal_catalogs/hallmark.tsv';hall=pd.read_csv(hallpath,sep='\t');out=[]
for invfile in ['panel_inventory.tsv','full_inventory.tsv']:
 inv=pd.read_csv(O/'config'/invfile,sep='\t')
 for panel in inv.panel:
  p=load(panel)
  for mode in ['raw','sample_centered']:
   mm=pd.read_csv(O/'results/wgcna'/f'{panel}__{mode}'/'members.tsv',sep='\t');feat=p['feature'];universe=set(feat.loc[feat.index.intersection(mm.feature)].symbols.dropna().loc[lambda s:~s.str.contains(';')])
   for module,a in mm[mm.module.ne('grey')].groupby('module'):
    genes=set(feat.reindex(a.feature).symbols.dropna().loc[lambda s:~s.str.contains(';')])
    for _,h in hall.iterrows():
     target=set(str(h.genes).split(';'))&universe
     if len(target)<10:continue
     k=len(genes&target);odds,pv=stats.fisher_exact([[k,len(genes-target)],[len(target-genes),len(universe-genes-target)]],alternative='greater')
     out.append(dict(panel=panel,mode=mode,module=module,pathway=h.gs_name,n_targets=len(target),overlap=k,odds=odds,p=pv,overlap_genes=';'.join(sorted(genes&target))))
d=pd.DataFrame(out);d['q_all_tests']=bh(d.p);d.to_csv(O/'results/module_Hallmark_enrichment.tsv.gz',sep='\t',index=False)
links=pd.read_csv(O/'results/own5_full_original_module_links.tsv',sep='\t');links=links[links.q_all_tests.lt(.05)&links.complete]
sets=json.loads((O/'config/fixed_sets.json').read_text())['mammal'];rr=[]
for a in links.itertuples():
 p=load(a.panel+'_full');src=R/'nonbird_regulatory_followup_2026-10-21_v1/results/ma2016_admitted_M.tsv.gz' if a.panel=='ma2016' else R/'nonbird_common_state_2026-10-21_v1/results'/f'{a.panel}_admission_before_gene_filter_M.tsv.gz';M=pd.read_csv(src,sep='\t',index_col=0).loc[:,p['sp']];G=np.log2(M+1);G=G.sub(G.mean(axis=1),axis=0).div(G.std(axis=1),axis=0);mm=pd.read_csv(O/'results/wgcna'/f'{a.panel}_full__{a.mode}'/'members.tsv',sep='\t');mm=mm[mm.module.eq(a.module)&~mm.feature.isin(set(sum(sets.values(),[])))]
 x=p['X'].loc[mm.feature]
 if a.mode=='sample_centered':x=x-p['X'].mean(axis=0)
 x=x.sub(x.mean(axis=1),axis=0).div(x.std(axis=1),axis=0);sc=mm.weight.to_numpy()@x.to_numpy();gg=sets[a.group]
 for g in gg:
  for typ,target in [('one_gene',[g]),('drop_one_gene',[v for v in gg if v!=g])]:
   y=G.loc[target].mean(axis=0).to_numpy();rr.append(dict(panel=a.panel,mode=a.mode,module=a.module,group=a.group,gene=g,check=typ,targets=';'.join(target),**fit(sc,y,p['C'],extra=np.column_stack([p['mass'],p['background']]),loo=True)))
d=pd.DataFrame(rr);d['q_followup']=bh(d.p);d.to_csv(O/'results/own5_module_gene_jackknife.tsv',sep='\t',index=False)
(O/'audit/Hallmark_and_gene_jackknife_execution.json').write_text(json.dumps(dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),hallmark_sha256=hashlib.sha256(hallpath.read_bytes()).hexdigest(),jackknife_scope='follow-up of three significant complete-group module links; all requested gene deletions exported; not an independent replication'),indent=2));print('DONE',len(out),len(rr),flush=True)
