from common import *
plan=json.loads((O/'config/analysis_plan.initial.json').read_text());sets=json.loads((O/'config/fixed_sets.json').read_text());inv=pd.read_csv(O/'config/panel_inventory.tsv',sep='\t')
groups=[];nodes=[];panelq=[]
for panel in inv.panel:
 p=load(panel);clade='bird' if panel.startswith('bird') else 'mammal';M=pd.read_csv(O/'inputs'/f'{panel}.linear.tsv.gz',sep='\t',index_col=0)
 # Family-level exact membership only when a family has one unambiguous symbol.
 for group,gg in sets[clade].items():
  f=p['feature'];ff=f.index[f.symbols.isin(gg)].intersection(p['X'].index);ss=set(f.loc[ff,'symbols']);complete=set(gg)<=ss
  if len(ff)<2:continue
  lx=M.loc[ff].sum(axis=0);summ=np.log2(lx+1).to_frame('v').join(p['meta'].set_index('sample').species).groupby('species',sort=False).v.mean().loc[p['sp']].to_numpy()
  z=p['X'].loc[ff];sc=z.sub(z.mean(axis=1),axis=0).div(z.std(axis=1),axis=0).mean(axis=0).to_numpy()
  for typ,y in [('log_sum_linear',summ),('mean_gene_Z',sc)]:
   groups.append(dict(panel=panel,group=group,score=typ,n_requested=len(gg),n_measured=len(ff),complete=complete,measured_symbols=';'.join(sorted(ss)),**fit(p['mass'],y,p['C'],loo=True)))
 for gene in list(dict.fromkeys(plan['candidate_TFs']+plan['non_PWM_nodes'])):
  if gene not in p['G'].index:
   nodes.append(dict(panel=panel,gene=gene,measured=False));continue
  mm=pd.read_csv(O/'results/wgcna'/f'{panel}__raw'/'members.tsv',sep='\t');ff=p['feature'].index[p['feature'].symbols.eq(gene)];module=';'.join(mm[mm.feature.isin(ff)].module);nodes.append(dict(panel=panel,gene=gene,measured=True,module=module,**fit(p['mass'],p['G'].loc[gene].to_numpy(),p['C'],loo=True)))
for name,rows in [('fixed_group_mass.tsv',groups),('growth_node_expression_mass.tsv',nodes)]:
 d=pd.DataFrame(rows);d['q_all_tests']=bh(d.p);d.to_csv(O/'results'/name,sep='\t',index=False)
d=pd.read_csv(O/'results/module_mass.tsv',sep='\t');d['q_within_panel_mode_sensitivity']=np.nan
for _,ix in d.groupby(['panel','mode']).groups.items():d.loc[ix,'q_within_panel_mode_sensitivity']=bh(d.loc[ix,'p'])
d.to_csv(O/'results/module_mass_q_sensitivity.tsv',sep='\t',index=False)
# An independent whitening check against base R lm is prepared without importing this fitter into R.
tests=[]
for panel,module in [('ma2016','red'),('bird_F12','orange'),('bird_E_label','white')]:
 p=load(panel);e=pd.read_csv(O/'results/wgcna'/f'{panel}__raw'/'eigengenes.tsv',sep='\t');e=e[e.module.eq(module)].merge(p['meta'][['sample','species']],on='sample').groupby('species').eigengene.mean().reindex(p['sp'])
 inp=O/'audit'/f'verify_{panel}.tsv';pd.DataFrame(dict(species=p['sp'],x=p['mass'],y=e.to_numpy())).to_csv(inp,sep='\t',index=False);pd.DataFrame(p['C']).to_csv(O/'audit'/f'verify_{panel}_C.tsv',sep='\t',index=False,header=False)
 tests.append(dict(panel=panel,module=module,**fit(p['mass'],e.to_numpy(),p['C'])))
pd.DataFrame(tests).to_csv(O/'audit/GLS_python.tsv',sep='\t',index=False)
(O/'audit/additional_execution.json').write_text(json.dumps(dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),node_scope='RNA expression only, not activation of protein signaling',fixed_groups='missing genes are not zero; complete flag required for full-group claim',within_panel_q='sensitivity only; global q across all networks remains prespecified primary'),indent=2))
