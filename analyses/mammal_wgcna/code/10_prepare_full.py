from common import *
V=R/'nonbird_regulatory_followup_2026-10-21_v1';A=R/'nonbird_common_state_2026-10-21_v1';plan=json.loads((O/'config/analysis_plan.initial.json').read_text());sets=json.loads((O/'config/fixed_sets.json').read_text())
anchors=set(sum(sets['mammal'].values(),[]))|set(plan['candidate_TFs'])|set(plan['non_PWM_nodes']);manifest=[]
for base in ['ma2016','gse105160','wagner']:
 panel=base+'_full';p=load(base);src=V/'results/ma2016_admitted_M.tsv.gz' if base=='ma2016' else A/'results'/f'{base}_admission_before_gene_filter_M.tsv.gz'
 M=pd.read_csv(src,sep='\t',index_col=0).loc[:,p['sp']];M=M.loc[np.isfinite(M).all(axis=1)&M.ge(0).all(axis=1)];M=M.loc[(M.ge(1).sum(axis=1)>=max(2,int(np.ceil(M.shape[1]*.5))))&(np.log2(M+1).std(axis=1)>1e-8)]
 M.to_csv(O/'inputs'/f'{panel}.linear.tsv.gz',sep='\t');np.log2(M+1).to_csv(O/'inputs'/f'{panel}.log.tsv',sep='\t');p['meta'].to_csv(O/'inputs'/f'{panel}.samples.tsv',sep='\t',index=False);pd.DataFrame(dict(feature=M.index,symbols=M.index,primary_symbol=M.index)).to_csv(O/'inputs'/f'{panel}.features.tsv',sep='\t',index=False);np.savez_compressed(O/'inputs'/f'{panel}.cov.npz',C=p['C'],species=np.array(p['sp'],dtype=str))
 manifest.append(dict(panel=panel,n_features=len(M),n_libraries=M.shape[1],n_species=M.shape[1],source_sha256=hashlib.sha256(src.read_bytes()).hexdigest()))
pd.DataFrame(manifest).to_csv(O/'config/full_inventory.tsv',sep='\t',index=False)
(O/'config/forced_anchors.json').write_text(json.dumps(sorted(anchors)))
(O/'config/full_scope_amendment.json').write_text(json.dumps(dict(reason='three-panel common universe removes 9 own5 genes due to Wagner annotation; full Ma/GSE own5 responses needed',network='top5000 MAD plus fixed own5 and a priori growth candidate anchors; selection traitblind but own5 historical selection remains',before_full_network_effect_tests=True,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2))
print(pd.DataFrame(manifest).to_string(index=False))
