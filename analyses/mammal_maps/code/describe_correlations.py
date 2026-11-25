"""Requested descriptive interpretation; retains phylogenetic results separately."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np,pandas as pd
from scipy import stats
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--results-root',type=Path,required=True);ap.add_argument('--genesets-root',type=Path,required=True);a=ap.parse_args();O=a.out
amendment=dict(user_request='Наша задача просто выделить млекопитающие связи; если что-то коррелирует с размером, это уже интересно.',primary_interpretation='Descriptive Pearson correlation of log10 linear normalized map sum with log10 body mass. Phylogeny is a separate diagnostic, not a gate.',timing='After initial full screen, directly requested by user; Pearson outcomes were already calculated and are not new discoveries from refitting.',multiple_testing='BH over all eligible maps within each panel and composition mode; replication shown as context, not an admission criterion.',minimum_measured_genes=5,human_review='not_completed',training_enabled=False)
(O/'config/descriptive_focus_amendment.json').write_text(json.dumps(amendment,ensure_ascii=False,indent=2))
T=pd.read_csv(O/'results/map_tests.tsv.gz',sep='\t');cat=pd.read_csv(O/'config/catalog.tsv.gz',sep='\t');hits=T[T.q_Pearson_panel.lt(.05)].copy().sort_values(['mode','panel','p_Pearson']);hits.to_csv(O/'results/Pearson_BH_hits.tsv',sep='\t',index=False)
ids=set(hits.set_id);score_rows=[];diag=[];contrib=[];err_r=err_p=0.;checks=0;pairrows=[]
def bh_manual(ps):
 p=np.asarray(ps);order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
for (mode,panel),z in T.groupby(['mode','panel']):
 design=np.load(O/'results'/f'{panel}_design.npz');sc=np.load(O/'results'/f'{mode}_{panel}_scores.npz');lookup={cat.iloc[k].set_id:i for i,k in enumerate(sc['catalog_index'])};x=design['mass_log10'];sp=list(design['species'])
 assert np.max(abs(bh_manual(z.p_Pearson)-z.q_Pearson_panel))<1e-10
 chosen=pd.concat([z[z.set_id.isin(ids)],z.nsmallest(10,'p_Pearson'),z.sample(25,random_state=2307)]).drop_duplicates('set_id')
 src=a.results_root/('nonbird_regulatory_followup_2026-10-21_v1/results/ma2016_admitted_M.tsv.gz' if panel=='ma2016' else f'nonbird_common_state_2026-10-21_v1/results/{panel}_admission_before_gene_filter_M.tsv.gz');M=pd.read_csv(src,sep='\t',index_col=0)[sp]
 for row in chosen.itertuples():
  y=sc['log10_sum'][lookup[row.set_id]];res=stats.pearsonr(x,y);err_r=max(err_r,abs(res.statistic-row.r_Pearson));err_p=max(err_p,abs(res.pvalue-row.p_Pearson));checks+=1
  if row.set_id not in ids:continue
  rr=[];pp=[];bb=[]
  for i,s in enumerate(sp):
   keep=np.arange(len(sp))!=i;fit=stats.linregress(x[keep],y[keep]);rr.append(fit.rvalue);pp.append(fit.pvalue);bb.append(fit.slope)
  gm=M.loc[row.genes.split(';')];share=gm.div(gm.sum(0),axis=1).mean(1).sort_values(ascending=False)
  diag.append(dict(mode=mode,panel=panel,set_id=row.set_id,n_genes=len(gm),loo_min_r=min(rr),loo_max_r=max(rr),loo_min_abs_r=min(map(abs,rr)),loo_sign_fraction=float(np.mean(np.sign(rr)==np.sign(row.r_Pearson))),loo_max_p=max(pp),top_gene=share.index[0],top_gene_mean_share=float(share.iloc[0]),top3_mean_share=float(share.iloc[:3].sum())))
  for i,s in enumerate(sp):score_rows.append(dict(mode=mode,panel=panel,set_id=row.set_id,species=s,log10_mass=x[i],log10_sum=y[i],linear_sum=10**y[i]))
  for gene in share.index:
   r,p=stats.pearsonr(x,np.log10(gm.loc[gene].to_numpy()+1));contrib.append(dict(mode=mode,panel=panel,set_id=row.set_id,gene=gene,mean_linear_share=float(share[gene]),individual_log1p_r=float(r),individual_log1p_p=float(p)))
pd.DataFrame(score_rows).to_csv(O/'results/Pearson_hit_species_scores.tsv',sep='\t',index=False);pd.DataFrame(diag).to_csv(O/'results/Pearson_hit_diagnostics.tsv',sep='\t',index=False);pd.DataFrame(contrib).to_csv(O/'results/Pearson_hit_genes.tsv.gz',sep='\t',index=False)
for mode,z in T.groupby('mode'):
 pv=z.pivot(index='set_id',columns='panel',values='p_Pearson');rv=z.pivot(index='set_id',columns='panel',values='r_Pearson')
 same=rv.gt(0).all(1)|rv.lt(0).all(1);complete=pv.notna().all(1);pmax=pv.max(1).where(same,1.);q=bh_manual(pmax[complete]);f=pd.DataFrame(dict(set_id=pv.index[complete],mode=mode,p_max_same_direction=pmax[complete],q_all3=q));pairrows.append(f)
pd.concat(pairrows).to_csv(O/'results/Pearson_all3_context.tsv.gz',sep='\t',index=False)
summary=dict(question='Какие суммы экспрессии полных измеренных частей KEGG/GO/Reactome коррелируют с массой млекопитающих?',primary_metric=amendment['primary_interpretation'],panels=[],common_composition_hits=hits[hits['mode'].eq('all3_common')][['panel','set_id','n_measured','r_Pearson','p_Pearson','q_Pearson_panel']].to_dict('records'),significant_main_terms=hits[hits['mode'].eq('panel_available')][['panel','set_id','n_measured','n_canonical','r_Pearson','beta_OLS','p_Pearson','q_Pearson_panel']].to_dict('records'),independent_checks=dict(Pearson_fits=checks,max_r_error=err_r,max_p_error=err_p,BH_every_row_verified=True),interpretation='18 whole measured-map associations in GSE105160; contextual differences across panels are visible and not used to dismiss those descriptive associations. Three additional Ma2016 associations refer to the smaller common composition, not the same complete sums.',human_review='not_completed',training_enabled=False)
for (mode,panel),z in T.groupby(['mode','panel']):summary['panels'].append(dict(mode=mode,panel=panel,n_species=int(z.n_species.iloc[0]),n_maps=len(z),nominal=int(z.p_Pearson.lt(.05).sum()),BH05=int(z.q_Pearson_panel.lt(.05).sum()),min_q=float(z.q_Pearson_panel.min())))
assert err_r<1e-9 and err_p<1e-9
(O/'results/descriptive_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));(O/'audit/Pearson_verification.json').write_text(json.dumps(summary['independent_checks'],indent=2))
# Preserve the source of derived GO TSV as well as the derived memberships.
src=a.genesets_root/'mammal_catalogs/msigdb.2026.1.Hs.C5.rds';(O/'audit/GO_source.json').write_text(json.dumps(dict(reference='mammal_catalogs/msigdb.2026.1.Hs.C5.rds',sha256=hashlib.sha256(src.read_bytes()).hexdigest(),bytes=src.stat().st_size,conversion='code/export_go.R'),indent=2))
print(json.dumps(dict(panels=summary['panels'],checks=summary['independent_checks']),ensure_ascii=False,indent=2))
