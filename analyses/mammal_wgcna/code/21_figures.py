from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(2,3,figsize=(14,8.5),constrained_layout=True)
mm=pd.read_csv(O/'results/wgcna/ma2016__raw/members.tsv',sep='\t');mm=mm[mm.module.eq('red')]
for ax,panel in zip(axs[0],['ma2016','gse105160','wagner']):
 p=load(panel)
 if panel=='ma2016':
  e=pd.read_csv(O/'results/wgcna/ma2016__raw/eigengenes.tsv',sep='\t');y=e[e.module.eq('red')].set_index('sample').loc[p['sp']].eigengene.to_numpy()
 else:
  a=mm[mm.feature.isin(p['X'].index)];x=p['X'].loc[a.feature];z=x.sub(x.mean(axis=1),axis=0).div(x.std(axis=1),axis=0);y=a.weight.to_numpy()@z.to_numpy()
 y=(y-y.mean())/y.std(ddof=1);x=p['mass'];f=fit(x,y,p['C']);X=np.column_stack([np.ones(len(x)),x]);L=np.linalg.cholesky(p['C']);coef=np.linalg.lstsq(linalg.solve_triangular(L,X,lower=True),linalg.solve_triangular(L,y,lower=True),rcond=None)[0]
 ax.scatter(x,y,s=36,color='#315c8e');xx=np.array([x.min(),x.max()]);ax.plot(xx,coef[0]+coef[1]*xx,color='#a84935');ax.set(title=f'{panel}: frozen Ma cell-cycle module',xlabel='log10 adult species mass (g)',ylabel='Module score (within-panel SD)')
 ax.text(.02,.97,f'n = {len(x)} species\nPGLS p = {f["p"]:.3g}',va='top',transform=ax.transAxes,bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
ax=axs[1,0];d=pd.read_csv(O/'results/wgcna_LOO_membership.tsv',sep='\t');u=pd.read_csv(O/'results/wgcna_full_LOO_membership.tsv',sep='\t');d=pd.concat([d,u]);g=d.groupby(['panel','module']).Jaccard.median().reset_index();pan=['ma2016','ma2016_full','bird_HG','bird_F12'];ax.boxplot([g[g.panel.eq(p)].Jaccard for p in pan],labels=['Ma\ncommon','Ma\nfull','Bird\nHG','Bird\nF12']);ax.set(title='Module membership under species deletion',ylabel='Median best Jaccard over omitted species',ylim=(0,1))
ax=axs[1,1];j=pd.read_csv(O/'results/own5_module_gene_jackknife.tsv',sep='\t');j=j[j['mode'].eq('raw')&j.module.eq('brown')&j.check.eq('drop_one_gene')];full=pd.read_csv(O/'results/own5_full_original_module_links.tsv',sep='\t');pv=full[full.panel.eq('gse105160')&full['mode'].eq('raw')&full.module.eq('brown')&full.group.eq('mammal:glycogen, down')].p.iloc[0]
vals=[pv]+[j[j.gene.eq(g)].p.iloc[0] for g in ['GRB10','GSK3A','SELENOS']];ax.bar(range(4),-np.log10(vals),color=['#315c8e','#7098ba','#bd5d47','#7098ba']);ax.axhline(-np.log10(.05),color='grey',ls='--');ax.set_xticks(range(4),['All 3','−GRB10','−GSK3A','−SELENOS'],rotation=25);ax.set(title='GSE: module link after dropping one gene',ylabel='−log10 nominal PGLS p');ax.set_ylim(0,6);ax.text(.02,.97,'Given mass + expression background\nAll story genes removed from module',va='top',transform=ax.transAxes,fontsize=8)
ax=axs[1,2];ax.axis('off')
text=("What the data support\n\nE2F-associated cell-cycle structure\nis recognizable across panels.\n\nIts mass association does not transfer.\n\nThe glycogen-trio module link\nneeds GSK3A.\n\nPWM = sequence compatibility;\nregulon edge ≠ demonstrated binding.")
ax.text(.04,.95,text,va='top',fontsize=11,linespacing=1.5)
fig.suptitle('Growth-regulator search: program structure, mass association and stability',fontsize=15)
for ext in ['png','pdf','svg']:fig.savefig(O/'figures'/('regulator_audit.'+ext),dpi=180)
plt.close(fig);print('FIGURE DONE')
