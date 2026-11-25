"""Render fixed, existing P-KLB measurements on the research host.

No discovery, gene selection by mass, new network fitting, or new motif scans.
The two expression versions and the candidate/broad motif families stay separate.
Only the plotted points, existing tests, verified sites and PWMs are exported.
"""
from pathlib import Path
import argparse, hashlib, json, platform, csv, re
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D

p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--later-results',type=Path,required=True);a=p.parse_args()
B=a.root;R=B/'results';O=a.out
for f in ['panels','data','audit']: (O/f).mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.8,'savefig.facecolor':'white'})
SP=['Zebra','Quail','Chicken','Duck','Peacock','Turkey','Emu']
LABEL=['Zebra finch','Quail','Chicken','Duck','Peafowl','Turkey','Emu']
MASS=np.array([11.93,96.6,912.5,1121.,4093.75,5811.,36200.]);X=np.log10(MASS)
IDS=['folate_12','carbon_7','amino_5','oxphos_16','epi_lipid_3']
NAMES=['Folate metabolism','Carbon metabolism','Amino-acid biosynthesis','Oxidative phosphorylation','Lipid-related set']
sources=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def register(p):
 area='legacy_project' if p.is_relative_to(B) else 'later_results';base=B if area=='legacy_project' else a.later_results
 rec={'area':area,'path':str(p.relative_to(base)),'sha256':sha(p),'bytes':p.stat().st_size}
 if rec not in sources:sources.append(rec)
 return p
def read(p,**kw):register(p);return pd.read_csv(p,sep='\t',**kw)
def export(fig,name):
 for ext in ['svg','png','pdf']:fig.savefig(O/'panels'/f'{name}.{ext}',dpi=240)
 plt.close(fig)
def save(name,d): (O/'data'/name).write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')

members=read(a.later_results/'sizes6_three_gaps_2026-10-19_v1/01_selection/v2/locked_memberships.tsv')
sets={k:list(g.gene) for k,g in members.groupby('set_id',sort=False)}
def matrix(path):
 d=read(path,index_col=0)[SP];d.index=d.index.str.upper();assert not d.index.duplicated().any();return d
legacy=R/'oxphos_species_thread'
matrices={
 'submitted':{'fib':matrix(legacy/'honest_gene_scan/M_fib_nearest2kb.tsv'),'epi':matrix(legacy/'story_epi_traits/M_epi_nearest2kb.tsv')},
 'corrected':{'fib':matrix(R/'body_mass_plan_2026-09-29/M_fib_verified40.tsv'),'epi':matrix(legacy/'epi_own_promoters/M_epi_own_nearest2kb.tsv')}
}
oldtests=read(legacy/'nine_set_slopes/nine_set_slopes.tsv')
allstats=[];points=[];scores={}
for version,mats in matrices.items():
 for tissue,m in mats.items():
  rows=[]
  for sid,genes in sets.items():
   assert set(genes)<=set(m.index),(version,tissue,sid)
   values=m.loc[genes].to_numpy(float);assert np.isfinite(values).all() and (values>=0).all()
   sums=values.sum(0);assert (sums>0).all();y=np.log10(sums);r,pval=stats.pearsonr(X,y)
   slope,intercept=np.polyfit(X,y,1);scores[version,tissue,sid]=y
   rows.append(dict(version=version,tissue=tissue,set=sid,n_genes=len(genes),n_species=7,r=float(r),p=float(pval),slope=float(slope),intercept=float(intercept)))
   for i,s in enumerate(SP):points.append(dict(version=version,tissue=tissue,set=sid,species=s,mass_g=float(MASS[i]),sum_linear=float(sums[i]),log10_mass=float(X[i]),log10_sum=float(y[i])))
  ps=np.array([v['p'] for v in rows]);order=np.argsort(ps);q=np.empty(len(ps));q[order]=np.minimum(1,np.minimum.accumulate((ps[order]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1])
  for z,qq in zip(rows,q):z['q_BH_9_within_state']=float(qq)
  if version=='submitted':
   for z in rows:
    old=oldtests[(oldtests['set']==z['set'])&(oldtests.tissue==tissue)].iloc[0]
    assert abs(z['r']-old.r)<1e-11,(z,old.to_dict())
    assert abs(z['q_BH_9_within_state']-old.bh)<1e-11
  allstats.extend(rows)
pd.DataFrame(points).to_csv(O/'data/bird_points.tsv',sep='\t',index=False)
pd.DataFrame(allstats).to_csv(O/'data/bird_tests.tsv',sep='\t',index=False)
members.to_csv(O/'data/frozen_gene_sets.tsv',sep='\t',index=False)

def annotate(ax,x,y,labels,fontsize=10):
 """Greedy screen-space label placement, preserving all observations."""
 fig=ax.figure;fig.canvas.draw();renderer=fig.canvas.get_renderer();used=[];bounds=ax.get_window_extent();avoid=None
 for xx,yy,label in zip(x,y,labels):
  best=None
  for dx,dy,ha in [(6,7,'left'),(-6,7,'right'),(6,-16,'left'),(-6,-16,'right'),(6,21,'left'),(-6,21,'right'),(6,-29,'left'),(-6,-29,'right')]:
   t=ax.annotate(label,(xx,yy),xytext=(dx,dy),textcoords='offset points',ha=ha,va='bottom',fontsize=fontsize,color='#222222')
   bb=t.get_window_extent(renderer).expanded(1.04,1.12)
   overlap=sum(bb.overlaps(b) for b in used)
   outside=max(bounds.x0-bb.x0,0)+max(bb.x1-bounds.x1,0)+max(bounds.y0-bb.y0,0)+max(bb.y1-bounds.y1,0)
   rank=overlap*1e5+outside*100+abs(dy)
   t.remove()
   if best is None or rank<best[0]:best=(rank,dx,dy,ha,bb)
  _,dx,dy,ha,bb=best;used.append(bb)
  ax.annotate(label,(xx,yy),xytext=(dx,dy),textcoords='offset points',ha=ha,va='bottom',fontsize=fontsize,color='#222222',arrowprops={'arrowstyle':'-','lw':.45,'color':'#777777','shrinkA':1,'shrinkB':3})

def bird_axis(ax,version,tissue,sid,ylim,title=None,individual=False):
 y=scores[version,tissue,sid];z=next(z for z in allstats if (z['version'],z['tissue'],z['set'])==(version,tissue,sid))
 ax.plot([.65,4.75],z['intercept']+z['slope']*np.array([.65,4.75]),color='#2864A0',lw=1.5,zorder=1)
 ax.scatter(X,y,s=39,color='#161616',edgecolors='white',linewidth=.5,zorder=3)
 ax.set(xlim=(.65,4.8),ylim=ylim);ax.set_xticks([1,2,3,4]);ax.grid(alpha=.16,lw=.6)
 annotate(ax,X,y,LABEL,fontsize=9.3 if not individual else 10)
 ax.tick_params(labelsize=10 if not individual else 11)
 q=z['q_BH_9_within_state'];qt=f'{q:.3g}'
 ax.set_title((title+'\n' if title else '')+f"n = {z['n_genes']} genes   r = {z['r']:+.3f}   BH q = {qt}",fontsize=11 if not individual else 13,pad=12,loc='left')
 if individual:ax.set_xlabel('log10 adult body mass (g)',fontsize=12);ax.set_ylabel('log10 summed CAGE signal',fontsize=12)

for version in matrices:
 fig,axes=plt.subplots(2,5,figsize=(29.2,8.55));fig.subplots_adjust(left=.045,right=.993,bottom=.12,top=.86,wspace=.32,hspace=.48)
 for j,(sid,name) in enumerate(zip(IDS,NAMES)):
  yy=np.r_[scores[version,'epi',sid],scores[version,'fib',sid]];pad=(yy.max()-yy.min())*.28+.025;lim=(float(yy.min()-pad),float(yy.max()+pad))
  for i,tissue in enumerate(['epi','fib']):
   bird_axis(axes[i,j],version,tissue,sid,lim,name if i==0 else None)
   if j==0:axes[i,j].set_ylabel(('Epiblast' if i==0 else 'Fibroblasts')+'\nlog10 summed CAGE signal',fontsize=15)
   f,ax=plt.subplots(figsize=(5.85,3.92));f.subplots_adjust(left=.18,right=.97,bottom=.20,top=.80)
   bird_axis(ax,version,tissue,sid,lim,name,True);export(f,f'F2_{version}_{tissue}_{sid}')
 fig.supxlabel('log10 adult body mass (g)',y=.045,fontsize=16)
 label='Submitted-analysis source matrices' if version=='submitted' else 'Corrected fibroblast40 + own-processed epiblast matrices'
 fig.suptitle(label,fontsize=19,y=.99)
 fig.text(.045,.005,'Seven species. Fixed historical gene sets; BH across all nine sets within each state. Estimates are descriptive after historical selection.',fontsize=12,color='#444444')
 export(fig,f'F2_bird_grid_{version}')

# Existing sensitivity example: what contributes to the folate sum?
fig,axs=plt.subplots(1,2,figsize=(11.6,4.4));fig.subplots_adjust(left=.09,right=.98,top=.82,bottom=.23,wspace=.3)
fol=sets['folate_12'];summaries=[]
for tissue,color,marker,label in [('epi','#3268A6','o','Epiblast'),('fib','#A54678','s','Fibroblasts')]:
 m=matrices['corrected'][tissue];s=m.loc[fol].sum(0).to_numpy(float);share=m.loc['AHCY'].to_numpy(float)/s
 axs[0].plot(X,share*100,marker=marker,color=color,label=label,lw=1.2)
 without=np.log10(m.loc[[g for g in fol if g!='AHCY']].sum(0).to_numpy(float));r,pv=stats.pearsonr(X,without)
 axs[1].scatter(X,without,marker=marker,color=color,label=f'{label}: r = {r:+.3f}',s=40)
 coef=np.polyfit(X,without,1);axs[1].plot(X,np.polyval(coef,X),color=color,lw=1.1)
 for i,sp in enumerate(SP):summaries.append(dict(tissue=tissue,species=sp,AHCY_fraction=float(share[i]),log10_folate_without_AHCY=float(without[i])))
axs[0].set(ylabel='AHCY share of the folate-set signal (%)',xlabel='log10 adult body mass (g)',title='Contribution of one abundant component')
axs[1].set(ylabel='log10 CAGE sum of the other 11 genes',xlabel='log10 adult body mass (g)',title='The fixed set with AHCY omitted')
for ax in axs:ax.grid(alpha=.16);ax.legend(frameon=False,fontsize=10)
fig.suptitle('Folate12: component-level interpretation',fontsize=17)
fig.text(.09,.025,'Sensitivity analysis on the same seven species; no genes selected by their mass correlation in this display.',fontsize=10,color='#444444')
export(fig,'F3_folate_component_example');pd.DataFrame(summaries).to_csv(O/'data/folate_component_points.tsv',sep='\t',index=False)

# Existing candidate-matched and broad motif tests remain separate families.
reg=R/'regulatory_mechanism_2026-09-29'
motif_specs=[('matched21','promoter_motif_matched_tests.tsv','q_all_matched_tests',['ESRRA','ESRRB','ESRRG','NR5A2']),('broad816','all_tf_promoter_refined.tsv','q_all_discovery_tests',['ESRRA','ESRRB','ESRRG','NR5A2','NR2F2','PPARA'])]
for name,file,qcol,tfs in motif_specs:
 d=read(reg/file);d=d[(d['set']=='oxphos_16')&d.tf.isin(tfs)].copy();d.to_csv(O/'data'/f'motif_{name}.tsv',sep='\t',index=False)
 fig,ax=plt.subplots(figsize=(13.8,3.8));fig.subplots_adjust(left=.12,right=.90,bottom=.22,top=.87)
 for i,sp in enumerate(SP):
  for j,tf in enumerate(tfs):
   z=d[(d.species==sp)&(d.tf==tf)];assert len(z)==1;(z,)=z.to_dict('records');q=float(z[qcol]);hits=int(z.get('n_hit',z.get('pp',0)));n=int(z.get('n_genes',z.get('n_pos',16)))
   ax.scatter(j,i,s=80+260*hits/max(n,1),c=[-np.log10(max(q,1e-6))],vmin=0,vmax=3,cmap='Purples',edgecolors='#3A2451' if q<.05 else '#B8BCC1',linewidths=1.2 if q<.05 else .7)
   ax.text(j,i,f'{hits}/{n}',ha='center',va='center',fontsize=9,color='white' if q<.015 else '#252525')
 ax.set(xticks=range(len(tfs)),xticklabels=tfs,yticks=range(7),yticklabels=LABEL,xlim=(-.6,len(tfs)-.4),ylim=(6.6,-.6));ax.tick_params(length=0,labelsize=12)
 for spine in ax.spines.values():spine.set_visible(False)
 sm=plt.cm.ScalarMappable(norm=plt.Normalize(0,3),cmap='Purples');cb=fig.colorbar(sm,ax=ax,fraction=.035,pad=.035);cb.set_label('-log10 BH q',fontsize=11)
 fig.suptitle('OXPHOS16 promoter motif evidence',fontsize=17)
 family='21 candidate models; BH across 1,617 matched tests' if name=='matched21' else '816-model screen; BH across 11,424 tests'
 fig.text(.12,.025,f'{family}. Numbers: motif-positive genes / tested genes. Dark outline: q < 0.05.',fontsize=10,color='#444444')
 export(fig,f'F4A_motif_evidence_{name}')

# Actual JASPAR probability matrices, not illustrative consensus sequences.
meme=register(legacy/'module_consistency/fimo_chicken/JASPAR2024_vertebrates.meme');lines=meme.read_text().splitlines();pwms={}
for i,line in enumerate(lines):
 if line.startswith('MOTIF '):
  bits=line.split();mid,tf=bits[1:3];tf=tf.upper()
  if tf not in ['ESRRA','ESRRB','ESRRG','NR5A2']:continue
  j=i+1
  while not lines[j].startswith('letter-probability'):j+=1
  n=int(re.search(r'w=\s*(\d+)',lines[j]).group(1));arr=np.array([[float(v) for v in lines[k].split()] for k in range(j+1,j+1+n)])
  pwms[tf]=(mid,arr)
assert len(pwms)==4
letters=['A','C','G','T'];colors={'A':'#238E55','C':'#2864A0','G':'#DAA124','T':'#C44149'};font=FontProperties(family='DejaVu Sans',weight='bold')
fig,axs=plt.subplots(1,4,figsize=(13.8,2.3));fig.subplots_adjust(left=.035,right=.995,bottom=.25,top=.72,wspace=.24)
logo_export={}
for ax,tf in zip(axs,['ESRRA','ESRRB','ESRRG','NR5A2']):
 mid,prob=pwms[tf];prob=prob/prob.sum(1,keepdims=True);info=2+np.sum(prob*np.log2(np.maximum(prob,1e-300)),axis=1);hh=prob*info[:,None]
 for pos,row in enumerate(hh):
  y=0
  for idx in np.argsort(row):
   h=float(row[idx]);letter=letters[idx]
   if h<.005:continue
   tp=TextPath((0,0),letter,size=1,prop=font);bb=tp.get_extents();trans=Affine2D().translate(-bb.x0,-bb.y0).scale(.86/bb.width,h/bb.height).translate(pos+.07,y)
   ax.add_patch(PathPatch(tp,transform=trans+ax.transData,color=colors[letter],lw=0));y+=h
 ax.set(xlim=(0,len(prob)),ylim=(0,2.05),xticks=[.5,len(prob)/2,len(prob)-.5],xticklabels=['1',str(round(len(prob)/2)),str(len(prob))],yticks=[0,1,2]);ax.tick_params(labelsize=8);ax.set_title(f'{tf}  {mid}',fontsize=12);ax.set_xlabel('Motif position',fontsize=9)
 logo_export[tf]={'motif_id':mid,'probability_matrix':prob.tolist()}
axs[0].set_ylabel('bits',fontsize=10);export(fig,'F4B_actual_NR_motif_logos');save('motif_probability_matrices.json',logo_export)

# Sequence-coordinate-verified examples; these are candidate windows, not proven TSSs.
win=read(reg/'promoter_candidate_windows.tsv');sites=read(reg/'sequence_verified_candidate_sites.tsv')
sites=sites[sites.tf.isin(['ESRRA','ESRRB','ESRRG','NR5A2'])&sites.sequence_matches_fimo.astype(str).str.lower().eq('true')]
coverage=sites.groupby('gene').species.nunique().sort_values(ascending=False)
examples=[g for g in ['ATP5F1E','NDUFB8','COX6C','NDUFB5'] if g in coverage.index][:2];assert len(examples)==2
selected=[]
for eno,gene in enumerate(examples,1):
 w=win[win.gene==gene];ss=sites[sites.gene==gene];species=[sp for sp in SP if sp in set(w.species)]
 fig,ax=plt.subplots(figsize=(13.8,3.4));fig.subplots_adjust(left=.10,right=.73,bottom=.24,top=.80)
 for i,sp in enumerate(species):
  wr=w[w.species==sp].iloc[0];peak=int(str(wr.peak_id).rsplit('_',1)[1]);sgn=1 if wr.gene_strand=='+' else -1
  lo=(int(wr.start_0based)-peak)*sgn;hi=(int(wr.end_exclusive)-peak)*sgn;lo,hi=sorted([lo,hi]);ax.hlines(i,lo,hi,color='#89929A',lw=1.1);ax.scatter([0],[i],marker='|',s=160,color='#222222',zorder=4)
  z=ss[(ss.species==sp)&(ss.peak_id==wr.peak_id)];intervals=[]
  for (start,end),rows in z.groupby(['start_0based','end_exclusive']):
   assert int(start)>=int(wr.start_0based) and int(end)<=int(wr.end_exclusive),(gene,sp,start,end,wr.to_dict())
   if sgn==1:left,right=int(start)-peak,int(end)-peak
   else:left,right=peak-int(end),peak-int(start)
   ax.add_patch(Rectangle((left,i-.16),right-left,.32,facecolor='#964B9A',edgecolor='white',lw=.4,zorder=3))
   intervals.append((float(rows.iloc[0].get('score',0)),rows))
   selected.extend({k:None if pd.isna(v) else v for k,v in zrow.items() if k!='genome'} for zrow in rows.to_dict('records'))
  seq=str(z.iloc[0].sequence) if len(z) else 'No verified hit for these models'
  ax.text(1.025,i,seq,transform=ax.get_yaxis_transform(),ha='left',va='center',fontsize=9.5,family='monospace',color='#53345C')
 allrel=[(float(v)-int(str(wr.peak_id).rsplit('_',1)[1]))*(1 if wr.gene_strand=='+' else -1) for _,wr in w.iterrows() for v in [wr.start_0based,wr.end_exclusive]]
 ax.set(xlim=(min(allrel)-25,max(allrel)+25),ylim=(len(species)-.5,-.5),yticks=range(len(species)),yticklabels=[LABEL[SP.index(sp)] for sp in species],xticks=[-200,-100,0,100,200],xlabel='Position relative to the CAGE-assigned reference peak (bp)')
 ax.set_title(f'{gene}: verified nuclear-receptor-like motif occurrences',fontsize=16,loc='left');ax.grid(axis='x',alpha=.15);ax.tick_params(axis='y',length=0,labelsize=10);ax.text(1.025,1.065,'Example matched sequence',transform=ax.transAxes,fontsize=10)
 fig.text(.10,.035,'Candidate sequence windows (unequal lengths). Rows are not a nucleotide alignment; motif similarity does not establish orthologous-site conservation or TF binding.',fontsize=9,color='#444444')
 export(fig,f'F4_promoter_example_{eno}_{gene}')
pd.DataFrame(selected).drop_duplicates().to_csv(O/'data/verified_NR_example_sites.tsv',sep='\t',index=False)
win[win.gene.isin(examples)].drop(columns=['genome']).to_csv(O/'data/verified_NR_example_windows.tsv',sep='\t',index=False)

save('source_manifest.json',{'host':'host1','inputs':sources,'method':'render fixed measurements and existing tests; no new gene search, FIMO or network fitting'})
checks={'submitted_correlations_and_BH_match_archived_nine_set_table':True,'n_exported_scatter_points':len(points),'all_gene_sets_complete':True,'frozen_sets':len(sets),'motif_examples':examples,'new_association_estimates':'Recomputed only for fixed existing sets and the already-described AHCY sensitivity; no new discovery.','version_boundary':'submitted matrices and corrected40/own-epiblast are separate outputs','limits':['Historical gene-set selection is not corrected by the displayed nine-test BH.','Not a phylogenetically adjusted or independent replication.','Verified example sequences do not establish active TSS, TF occupancy or orthologous motif-site conservation.']}
save('verification.json',checks)
(O/'audit/environment.json').write_text(json.dumps({'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':__import__('scipy').__version__,'matplotlib':matplotlib.__version__},indent=2)+'\n')
print(json.dumps(checks,indent=2),flush=True)
