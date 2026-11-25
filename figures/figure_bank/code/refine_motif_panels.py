"""Re-render only the two motif tables exported by draw_scientific_panels.py.
Replace small dots with legible numeric cells; input values and tests unchanged.
"""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();O=a.out
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42})
SP=['Zebra','Quail','Chicken','Duck','Peacock','Turkey','Emu'];LABEL=['Zebra finch','Quail','Chicken','Duck','Peafowl','Turkey','Emu']
for name,qcol,tfs,family in [('matched21','q_all_matched_tests',['ESRRA','ESRRB','ESRRG','NR5A2'],'21 candidate models; BH across 1,617 matched tests'),('broad816','q_all_discovery_tests',['ESRRA','ESRRB','ESRRG','NR5A2','NR2F2','PPARA'],'816-model screen; BH across 11,424 tests')]:
 d=pd.read_csv(O/'data'/f'motif_{name}.tsv',sep='\t');fig,ax=plt.subplots(figsize=(13.8,3.5));fig.subplots_adjust(left=.12,right=.91,bottom=.22,top=.86)
 for i,sp in enumerate(SP):
  for j,tf in enumerate(tfs):
   z=d[(d.species==sp)&(d.tf==tf)].iloc[0];q=float(z[qcol]);level=-np.log10(max(q,1e-6));color=plt.cm.Purples(min(level/3,1));lum=.2126*color[0]+.7152*color[1]+.0722*color[2]
   ax.add_patch(Rectangle((j-.465,i-.44),.93,.88,facecolor=color,edgecolor='#34213E' if q<.05 else '#D0D0D0',lw=1.4 if q<.05 else .5))
   ax.text(j,i,f'{int(z.n_hit)}/{int(z.n_genes)}',ha='center',va='center',fontsize=11,color='white' if lum<.55 else '#202020')
 ax.set(xticks=range(len(tfs)),xticklabels=tfs,yticks=range(7),yticklabels=LABEL,xlim=(-.5,len(tfs)-.5),ylim=(6.5,-.5));ax.tick_params(length=0,labelsize=12)
 for spine in ax.spines.values():spine.set_visible(False)
 sm=plt.cm.ScalarMappable(norm=plt.Normalize(0,3),cmap='Purples');cb=fig.colorbar(sm,ax=ax,fraction=.025,pad=.025);cb.set_label('-log10 BH q',fontsize=11)
 fig.suptitle('OXPHOS16 promoter motif evidence',fontsize=17)
 fig.text(.12,.03,f'{family}. Counts: motif-positive genes / 16. Dark outline: q < 0.05.',fontsize=10,color='#444444')
 for ext in ['svg','png','pdf']:fig.savefig(O/'panels'/f'F4A_motif_evidence_{name}.{ext}',dpi=240,facecolor='white')
 plt.close(fig)
(O/'audit/motif_refinement.json').write_text(json.dumps({'change':'numeric cells replace dots to ensure readable counts','data_changed':False,'new_tests':False,'source':'data/motif_matched21.tsv and data/motif_broad816.tsv'},indent=2)+'\n')
