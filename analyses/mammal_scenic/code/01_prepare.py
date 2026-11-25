from pathlib import Path
import os,json,hashlib,sys,re,gzip,math,importlib.metadata as md
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];P=R.parent;H=Path('/path/to/user-home/chicken/kleiber_birds_gemcat');inputs=[]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<22),b''):h.update(b)
 return h.hexdigest()
def use(p):
 p=Path(p);inputs.append(dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size));return p
for d in ['inputs','config','reference','results','audit','logs','figures']:(R/d).mkdir(exist_ok=True)
plan=dict(run='ma2016_pyscenic_2026-10-23_v1',authorization='User: Ну так в чем проблема это сделать?',experiment='Ma2016 primary mammalian fibroblasts, same source as frozen WGCNA',technical_pairs=[['RedSquirrel.1A','RedSquirrel.1B'],['DeerMouse.2A','DeerMouse.2B']],normalization='sum raw counts within explicit technical pairs; CPM from all source gene counts; log2(CPM+1)',filter='CPM >=1 in >=20% of26 collapsed sample-label groups, nonconstant; no body mass or WGCNA selection',orthology='unique one-to-one mouse-to-human mappings, unique human symbols',seeds=[20261023,20261024,20261025],primary_seed=20261023,workers=24,regressor='Official arboreto GRNBoost2 defaults',cisTarget='Official custom database builder, cluster-buster, JASPAR2024 one-per-symbol816 motifs; standard downloads and official mirror refused connection',reference_windows={'500up100down':[-500,100],'10kb_each_side':[-10000,10000]},reference='Human GRCh38.99 longest exonic transcript on primary chromosomes; protein-coding genes plus expression-input symbols; unique gene symbols',pruning=dict(rank_threshold=1500,auc_threshold=.05,nes_threshold=3.,annotation='direct JASPAR TF label, heterodimer components expanded, HGNC unambiguous aliases normalized'),overlap='all inferred primary regulons x4 fixed WGCNA modules; measured/ranked WGCNA-network gene background; one-sided hypergeometric, BH across all tested module-regulon pairs',activity='AUCell continuous per donor-labelled group; species means for mass correlation; BH across primary regulons; exploratory same-study analysis',stability='two additional predeclared full GRN/prune runs, target Jaccard and TF recovery; not independent data',wgcna_rebuilt=False,slides_edited=False)
(R/'config/plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2)+'\n')
raw=pd.read_excel(use(H/'data/file_panels/ma2016/elife19130_supp1.xlsx'),sheet_name='A. raw counts',index_col=0)
meta0=pd.read_excel(use(H/'data/file_panels/ma2016/elife19130_fig1.xlsx'),sheet_name='B. RNAseq & Mapping to ortholog')
assert raw.index.is_unique and raw.columns.is_unique and raw.shape==(16816,28)
map0=pd.read_csv(use(H/'data/ensembl_orthologs/mouse_to_human.tsv'),sep='\t').dropna();map0=map0[map0['Human homology type'].eq('ortholog_one2one')]
map0=map0[['Gene stable ID','Human gene name']].drop_duplicates();map0=map0[~map0['Gene stable ID'].duplicated(keep=False)&~map0['Human gene name'].duplicated(keep=False)]
map0['Human gene name']=map0['Human gene name'].str.upper()
# Normalize unambiguous historic symbols to HGNC accepted names once, before matching.
hgnc=pd.read_csv(use(P/'sizes5_execution_2026-10-16_v1/sources/hgnc_complete_set.txt'),sep='\t',dtype=str).fillna('');official=set(hgnc.symbol.str.upper());aliases={}
for row in hgnc.itertuples():
 for col in ['alias_symbol','prev_symbol']:
  for a in str(getattr(row,col,'')).upper().split('|'):
   if a:aliases.setdefault(a,set()).add(row.symbol.upper())
rename={a:next(iter(s)) for a,s in aliases.items() if len(s)==1 and a not in official}
map0['Human gene name']=map0['Human gene name'].replace(rename);map0=map0[~map0['Human gene name'].duplicated(keep=False)]
meta=meta0.set_index('Sample ID');cols={s:s for s in raw.columns}
for aa,bb in plan['technical_pairs']:cols[bb]=aa[:-1];cols[aa]=aa[:-1]
X=raw.T.groupby(pd.Index([cols[s] for s in raw.columns]),sort=False).sum().T
lib=X.sum(axis=0);CPM=X.div(lib,axis=1)*1e6
m=map0.set_index('Gene stable ID')['Human gene name'];ids=[g for g in X.index if g in m]
CPM=CPM.loc[ids];CPM.index=m.loc[ids].tolist();assert CPM.index.is_unique
keep=(CPM.ge(1).sum(axis=1)>=math.ceil(.2*CPM.shape[1])) & (np.log2(CPM+1).std(axis=1)>1e-8)
CPM=CPM.loc[keep];E=np.log2(CPM+1).T;assert len(E)==26 and E.columns.is_unique
E.to_pickle(R/'inputs/expression.pkl');E.to_csv(R/'inputs/expression.tsv.gz',sep='\t');CPM.to_csv(R/'inputs/CPM.tsv.gz',sep='\t');map0.to_csv(R/'inputs/orthology_used.tsv',sep='\t',index=False)
tfs=sorted(set(use(R/'reference/hs_hgnc_tfs.txt').read_text().splitlines()));tfs=sorted({rename.get(t.upper(),t.upper()) for t in tfs}&set(E.columns));(R/'inputs/tf_names.txt').write_text('\n'.join(tfs)+'\n')
rows=[]
label={'Shrew':'shrew','BigBrownBat':'big brown bat','LittleBrownBat':'little brown bat','GuineaPig':'guinea pig','Porcupine':'porcupine','Chinchilla':'chinchilla','Chipmunk':'chipmunk','RedSquirrel':'red squirrel','FoxSquirrel':'fox squirrel','Beaver':'beaver','AfricanGrassRat':'grass rat','MeadowVole':'meadow vole','CottonRat':'cotton rat','WhiteFootedMouse':'white-footed mouse','DeerMouse':'deer mouse'}
for s in E.index:
 original=[k for k,v in cols.items() if v==s];sp=s.split('.')[0]
 rows.append(dict(sample=s,species=label[sp],source_columns=';'.join(original),technical_libraries=len(original),source_library_sum=int(lib[s])))
T=pd.DataFrame(rows);T.to_csv(R/'inputs/samples.tsv',sep='\t',index=False);assert T.species.nunique()==15
memb=pd.read_csv(use(P/'growth_regulators_2026-10-21_v1/results/wgcna/ma2016_full__raw/members.tsv'),sep='\t');memb.to_csv(R/'inputs/frozen_WGCNA_members.tsv',sep='\t',index=False)
# Prepare all uniquely named protein-coding human genes plus all expressed symbols, independent of mass/module membership.
gtf=use(P/'sizes5_execution_2026-10-16_v1/sources/Homo_sapiens.GRCh38.99.gtf.gz');genome=use(Path('/path/to/storage/human_cs6_embryo/track2/ref/GRCh38.fa'));use(Path(str(genome)+'.fai'))
genes={};trans={};exlen={};biotype={}
with gzip.open(gtf,'rt') as f:
 for line in f:
  if line.startswith('#'):continue
  c=line.rstrip().split('\t')
  if len(c)<9 or c[2] not in ['gene','transcript','exon']:continue
  d=dict(re.findall(r'(\w+) "([^"]+)"',c[8]));name=d.get('gene_name','').upper();name=rename.get(name,name);gid=d.get('gene_id','');tid=d.get('transcript_id')
  if not name:continue
  genes.setdefault(name,set()).add(gid)
  if c[2]=='gene':biotype[name]=d.get('gene_biotype','')
  if c[2]=='transcript' and tid:trans[tid]=(name,c[0] if c[0].startswith('chr') else 'chr'+c[0],int(c[3]),int(c[4]),c[6])
  if c[2]=='exon' and tid:exlen[tid]=exlen.get(tid,0)+int(c[4])-int(c[3])+1
best={}
for tid,(name,chrom,start,end,strand) in trans.items():
 if len(genes[name])!=1 or (biotype.get(name)!='protein_coding' and name not in E.columns):continue
 if not re.match(r'^chr([1-9]|1[0-9]|2[0-2]|X|Y)$',chrom):continue
 length=exlen.get(tid,0);v=(length,tid,chrom,start,end,strand)
 if name not in best or (length,tid)>(best[name][0],best[name][1]):best[name]=v
fai={c[0]:list(map(int,c[1:5])) for c in (line.split('\t') for line in Path(str(genome)+'.fai').read_text().splitlines())};comp=str.maketrans('ACGTN','TGCAN');genomef=genome.open('rb');regions=[]
for tag,(up,down) in plan['reference_windows'].items():
 n=0
 with (R/'reference'/f'human_{tag}.fa').open('w') as out:
  for name,(length,tid,chrom,start,end,strand) in sorted(best.items()):
   lo,hi=(start-1+up,start-1+down) if strand=='+' else (end-down,end-up)
   if chrom not in fai or lo<0 or hi>fai[chrom][0]:continue
   size,offset,bases,width=fai[chrom];p0=offset+lo//bases*width+lo%bases;p1=offset+(hi-1)//bases*width+(hi-1)%bases+1;genomef.seek(p0);seq=genomef.read(p1-p0).replace(b'\n',b'').replace(b'\r',b'').decode().upper();assert len(seq)==down-up
   if seq.count('N')>.2*len(seq):continue
   if strand=='-':seq=seq.translate(comp)[::-1]
   out.write('>'+name+'\n'+seq+'\n');n+=1;regions.append(dict(window=tag,gene=name,transcript=tid,chrom=chrom,start0=lo,end=hi,strand=strand,N_fraction=seq.count('N')/len(seq)))
 print('FASTA',tag,n,flush=True)
genomef.close();pd.DataFrame(regions).to_csv(R/'reference/promoter_regions.tsv.gz',sep='\t',index=False)
meme=use(P/'size_program_regulons_2026-10-21_v1/reference/JASPAR2024_one_per_symbol.meme');lines=meme.read_text().splitlines();motifs=[];annot=[];(R/'reference/motifs_cb').mkdir(exist_ok=True)
i=0
while i<len(lines):
 if lines[i].startswith('MOTIF '):
  mid,tf=lines[i].split()[1:3];i+=1
  while i<len(lines) and 'letter-probability matrix:' not in lines[i]:i+=1
  length=int(re.search(r'w=\s*(\d+)',lines[i]).group(1));mat=np.array([[float(x) for x in lines[i+j+1].split()] for j in range(length)]);counts=np.rint(mat*10000).astype(int)
  (R/'reference/motifs_cb'/f'{mid}.cb').write_text('>'+mid+'\n'+'\n'.join('\t'.join(map(str,row)) for row in counts)+'\n');motifs.append(mid)
  for gene in tf.upper().split('::'):
   gene=rename.get(gene,gene);annot.append({'#motif_id':mid,'gene_name':gene,'motif_similarity_qvalue':0.,'orthologous_identity':1.,'description':'JASPAR2024 direct TF-label assignment; heterodimer component if applicable; no native-mammal occupancy measurement'})
 i+=1
assert len(motifs)==816;pd.DataFrame(annot).to_csv(R/'reference/motif_annotations.tsv',sep='\t',index=False);(R/'reference/motifs.txt').write_text('\n'.join(motifs)+'\n')
out=dict(n_profiles=E.shape[0],n_species=T.species.nunique(),n_genes=E.shape[1],n_TFs=len(tfs),filter_min_samples=math.ceil(.2*E.shape[0]),E2F2_present='E2F2' in E.columns,E2F4_present='E2F4' in E.columns,motif_count=len(motifs),selected_WGCNA_coverage={mod:{'total':int(memb.module.eq(mod).sum()),'RNA_input':int(memb[memb.module.eq(mod)].feature.isin(E.columns).sum())} for mod in ['orange','violet','cyan','skyblue3']},code_sha256=sha(Path(__file__)),inputs=inputs,versions={n:md.version(n) for n in ['numpy','pandas','scipy','pyscenic','arboreto','ctxcore','scikit-learn','pyarrow']})
(R/'audit/preparation.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k not in ['inputs']},indent=2),flush=True)
