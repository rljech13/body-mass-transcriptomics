from pathlib import Path
import sys,json,hashlib,subprocess,re,collections
import pandas as pd,numpy as np
from native_bridge_lib import Fasta,parse_annotation,quantify,sha
O=Path(__file__).resolve().parents[1];R=O.parent;H=Path('/path/to/user-home/chicken/kleiber_birds_gemcat');B=R/'sizes6_three_gaps_2026-10-19_v1';ASM=R.parent/'ref/source_assemblies'
plan=json.loads((O/'config/analysis_plan.initial.json').read_text());print("plan keys",list(plan),flush=True)
# Exact initial TF list, kept explicit to make compatibility with plan schema auditable.
candidates='STAT5A STAT5B FOXO1 FOXO3 TEAD1 TEAD2 TEAD3 TEAD4 SMAD2 SMAD3 SMAD4 MYC MYCN MAX SREBF1 SREBF2 TFEB TFE3 E2F1 E2F2 E2F3 E2F4 PLAGL1 RUNX2 SOX9 THRA THRB ESR1 ESR2 ATF4 CEBPG'.split()
motifsrc=H/'results/oxphos_species_thread/story_cre_motifs/method_compare/JASPAR2024_one_per_symbol.meme'
s=motifsrc.read_text();parts=re.split(r'(?=^MOTIF )',s,flags=re.M);kept=[];mm=[]
for block in parts[1:]:
 head=block.splitlines()[0].split();mid,tf=head[1:3]
 if set(tf.split('::')) & set(candidates):kept.append(block);mm.append(dict(motif_id=mid,TF=tf))
motif=O/'reference/candidate_JASPAR2024.meme';motif.write_text(parts[0]+''.join(kept));pd.DataFrame(mm).to_csv(O/'reference/candidate_motifs.tsv',sep='\t',index=False)
reg=pd.read_csv(B/'02_measurements/sample_registry.tsv',sep='\t');mapping=pd.read_csv(B/'02_measurements/native_ortholog7_mapping.tsv',sep='\t',dtype=str)
meta=[];receipts=[];allfa=O/'reference/bird_native_promoters.fa';f=allfa.open('w')
for sp in ['Zebra','Quail','Chicken','Duck','Peacock','Turkey','Emu']:
 print('SEQUENCE',sp,flush=True);rawpath=H/'data/dropbox_20260918/exp_tab'/f'{sp}_fib_epi_exp.txt';raw=pd.read_csv(rawpath,sep='\t');cols=reg[reg.species.eq(sp)]['column'].tolist()
 if sp in ['Peacock','Turkey','Emu']:
  ap=B/'02_measurements/native_bridge'/sp/'source_TSS_anchors.tsv';anchors=pd.read_csv(ap,sep='\t').to_dict('records');ownerspath=B/'02_measurements/native_bridge'/sp/'peak_owners.tsv'
  owners=pd.read_csv(ownerspath,sep='\t',dtype=str);reasons={'reused_native_bridge_owners':len(owners)}
 else:
  ap=ASM/sp/'genomic.gff';anchors,inventory=parse_annotation(ap)
  out,owners0,reasons=quantify(raw,anchors,cols);owners=pd.DataFrame(owners0).astype(str)
  old=pd.read_csv(B/'02_measurements'/f'native_TSS_{sp}_all_gene_by_library.tsv',sep='\t',index_col=0);old.index=old.index.astype(str)
  assert set(out.index)==set(old.index),(sp,'gene registry mismatch')
  assert np.allclose(out.sort_index().to_numpy(),old.sort_index().to_numpy(),rtol=1e-9,atol=1e-7),(sp,'native expression replay mismatch')
 own=owners.merge(raw,on='peak_id' if 'peak_id' in raw else 'id') if False else owners.merge(raw,left_on='peak_id',right_on='id',validate='one_to_one')
 fa=ASM/sp/('genome_localnames.fa' if sp=='Peacock' else 'genome.fa');ix=O/'reference'/f'{sp}.source.fai';source=Fasta(fa,ix);canonical={k.split('|')[-1]:k for k in source.ix}
 gid2symbol=mapping.set_index(sp+'_GeneID').symbol.to_dict()
 for context in ['fib_high_glucose','fib_F12','epi_source_label']:
  cc=reg[reg.species.eq(sp)&reg.context.eq(context)]['column'].tolist()
  if not cc:continue
  z=own[own.gene.isin(gid2symbol)].copy();z['mean_peak_expression']=z[cc].mean(axis=1);z=z.sort_values(['gene','mean_peak_expression','peak_id'],ascending=[True,False,True]).drop_duplicates('gene')
  for row in z.itertuples():
   symbol=gid2symbol[row.gene];seq=str(row.seqnames).split('|')[-1];center=int(np.floor((row.start+row.end)/2))-1
   a,b=(center-300,center+100) if row.strand=='+' else (center-99,center+301)
   if seq not in canonical or a<0:continue
   dna=source.seq(canonical[seq],a,b)
   if row.strand=='-':dna=dna.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
   if len(dna)!=400 or dna.count('N')>4:continue
   name=sp+'__'+context+'__'+symbol;f.write('>'+name+'\n'+dna+'\n')
   meta.append(dict(sequence_name=name,species=sp,context=context,gene=symbol,gene_id=row.gene,peak_id=row.peak_id,seq=seq,start0=a,end0=b,strand=row.strand,length=len(dna),GC=(dna.count('G')+dna.count('C'))/len(dna),CpG=dna.count('CG')/(len(dna)-1),N_fraction=dna.count('N')/len(dna),mean_peak_expression=row.mean_peak_expression,sequence_sha256=hashlib.sha256(dna.encode()).hexdigest()))
 receipts.append(dict(species=sp,source_fasta=str(fa),source_fasta_sha256=sha(fa),annotation=str(ap),annotation_sha256=sha(ap),raw_sha256=sha(rawpath),assignment_counts=reasons))
f.close();pd.DataFrame(meta).to_csv(O/'reference/bird_native_promoters.tsv',sep='\t',index=False)
(O/'audit/native_promoters.json').write_text(json.dumps(dict(inputs=receipts,motif_source_sha256=sha(motifsrc),motif_subset_sha256=sha(motif),n_motifs=len(mm),n_sequences=len(meta),script_sha256=sha(Path(__file__)),window='400nt around dominant CAGE cluster midpoint; origin1 for sequence extraction, origin0/1 agreement for gene assignment; not a validated single-base TSS'),indent=2))
cmd=['/path/to/user-home/miniconda3/envs/meme/bin/fimo','--text','--thresh','1e-4',str(motif),str(allfa)]
with (O/'reference/bird_candidate_fimo.tsv').open('w') as out,(O/'logs/fimo_bird.log').open('w') as err:r=subprocess.run(cmd,stdout=out,stderr=err)
(O/'audit/fimo_bird_execution.json').write_text(json.dumps(dict(command=cmd,returncode=r.returncode,fasta_sha256=sha(allfa),motif_sha256=sha(motif),script_sha256=sha(Path(__file__))),indent=2));assert r.returncode==0
print('FIMO DONE',len(meta),len(mm),flush=True)
