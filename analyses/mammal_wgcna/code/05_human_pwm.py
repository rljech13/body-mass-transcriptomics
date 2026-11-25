from pathlib import Path
import json,hashlib,time,subprocess
import numpy as np,pandas as pd,requests
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
O=Path(__file__).resolve().parents[1];R=O.parent;V=R/'nonbird_regulatory_followup_2026-10-21_v1/reference/human_promoters';P=O/'reference/human';P.mkdir(exist_ok=True);cache=P/'cache';cache.mkdir(exist_ok=True)
trackp=V/'ucsc_ncbiRefSeqSelect.json';track=json.loads(trackp.read_text());data=track['ncbiRefSeqSelect']
if isinstance(data,dict):data=[a for li in data.values() for a in li]
ann=pd.DataFrame(data);ann=ann[ann.chrom.str.match(r'^chr([1-9]|1[0-9]|2[0-2]|X|Y)$')].drop_duplicates(['name2','chrom','strand','txStart','txEnd']);ann=ann[~ann.name2.duplicated(keep=False)].set_index('name2')
genes=set()
for panel in ['ma2016','gse105160','wagner']:genes.update(pd.read_csv(O/'inputs'/f'{panel}.log.tsv',sep='\t',usecols=[0]).iloc[:,0])
metadata=[];fail=[];lock=Lock();last=[0.]
for g in sorted(genes&set(ann.index)):
 t=ann.loc[g];sgn=1 if t.strand=='+' else -1;tss=int(t.txStart)+1 if sgn==1 else int(t.txEnd);a,b=(tss-800,tss+199) if sgn==1 else (tss-199,tss+800)
 if a<1:continue
 metadata.append(dict(gene=g,chrom=t.chrom[3:],start=a,end=b,strand=sgn,region=f'{t.chrom[3:]}:{a}..{b}:{sgn}',reference_transcript=t['name']))
def one(z):
 p=cache/(hashlib.sha256(z['region'].encode()).hexdigest()+'.json');old=V/'sequence_cache'/p.name
 if old.exists():rec=json.loads(old.read_text());origin='reused_previous_receipt_cache'
 elif p.exists():rec=json.loads(p.read_text());origin='new_cache'
 else:
  rec=None
  for attempt in range(3):
   try:
    with lock:
     delay=.14-(time.monotonic()-last[0])
     if delay>0:time.sleep(delay)
     last[0]=time.monotonic()
    r=requests.get('https://api.genome.ucsc.edu/getData/sequence',params={'genome':'hg38','chrom':'chr'+z['chrom'],'start':z['start']-1,'end':z['end']},timeout=(10,30));r.raise_for_status();rec=r.json();assert len(rec['dna'])==1000;p.write_text(r.text);origin='new_UCSC_API';break
   except Exception as e:
    if attempt==2:return z,None,str(e)[:200]
    time.sleep(1)
 s=rec['dna'].upper()
 if z['strand']==-1:s=s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
 if len(s)!=1000 or s.count('N')>10:return z,None,'sequence_quality'
 z.update(sequence_name=z['gene'],GC=(s.count('G')+s.count('C'))/len(s),CpG=s.count('CG')/999,N_fraction=s.count('N')/len(s),length=len(s),sequence_sha256=hashlib.sha256(s.encode()).hexdigest(),retrieval=origin)
 return z,s,None
rows=[];fasta=P/'promoters.fa'
with fasta.open('w') as out,ThreadPoolExecutor(max_workers=10) as pool:
 for i,(z,s,why) in enumerate(pool.map(one,metadata)):
  if s is None:fail.append(dict(gene=z['gene'],reason=why))
  else:rows.append(z);out.write('>'+z['gene']+'\n'+s+'\n')
  if i%200==0:print('HUMAN',i+1,'/',len(metadata),'fail',len(fail),flush=True)
pd.DataFrame(rows).to_csv(P/'promoters.tsv',sep='\t',index=False);pd.DataFrame(fail).to_csv(P/'failures.tsv',sep='\t',index=False)
motif=O/'reference/candidate_JASPAR2024.meme';cmd=['/path/to/user-home/miniconda3/envs/meme/bin/fimo','--text','--thresh','1e-4',str(motif),str(fasta)]
with (P/'fimo.tsv').open('w') as out,(O/'logs/fimo_human.log').open('w') as err:r=subprocess.run(cmd,stdout=out,stderr=err)
assert r.returncode==0
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
(O/'audit/human_promoters.json').write_text(json.dumps(dict(n_requested=len(genes),n_unambiguous_annotation=len(metadata),n_sequences=len(rows),n_failures=len(fail),window='hg38 RefSeqSelect transcript TSS -800/+200, 1000nt',scope='human reference for mammal gene universe; no native interspecies conservation or active fibroblast TSS claim',annotation_sha256=sha(trackp),sequence_fasta_sha256=sha(fasta),motif_sha256=sha(motif),script_sha256=sha(Path(__file__)),command=cmd,returncode=r.returncode),indent=2));print('HUMAN DONE',len(rows),flush=True)
