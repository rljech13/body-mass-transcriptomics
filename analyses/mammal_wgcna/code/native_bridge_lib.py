from pathlib import Path
import numpy as np,pandas as pd,re,collections,json,hashlib,subprocess
from urllib.parse import unquote
def save(p,z):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
def makeindex(fa,out):
 d={};name=None;length=0;offset=0;bases=0;width=0
 with fa.open('rb') as f:
  while True:
   pos=f.tell();l=f.readline()
   if not l:
    if name is not None:d[name]=(length,offset,bases,width)
    break
   if l.startswith(b'>'):
    if name is not None:d[name]=(length,offset,bases,width)
    name=l[1:].split()[0].decode();length=0;offset=f.tell();bases=0;width=0
   else:
    z=len(l.rstrip(b'\r\n'));length+=z
    if not bases:bases=z;width=len(l)
 out.write_text(''.join(k+'\t'+'\t'.join(map(str,v))+'\n' for k,v in d.items()));return d
class Fasta:
 def __init__(self,p,ix):
  self.p=p;self.f=p.open('rb');self.ix={l.split('\t')[0]:tuple(map(int,l.rstrip().split('\t')[1:])) for l in ix.read_text().splitlines()} if ix.exists() else makeindex(p,ix)
 def seq(self,s,a,b):
  n,off,lb,lw=self.ix[s];a=max(0,int(a));b=min(int(b),n);out=[]
  while a<b:
   c=min(b-a,lb-a%lb);self.f.seek(off+(a//lb)*lw+a%lb);out.append(self.f.read(c));a+=c
  return b''.join(out).decode().upper()
def parse_annotation(p,gtf=False):
 genes={};counts=collections.Counter();tx=[]
 with p.open() as f:
  for l in f:
   if l.startswith('#'):continue
   a=l.rstrip().split('\t')
   if len(a)<9:continue
   at=dict(re.findall(r'(\w+) "([^"]*)"',a[8])) if gtf else {k:unquote(v) for k,v in (x.split('=',1) for x in a[8].split(';') if '=' in x)}
   if gtf:
    gid=at.get('gene_id','')
    if a[2]=='gene':counts[gid]+=1
    if a[2]!='transcript':continue
    eligible=at.get('gene_biotype')=='protein_coding' and at.get('transcript_biotype')=='protein_coding' and 'cds_start_NF' not in a[8] and 'mRNA_start_NF' not in a[8]
   else:
    if a[2] in ['gene','pseudogene']:
     m=re.search(r'GeneID:(\d+)',at.get('Dbxref',''))
     if m:genes[at.get('ID','')]=(m.group(1),at.get('gene_biotype',''));counts[m.group(1)]+=1
     continue
    if at.get('Parent','') not in genes:continue
    gid,biotype=genes[at['Parent']];eligible=a[2]=='mRNA' and biotype=='protein_coding' and at.get('partial')!='true'
   tx.append(dict(gene=gid,seq=a[0],strand=a[6],five=int(a[3])-1 if a[6]=='+' else int(a[4])-1,eligible=bool(eligible)))
 groups={}
 for z in tx:
  z['eligible'] &= counts[z['gene']]==1;k=(z['gene'],z['seq'],z['strand'],z['five'])
  if k in groups:groups[k]['eligible'] |= z['eligible']
  else:groups[k]=z
 return list(groups.values()),dict(gene_features=sum(counts.values()),unique_gene_ids=len(counts),duplicate_gene_ids=sum(v>1 for v in counts.values()),transcripts=len(tx),unique_gene_TSS=len(groups))
def paf_hits(p,queries,min_identity=.95,min_cov=.90):
 hits=collections.defaultdict(list)
 for l in p.read_text().splitlines():
  a=l.split('\t');tags={x[:2]:x[5:] for x in a[12:] if len(x)>5};q=a[0];info=queries[q];frac=(int(a[3])-int(a[2]))/int(a[1]);ident=int(a[9])/max(1,int(a[10]));center=info['center']
  if tags.get('tp')!='P' or int(a[11])<40 or ident<min_identity or frac<min_cov or not(int(a[2])<=center<int(a[3])):continue
  offset=center-int(a[2]) if a[4]=='+' else int(a[3])-1-center;qp=0;tp=int(a[7]);mapped=None
  for num,op in re.findall(r'(\d+)([MIDNSHP=X])',tags.get('cg','')):
   n=int(num)
   if op in ['M','=','X']:
    if qp<=offset<qp+n:mapped=tp+offset-qp;break
    qp+=n;tp+=n
   elif op in ['I','S']:
    if qp<=offset<qp+n:break
    qp+=n
   elif op in ['D','N']:tp+=n
  if mapped is None:continue
  strand=info['strand'] if a[4]=='+' else ('-' if info['strand']=='+' else '+');hits[q].append(dict(target=a[5],five=mapped,strand=strand,mapq=int(a[11]),identity=ident,coverage=frac))
 return {q:h[0] for q,h in hits.items() if len(h)==1}
def run_map(reference,query,paf,log):
 cmd=['/usr/bin/minimap2','-t','4','-x','asm10','-k','15','-w','10','-A','1','-B','4','-O','6,26','-E','2,1','-s','100','-m','100','-c','--secondary=yes','-N','5',str(reference),str(query)]
 with paf.open('w') as out,log.open('w') as err:r=subprocess.run(cmd,stdout=out,stderr=err)
 save(Path(str(paf)+'.execution.json'),dict(command=cmd,returncode=r.returncode,query_sha256=sha(query)));assert r.returncode==0
 return cmd
def quantify(raw,anchors,cols):
 ix={}
 for z in anchors:
  ix.setdefault((z['seq'],z['strand']),[]).append(z)
 for k,rr in ix.items():rr.sort(key=lambda z:z['five']);ix[k]=(np.array([z['five'] for z in rr]),rr)
 reason=collections.Counter();acc={};owners=[];X=raw[cols].to_numpy(float);assert np.isfinite(X).all() and (X>=0).all()
 def win(seq,st,m):
  pair=ix.get((seq,st))
  if pair is None:return None,'no_projected_TSS'
  xx,rr=pair;j=np.searchsorted(xx,m);ds=min(abs(m-xx[j-1]) if j else np.inf,abs(m-xx[j]) if j<len(xx) else np.inf)
  if ds>2000:return None,'outside_2000bp'
  near=np.flatnonzero(abs(abs(xx-m)-ds)<1e-9);gs={rr[k]['gene'] for k in near}
  if len(gs)!=1:return None,'different_gene_tie'
  return next(iter(gs)),'accepted' if any(rr[k]['eligible'] for k in near) else 'nearest_ineligible_no_reassignment'
 for i,z in enumerate(raw.itertuples()):
  seq=str(z.seqnames).split('|')[-1];m=(z.start+z.end)/2.;g,why=win(seq,z.strand,m);g1,w1=win(seq,z.strand,m-1)
  if z.start<1:why='source_boundary_quarantine'
  elif (g,why)!=(g1,w1):why='coordinate_origin_discordant'
  reason[why]+=1
  if why=='accepted':acc.setdefault(g,np.zeros(len(cols)));acc[g]+=X[i];owners.append(dict(peak_id=z.id,gene=g))
 out=pd.DataFrame.from_dict(acc,orient='index',columns=cols);out.index.name='gene_id';assert len({z['peak_id'] for z in owners})==len(owners);assert np.isclose(out.to_numpy().sum(),X[raw.id.isin([z['peak_id'] for z in owners])].sum(),rtol=1e-10,atol=1e-5)
 return out,owners,dict(reason)
