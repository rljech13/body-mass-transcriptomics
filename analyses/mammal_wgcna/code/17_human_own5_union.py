from common import *
import gzip
V=R/'nonbird_regulatory_followup_2026-10-21_v1/reference/human_promoters';P=O/'reference/human';Q=O/'reference/human_plus_own5';Q.mkdir(exist_ok=True);genes=set(sum(json.loads((O/'config/fixed_sets.json').read_text())['mammal'].values(),[]))
def fa(p):
 d={};k=None
 for l in p.read_text().splitlines():
  if l.startswith('>'):k=l[1:].split()[0];d[k]=''
  elif k:d[k]+=l.strip()
 return d
t=pd.read_csv(P/'promoters.tsv',sep='\t');u=pd.read_csv(V/'promoters.tsv',sep='\t');seq=fa(P/'promoters.fa');old=fa(V/'promoters.fa')
for g in set(seq)&set(old):assert seq[g]==old[g],g
missing=genes-set(t.gene);extra=u[u.gene.isin(missing)].copy()
for g in extra.gene:seq[g]=old[g]
extra['CpG']=[old[g].count('CG')/999 for g in extra.gene];extra['length']=1000;extra['sequence_name']=extra.gene;extra['retrieval']='reused_exact_previous_hg38_own5'
tt=pd.concat([t,extra[t.columns.intersection(extra.columns)]],ignore_index=True);tt.to_csv(Q/'promoters.tsv',sep='\t',index=False)
with (Q/'promoters.fa').open('w') as f:
 for g in tt.gene:f.write('>'+g+'\n'+seq[g]+'\n')
F=pd.read_csv(P/'fimo.tsv',sep='\t',comment='#');G=pd.read_csv(V/'fimo.tsv.gz',sep='\t',comment='#');mi=set(pd.read_csv(O/'reference/candidate_motifs.tsv',sep='\t').motif_id);G=G[G.sequence_name.isin(extra.gene)&G.motif_id.isin(mi)];pd.concat([F,G],ignore_index=True).to_csv(Q/'fimo.tsv',sep='\t',index=False)
oldreceipt=json.loads((V.parent.parent/'audit/motif_receipt.json').read_text());motifhash=json.loads((O/'audit/native_promoters.json').read_text())['motif_source_sha256'];assert motifhash in oldreceipt['input_sha256'].values()
(O/'audit/human_own5_union.json').write_text(json.dumps(dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),missing_from_shared=int(len(missing)),added_own5_genes=extra.gene.tolist(),all20_own5_present=genes<=set(tt.gene),shared_sequences_identical=True,frozen_motif_source_hash_agrees=True,source_promoters_sha256=hashlib.sha256((V/'promoters.tsv').read_bytes()).hexdigest(),source_fimo_sha256=hashlib.sha256((V/'fimo.tsv.gz').read_bytes()).hexdigest(),n_sequences=len(tt),scope='human reference common-annotation background plus exact cached own5 promoters; not complete native mammal CRE'),indent=2));print('OWN5 UNION',extra.gene.tolist(),flush=True)
