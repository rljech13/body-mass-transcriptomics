from pathlib import Path
import csv,gzip,json,hashlib,statistics,math

w=Path(__file__).resolve().parents[1];root=w.parents[1]
source=root/'outputs/mammal_SCENIC_2026-10-23/selected'
targets=root/'outputs/mammal_SCENIC_regulon_screen_2026-10-23/data'
inputs=[]
def read(p,compressed=False):
    inputs.append(dict(reference=str(p.relative_to(root)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
    with (gzip.open(p,'rt') if compressed else p.open()) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rows):
    with (w/'data'/name).open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)

selected=[]
for tf in ['ERG','FOXP4']:
    edges=read(targets/f'{tf}_all_targets.tsv')
    assert all(x['TF']==tf and x['target']!=tf for x in edges)
    edges.sort(key=lambda x:(-float(x['importance']),x['target']))
    for rank,x in enumerate(edges[:10],1):selected.append(dict(TF=tf,gene=x['target'],rank=rank,importance=float(x['importance']),row_id=tf+'::'+x['target']))
assert len(selected)==20 and len({x['row_id'] for x in selected})==20
samples=read(source/'inputs/samples.tsv');meta={x['sample']:x for x in samples}
traits=read(source/'inputs/species_traits.tsv');traits.sort(key=lambda x:(float(x['mass_g']),x['species']))
assert len(traits)==15 and len(meta)==26
expression=read(source/'inputs/expression.tsv.gz',compressed=True)
sample_field=list(expression[0])[0]
assert len(expression)==26 and set(x[sample_field] for x in expression)==set(meta)
genes={x['gene'] for x in selected};assert genes<=set(expression[0])
sample_data=[];long=[];row_data=[]
for sel in selected:
    gene=sel['gene'];values=[]
    for sample in expression:
        name=sample[sample_field];v=float(sample[gene]);assert math.isfinite(v) and v>=0
        sample_data.append(dict(TF=sel['TF'],gene=gene,rank=sel['rank'],sample=name,species=meta[name]['species'],log2_CPM_plus1=v))
    for tr in traits:
        sp=tr['species'];v=[float(x[gene]) for x in expression if meta[x[sample_field]]['species']==sp]
        assert v;values.append(statistics.mean(v))
    mean=statistics.mean(values);sd=statistics.stdev(values);assert sd>1e-12
    row_data.append(dict(sel,row_mean=mean,row_sample_sd=sd))
    for pos,(tr,v) in enumerate(zip(traits,values),1):
        long.append(dict(sel,species=tr['species'],species_order=pos,mass_g=float(tr['mass_g']),log10_mass=math.log10(float(tr['mass_g'])),n_sample_groups=sum(x['species']==tr['species'] for x in samples),species_mean_log2_CPM_plus1=v,row_z=(v-mean)/sd))
assert len(long)==300 and len(sample_data)==520
for sel in selected:
    z=[x['row_z'] for x in long if x['row_id']==sel['row_id']]
    assert abs(statistics.mean(z))<1e-12 and abs(statistics.stdev(z)-1)<1e-12
write('selected_targets.tsv',row_data);write('sample_gene_expression.tsv',sample_data);write('species_gene_expression.tsv',long)
write('species_order.tsv',[dict(species=x['species'],mass_g=float(x['mass_g']),log10_mass=math.log10(float(x['mass_g']))) for x in traits])
(w/'audit/inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
plan=dict(date='2026-10-24',primary_seed=20261023,study='Ma2016',selection='Top ten nonself targets per ERG/FOXP4 by saved motif-pruned GRNBoost2 importance; descending importance, alphabetical gene as tie breaker; no body-mass or gene correlation filtering.',normalization='Frozen SCENIC input log2(CPM+1), arithmetic mean across sample-label groups within each species, then row z over 15 equally weighted species using sample SD (ddof=1).',ordering='Species by ascending mass; rows by saved importance within fixed ERG then FOXP4 blocks. No clustering, sign flipping, winsorization or clipping.',rows=20,columns=15,cells=300,sample_values=520,max_abs_z=max(abs(x['row_z']) for x in long),minimum_z=min(x['row_z'] for x in long),maximum_z=max(x['row_z'] for x in long),new_SCENIC=False,new_association_tests=False,figure_approval='pending',presentation_modified=False,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(w/'audit/preparation.json').write_text(json.dumps(plan,indent=2)+'\n');print(json.dumps(plan,indent=2))
