from pathlib import Path
import json,hashlib,tarfile,subprocess,importlib.metadata as md
R=Path(__file__).resolve().parents[1];O=R/'selected_export';O.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();files=[]
for sub in ['code','config','audit']:
 files.extend(p for p in (R/sub).glob('*') if p.is_file() and p.suffix not in ['.tar','.gz'])
files.extend(R/'inputs'/n for n in ['expression.pkl','expression.tsv.gz','samples.tsv','species_traits.tsv','tf_names.txt','orthology_used.tsv'])
for seed in [20261023,20261024,20261025]:
 for n in ['GRN_before.json','GRN_summary.json','SCENIC_summary.json','prune_before.json','coexpression_modules.tsv','regulon_inventory.tsv','regulon_edges.tsv.gz','regulons.pkl','motif_enrichment.csv','AUCell.tsv']:files.append(R/'results'/str(seed)/n)
files.extend(p for p in (R/'results/analysis').glob('*') if p.is_file())
for n in ['summary.json','module_regulon_results.tsv','regulon_mass_associations.tsv']:files.append(R/'results/analysis_v1_filtered_family'/n)
files.extend(R/'logs'/n for n in ['independent_R_v1.log','independent_R.log'])
prefixes=[('/path/to/storage/kleiber_birds_gemcat/results/','run:'),('/path/to/storage/chicken/Dosage_compensation_proj/','owner:Dosage_compensation_proj/'),('/path/to/storage/human_cs6_embryo/','owner:human_cs6_embryo/'),('/path/to/user-home/chicken/kleiber_birds_gemcat/','owner:legacy_Kleiber/'),('/path/to/user-home/miniconda3/','owner:conda/'),('/path/to/storage/','owner:data/'),('/path/to/user-home/','owner:home/')]
manifest=[]
for p in sorted(set(files)):
 if not p.exists():continue
 rel=p.relative_to(R);dst=O/rel;dst.parent.mkdir(parents=True,exist_ok=True);b=p.read_bytes();original=sha(p)
 if rel.parts[0]!='code' and p.suffix in ['.json','.txt','.log','.sha256']:
  s=b.decode()
  for old,new in prefixes:s=s.replace(old,new)
  b=s.encode()
 dst.write_bytes(b);manifest.append(dict(name=str(rel),source_sha256=original,preserved_sha256=sha(dst),bytes=len(b),sanitized=sha(dst)!=original))
remote=[]
for sub in ['inputs','reference','results']:
 for p in sorted((R/sub).rglob('*')):
  if p.is_file():remote.append(dict(locator='run:ma2016_pyscenic_2026-10-23_v1/'+str(p.relative_to(R)),sha256=sha(p),bytes=p.stat().st_size))
(O/'audit/owner_retained_file_manifest.json').write_text(json.dumps(remote,indent=2))
for p in [R/'config/existing_reference.json']:
 cfg=json.loads(p.read_text());refs=[]
 for key in ['database','annotation']:
  z=Path(cfg[key]);locator=str(z)
  for old,new in prefixes:locator=locator.replace(old,new)
  refs.append(dict(locator=locator,sha256=sha(z),bytes=z.stat().st_size))
(O/'audit/existing_reference_hashes.json').write_text(json.dumps(refs,indent=2))
cmd=subprocess.run(['git','-C',str(R.parents[1]),'rev-parse','HEAD'],text=True,capture_output=True);(O/'audit/runtime_and_git.json').write_text(json.dumps(dict(versions={n:md.version(n) for n in ['numpy','pandas','scipy','pyscenic','arboreto','ctxcore','scikit-learn','pyarrow','dask','distributed']},git_HEAD=cmd.stdout.strip() if cmd.returncode==0 else None,is_git_repository=cmd.returncode==0,source_policy='Code snapshots retain exact bytes; provider and run audit paths sanitized; original/preserved hashes both listed'),indent=2))
(O/'audit/export_manifest.json').write_text(json.dumps(manifest,indent=2));tar=R/'selected_export.tar.gz'
with tarfile.open(tar,'w:gz') as t:
 for p in sorted(O.rglob('*')):
  if p.is_file():t.add(p,arcname=str(p.relative_to(O)))
print(json.dumps(dict(snapshot_files=len(manifest),archive_sha256=sha(tar),archive_bytes=tar.stat().st_size)))
