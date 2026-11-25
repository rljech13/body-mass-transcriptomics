from pathlib import Path
import subprocess,json,hashlib,time,datetime
O=Path(__file__).resolve().parents[1];py='/path/to/user-home/chicken/kleiber_birds_gemcat/results/celloracle_pilot_2026-09-29/venv/bin/python'
while not (O/'results/wgcna_full_runs.tsv').exists():time.sleep(5)
receipts=[]
for name in ['12_statistics_full.py','13_own5_full_original.py']:
 p=O/'code'/name;before=hashlib.sha256(p.read_bytes()).hexdigest();start=datetime.datetime.now(datetime.timezone.utc).isoformat();cmd=[py,str(p)]
 with (O/'logs'/f'{p.stem}.log').open('w') as out:r=subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT)
 receipts.append(dict(file=name,command=cmd,started_utc=start,ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sha256_before=before,sha256_after=hashlib.sha256(p.read_bytes()).hexdigest(),returncode=r.returncode,common_dependency_sha256=hashlib.sha256((O/'code/common.py').read_bytes()).hexdigest()))
 (O/'audit/full_pipeline_execution.json').write_text(json.dumps(receipts,indent=2));assert r.returncode==0 and receipts[-1]['sha256_before']==receipts[-1]['sha256_after']
