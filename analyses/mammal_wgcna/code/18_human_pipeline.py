from pathlib import Path
import subprocess,time,json,hashlib,datetime
O=Path(__file__).resolve().parents[1];py='/path/to/user-home/chicken/kleiber_birds_gemcat/results/celloracle_pilot_2026-09-29/venv/bin/python'
while not (O/'audit/human_promoters.json').exists():time.sleep(5)
rr=[]
for args in [['17_human_own5_union.py'],['07_pwm_statistics.py','human'],['07_pwm_statistics.py','human_full']]:
 p=O/'code'/args[0];before=hashlib.sha256(p.read_bytes()).hexdigest()
 with (O/'logs'/('human_pipeline_'+'_'.join(args)+'.log')).open('w') as out:r=subprocess.run([py,str(p)]+args[1:],stdout=out,stderr=subprocess.STDOUT)
 rr.append(dict(script=args[0],args=args[1:],before_sha256=before,after_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),returncode=r.returncode,ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 (O/'audit/human_pipeline_execution.json').write_text(json.dumps(rr,indent=2));assert r.returncode==0
