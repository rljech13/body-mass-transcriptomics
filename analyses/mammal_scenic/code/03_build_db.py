from pathlib import Path
import subprocess,json,hashlib,time
R=Path(__file__).resolve().parents[1];P=R.parent/'planaria_pyscenic_2026-10-21_v1/reference';tool=P/'create_cistarget_tools/create_cistarget_motif_databases.py';py='/path/to/user-home/miniconda3/envs/create_cistarget_databases/bin/python'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for tag in ['500up100down','10kb_each_side']:
 t=time.time();prefix=R/'reference'/f'human_{tag}';fa=Path(str(prefix)+'.fa');args=[py,str(tool),'-f',str(fa),'-M',str(R/'reference/motifs_cb'),'-m',str(R/'reference/motifs.txt'),'-o',str(prefix),'-c',str(P/'cbust'),'-t','24','-s','20261023']
 (R/'audit'/f'build_{tag}.json').write_text(json.dumps(dict(command=args,tool_sha256=sha(tool),cbust_sha256=sha(P/'cbust'),fasta_sha256=sha(fa),code_sha256=sha(Path(__file__))),indent=2))
 with (R/'logs'/f'build_{tag}.log').open('w') as f:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
 print('DATABASE_DONE',tag,time.time()-t,flush=True)
