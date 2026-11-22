#!/usr/bin/env bash
# Copyright 2026 Roby Joehanes; GNU GPL version 3.
# Regenerate synthetic fixtures only; requires bcftools and Python 3.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
destination=${1:?Supply a new output directory}
mkdir "$destination"
python3 - "$root" "$destination" <<'PY'
import base64,hashlib,json,pathlib,subprocess,sys
root,out=map(pathlib.Path,sys.argv[1:])
source=(root/'test/resources/variant-reference/sparse-dictionary.vcf').read_bytes().replace(b'\r\n',b'\n')
proof=[]
for name,data in [('sparse-dictionary',source),('sparse-contigs',source.replace(b'ID=1,length=1000,IDX=0',b'ID=1,length=1000,IDX=3').replace(b'ID=2,length=1000,IDX=1',b'ID=2,length=1000,IDX=7'))]:
    vcf=out/(name+'.vcf');bcf=out/(name+'.bcf');vcf.write_bytes(data)
    command=['bcftools','view','--no-version','-Ob','-o',str(bcf),str(vcf)]
    subprocess.run(command,check=True)
    (out/(name+'.bcf.base64')).write_bytes(base64.b64encode(bcf.read_bytes())+b'\n')
    proof.append(dict(name=name,command=command,sha256=hashlib.sha256(bcf.read_bytes()).hexdigest()))
(out/'generation.json').write_text(json.dumps(dict(bcftools=subprocess.check_output(['bcftools','--version'],text=True),files=proof),indent=2)+'\n')
PY
