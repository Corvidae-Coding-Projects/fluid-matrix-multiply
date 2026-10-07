"""Build standalone and preview HTML, plus an optional inline fragment, from computed fields."""
from pathlib import Path
import argparse,base64,gzip,html

ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--inline',type=Path);a=p.parse_args()
template=(ROOT/'viewer-template.html').read_text()
data=(ROOT/'results/viewer-data.json').read_text()
script=(ROOT/'viewer.js').read_text()
encoded=base64.b64encode(gzip.compress(data.encode(),mtime=0)).decode()
library=base64.b64encode(gzip.compress((ROOT/'vendor/three.min.js').read_bytes(),mtime=0)).decode()
loader="(async function(){async function unzip(s){const bytes=Uint8Array.from(atob(s),c=>c.charCodeAt(0));return await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).text();}const lib=document.createElement('script');lib.textContent=await unzip('"+library+"');document.head.append(lib);window.FLUIDIC_CFD=JSON.parse(await unzip('"+encoded+"'));\n"
fragment=template.replace('<!-- CFD_DATA -->','').replace('<!-- CFD_SCRIPT -->','<script>'+loader+script+'\n})().catch(e=>{console.error(e);document.getElementById("cfd-status").textContent="The computed field data could not be loaded in this browser."});</script>')
fragment=fragment.replace('<script src="https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.min.js"></script>','')
if a.inline:
    if len(fragment.encode())>1_000_000:raise ValueError('Inline fragment exceeds 1 MB')
    a.inline.write_text(fragment)
offline=fragment
shell=(ROOT/'standalone-template.html').read_text().replace('/* CFD_STYLES */',(ROOT/'vendor/visualize.css').read_text()).replace('<!-- CFD_FRAGMENT -->',offline)
(ROOT/'fluidic-matrix-multiplier-3d.html').write_text(shell)
(ROOT/'preview.html').write_text((ROOT/'preview-template.html').read_text().replace('__CFD_FRAGMENT_HTML_ATTRIBUTE__',html.escape(fragment,quote=True)))
print('Standalone viewer:',ROOT/'fluidic-matrix-multiplier-3d.html')
if a.inline:print('Inline fragment:',a.inline)
