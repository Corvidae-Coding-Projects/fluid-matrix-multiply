"""Check every packaged file against MANIFEST.sha256, using Python's standard library."""
from pathlib import Path
import hashlib
import sys

root = Path(__file__).resolve().parent
manifest = root / 'MANIFEST.sha256'
if not manifest.is_file():
    sys.exit('Missing MANIFEST.sha256. Extract the complete ZIP first.')

failures = []
checked = 0
for line in manifest.read_text(encoding='utf-8').splitlines():
    expected, name = line.split('  ', 1)
    path = root / name
    if not path.is_file():
        failures.append(f'MISSING: {name}')
        continue
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        failures.append(f'CHANGED OR DAMAGED: {name}')
    checked += 1
if failures:
    sys.exit('\n'.join(failures))
print(f'PASS: all {checked} packaged files match their SHA-256 checksums.')
