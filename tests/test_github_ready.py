"""Offline packaging checks: run locally or on a clean GitHub checkout."""
from pathlib import Path
import hashlib
import subprocess

ROOT = Path(__file__).resolve().parents[1]
def git(*args):
    return subprocess.run(['git', *args], cwd=ROOT, check=True, capture_output=True, text=True).stdout

# The only generated file published on purpose. It changes on every pipeline
# run, so it is deliberately absent from SHA256SUMS.txt.
GENERATED_RECEIPT = 'outputs/audit/REGRESSIONS.txt'

candidates = set(git('ls-files', '--cached', '--others', '--exclude-standard').splitlines())
assert candidates, 'No distributable files detected'
for name in candidates:
    p = Path(name)
    assert p.name != 'config.local.do', name
    assert not p.name.startswith('.env'), name
    assert p.suffix not in {'.dta', '.ster', '.gph', '.smcl', '.log', '.zip', '.pem', '.key'}, name
    if p.parts[0] == 'data':
        assert p.name == '.gitkeep', name
    if p.parts[0] == 'outputs':
        assert p.name == '.gitkeep' or name in (
            'outputs/README.md', GENERATED_RECEIPT), name
    if p.parts[0] == 'archive':
        assert name == 'archive/README.md', name
for folder in ['data/raw', 'data/raw/LAUS', 'data/raw/QCEW', 'data/raw/census',
               'data/interim', 'data/input', 'outputs/audit',
               'outputs/preprocessing/audits', 'outputs/preprocessing/panels',
               'outputs/preprocessing/logs', 'outputs/labor/data',
               'outputs/labor/tables', 'outputs/labor/figures',
               'outputs/labor/logs', 'outputs/labor/placebos',
               'outputs/migration', 'outputs/reports/graphs']:
    assert folder + '/.gitkeep' in candidates, folder
manifest = {}
for line in (ROOT / 'SHA256SUMS.txt').read_text().splitlines():
    expected, name = line.split('  ', 1)
    assert name in candidates, name
    assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
    manifest[name] = expected
assert set(manifest) == candidates - {'SHA256SUMS.txt', GENERATED_RECEIPT}, \
    'Source manifest and upload file list differ'
print(f'PASS: {len(candidates)} upload files; required placeholders retained; no private/generated files; hashes verified.')
