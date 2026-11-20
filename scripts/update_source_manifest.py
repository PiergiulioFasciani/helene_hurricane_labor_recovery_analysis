"""Refresh checksums for files Git would include, excluding the manifest itself."""
from pathlib import Path
import hashlib
import subprocess

root = Path(__file__).resolve().parents[1]
result = subprocess.run(
    ['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'],
    cwd=root, check=True, capture_output=True)
paths = sorted(set(p.decode() for p in result.stdout.split(b'\0') if p))
lines = []
for name in paths:
    if name in ('SHA256SUMS.txt', 'outputs/audit/REGRESSIONS.txt'):
        continue  # the manifest itself, and a generated receipt that changes every run
    path = root / name
    if not path.is_file():
        raise SystemExit(f'Missing upload file: {name}')
    lines.append(hashlib.sha256(path.read_bytes()).hexdigest() + '  ' + name)
(root / 'SHA256SUMS.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
print(f'Updated source manifest: {len(lines)} files.')
