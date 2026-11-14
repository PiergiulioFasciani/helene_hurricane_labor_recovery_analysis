"""Ensure stale estimates are never published when their stage failed."""
from pathlib import Path
import json
import subprocess
import sys
import tempfile
root=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as td:
    fixture=Path(td)
    reports=fixture/'outputs/reports';reports.mkdir(parents=True)
    audit=fixture/'outputs/audit';audit.mkdir(parents=True)
    (audit/'STAGE_STATUS.csv').write_text('stage,status,return_code\njoint_multioutcome,FAILED,459\njoint_placebos,FAILED,9\n')
    tables=fixture/'outputs/labor/tables';tables.mkdir(parents=True)
    (tables/'period_effects_WITH_HOLM.csv').write_text('model,outcome,window,display_beta\njoint,mc_yemp,immediate,987654321\n')
    subprocess.run([sys.executable,str(root/'src/reporting/make_report.py'),'--root',str(fixture)],check=True)
    report=(reports/'HELENE_COMPLETE_RESULTS.txt').read_text()
    manifest=json.loads((reports/'graphs/manifest.json').read_text())
    assert '987654321' not in report
    assert manifest['inputs']==[]
    assert len(manifest['figures'])==1
    assert 'not PASSED' in report
print('PASS: failed stages cannot publish stale estimates or numerical figures.')
