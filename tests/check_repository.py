#!/usr/bin/env python3
from pathlib import Path
import re
root=Path(__file__).resolve().parents[1]
needed=['run_all.do','src/build_panel_from_sources.do','src/reporting/generate_report.do',
        'src/preprocessing/06_geography.do',
        'src/preprocessing/07_final_es_panel.do',
        'src/labor/steps/02_fit_weights.do',
        'src/labor/steps/07_joint_placebo_refits.do',
        'src/migration/01_census_v2025.do']
for fn in needed: assert (root/fn).is_file(),fn
assert not (root/'modules/legacy-inference').exists()
assert not (root/'src/preprocessing/06_scm.do').exists()
assert not (root/'src/preprocessing/08_reference_audit.do').exists()
files=[p for p in root.rglob('*.do') if 'archive' not in p.relative_to(root).parts]
for p in files:
 text=p.read_text()
 assert 'mc_w_legacy' not in text,p
 assert 'scm_donor_weight' not in text,p
 assert 'mc_w_emp' not in text,p
 assert 'modules/legacy-inference' not in text,p
 assert 'mc_raw_urate' not in text,p
 for m in re.finditer(r'^\s*(?:gen|generate)\s+(?:(?:byte|int|long|double|float|str\d+)\s+)?([a-zA-Z_][a-zA-Z0-9_]*)\s*=',text,re.M):
  assert len(m.group(1))<=32,(p,m.group(1))
main=(root/'run_all.do').read_text()
assert 'src/build_panel_from_sources.do' in main
assert 'run_multicontrol.do' in main
assert 'run_joint_placebos.do' in main
assert 'run_inference.do' not in main
fit=(root/'src/labor/steps/02_fit_weights.do').read_text()
assert 'mc_fit_joint_only' in fit
assert 'mc_w_joint' in fit
assert '2023' in fit
assert 'mc_w_legacy' not in fit
assert not list((root/'src').rglob('*.dta'))
print('PASS: no legacy modules, estimators, weight vectors or bundled analysis panels; one source-first pipeline.')
