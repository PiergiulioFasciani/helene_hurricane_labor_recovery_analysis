from pathlib import Path
p=Path(__file__).resolve().parents[1]
main=(p/'run_all.do').read_text()
assert 'capture erase' in main and 'panel' in main
builder=(p/'src/build_panel_from_sources.do').read_text()
assert '06_geography.do' in builder and '07_final_es_panel.do' in builder
assert '06_scm.do' not in builder
geog=(p/'src/preprocessing/06_geography.do').read_text()
assert 'joint_training_n==24' in geog
assert 'donor_weight' not in geog
panel=(p/'src/preprocessing/07_final_es_panel.do').read_text()
assert 'joint_eligible_donor' in panel
assert 'es_donor_pool' in panel
labor=(p/'src/preprocessing/05_labor.do').read_text()
assert 'https://api.economicdata.dev/v1/laus/download/laus.csv' in labor
assert 'LAUS_mirror_all_counties.csv' in labor
rep=(p/'src/reporting/generate_report.do').read_text()
assert 'STALE FILE IGNORED' in rep
assert 'modules/legacy-inference' not in rep
print('PASS: raw-only rebuild, four-outcome-only donor screening and provenance.')
