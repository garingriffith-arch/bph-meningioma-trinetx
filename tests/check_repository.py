#!/usr/bin/env python3
"""Basic offline repository hygiene and version-presence check; no patient data."""
from pathlib import Path
import re, sys
root = Path(__file__).resolve().parents[1]
stages = [p for p in (root / 'R').glob('*.R') if re.match(r'^\d\d_', p.name)]
nums = {int(p.name[:2]) for p in stages}
missing = sorted(set(range(28)) - nums)
assert not missing, f'Missing stage numbers: {missing}'
for name in ['README.md', '.gitignore', 'docs/analysis-pipeline.md', 'reporting/methods-source-verified.md']:
    assert (root/name).is_file(), name
for path in root.rglob('*'):
    if not path.is_file(): continue
    if any(x in path.parts for x in ('.git', '__pycache__')): continue
    assert path.suffix.lower() not in {'.duckdb','.sqlite','.parquet','.rds','.rdata'}, str(path)
    assert path.name not in {'medication.csv','diagnosis.csv','encounter.csv','patient.csv'}, str(path)
    if path.suffix.lower() == '.csv':
        assert path.parent.name == 'source-provenance', f'Unreviewed CSV: {path}'
print('PASS: 28 numbered stages present; no common restricted-data file types detected.')
print('NOTE: this does not inspect script scientific validity or substitute for manual release review.')
