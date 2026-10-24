#!/usr/bin/env python3
"""Create a locally configured copy of historical R scripts; never alter source."""
import argparse
from pathlib import Path

OLD = {
    'C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R': 'script_dir',
    'C:/TriNetX/bph_m_work': 'work_dir',
    'C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/Data/BPH': 'raw_dir',
}


def norm(x):
    return str(Path(x).expanduser().resolve()).replace('\\','/')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw-dir', required=True, help='Authorized TriNetX input CSV directory')
    p.add_argument('--work-dir', required=True, help='Local DuckDB/work directory')
    p.add_argument('--out-dir', required=True, help='New directory for configured R scripts')
    a = p.parse_args()
    repo = Path(__file__).resolve().parents[1]
    out = Path(a.out_dir).resolve()
    if out == (repo/'R').resolve() or repo/'R' in out.parents:
        p.error('Output must be outside the archived R/ source directory')
    if out.exists() and any(out.iterdir()):
        p.error('Output directory must be new or empty to avoid silently overwriting a prior analysis')
    out.mkdir(parents=True,exist_ok=True)
    target = {'raw_dir':norm(a.raw_dir),'work_dir':norm(a.work_dir),'script_dir':norm(out)}
    mapping = {k: target[key] for k,key in OLD.items()}
    changes=[]
    for file in sorted((repo/'R').glob('*.R')):
        source=file.read_text(encoding='utf-8-sig')
        text=source
        hit=[]
        for old,new in mapping.items():
            count=text.count(old)
            if count:
                text=text.replace(old,new)
                hit.append(f'{old}: {count} occurrence(s)')
        (out/file.name).write_text(text,encoding='utf-8')
        if hit: changes.append((file.name,hit))
    note = out/'CONFIGURATION_NOTES.txt'
    note.write_text('Path-only configured copy of archived R scripts. This does not run the study.\n'+
        '\n'.join(f'{name}: '+', '.join(hits) for name,hits in changes)+
        '\nReview dependencies, file schemas and available licensed source data before execution.\n',encoding='utf-8')
    print(f'Prepared {len(list(out.glob("*.R")))} scripts at {out}')
    print(f'Configuration audit: {note}')

if __name__ == '__main__': main()
