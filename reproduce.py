"""Offline verification and regenerated numerical tables/figures; no QM launch."""
import argparse
import csv
import datetime as dt
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from verify_public_data import verify,load

def table(path,rows):
    with path.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser();p.add_argument('--figures',action='store_true')
    p.add_argument('--animations',action='store_true');p.add_argument('--tex',action='store_true')
    p.add_argument('--output',type=Path);a=p.parse_args()
    result=verify();print(json.dumps(result,indent=2))
    dest=a.output or ROOT/'reproduced'/dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    dest=dest.resolve();dest.mkdir(parents=True,exist_ok=False)
    (dest/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    table(dest/'electronic-method-comparison.csv',load('analysis/report-review-20261019/method-comparison-data.json')['rows'])
    table(dest/'pm7-scan.csv',[{k:r[k] for k in ['C1_C2_A','energy_hartree','relative_kcal_mol','origin']}
        for r in load('analysis/report-review-20261019/pm7-complete-audit.json')['joined_points']])
    table(dest/'thermochemistry.csv',[{'state':s,'E_hartree':r['E_hartree'],'G_hartree':r['G_corrected_hartree'],
        **r['relative_kcal_mol']} for s,r in load('analysis/thermochemistry-20261019/audit.json')['states'].items()])
    table(dest/'irc-spin.csv',[{'q_sqrt_amu_bohr':r['q_sqrt_amu_bohr'],'E_hartree':r['energy_hartree'],
        'relative_E_kcal_mol':r['relative_E_kcal_mol'],'S2':r['S2'],**{f'C{i+1}_spin':r['Mulliken_populations'][i]['spin'] for i in [0,1,2,8]},
        **r['distances_A']} for r in load('analysis/irc-20261019/spin-audit.json')['rows']])
    if a.figures or a.animations:
        from plot_results import render
        render(dest,a.animations)
    if a.tex:
        for _ in range(2):
            subprocess.run(['xelatex','-interaction=nonstopmode','-halt-on-error',f'-output-directory={dest}',
                            'radical_isomerization.tex'],cwd=ROOT/'experiment-report',check=True)
    print(f'Regenerated outputs: {dest}')

if __name__=='__main__':main()
