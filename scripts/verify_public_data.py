"""Offline verification of the published files and original scientific results."""
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from gaussian_data import stationary_data, fchk_field, BOHR_A, HARTREE_KCAL
from parse_gaussian import parse_text
from cc_parser import parse_output
from crosscheck_parser import populations, differences

ROOT = Path(__file__).resolve().parents[1]

def load(path): return json.loads((ROOT/path).read_text())
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def verify():
    manifest=load('data/manifest.json'); sources={r['source']:r for r in manifest['files']}
    for row in manifest['files']:
        path=ROOT/row['path']; assert path.is_relative_to(ROOT)
        assert path.stat().st_size==row['bytes'] and digest(path)==row['sha256'], row['path']
        if path.suffix in {'.json','.gjf','.inp','.out','.tex','.py'}:
            text=path.read_text()
            assert '/Users/yeqingyang' not in text and '/home/class-qc/qc26' not in text
    def source(path, original_hash=None):
        row=sources[path]
        if original_hash: assert row['original_sha256']==original_hash,path
        return ROOT/row['path']
    # All outputs: the original failed PM7 scan remains failed, not relabeled.
    outputs=[r for r in manifest['files'] if r['path'].endswith('calc.out')]
    for r in outputs:
        text=(ROOT/r['path']).read_text()
        if '/cc-' in r['path']:
            parse_output(text)
        elif '/pm7-scan/' in r['path']:
            assert 'Error termination' in text
        else:
            assert 'Normal termination of Gaussian' in text[-1500:] and 'Error termination' not in text,r['path']
    thermo=load('analysis/thermochemistry-20261019/audit.json')['states']
    geometries={}
    for s,r in thermo.items():
        out=source(r['source'],r['source_sha256'])
        fchk=out.with_suffix('.fchk'); d=stationary_data(fchk); p=parse_text(out.read_text())
        assert digest(fchk)==r['fchk_sha256']
        assert d['electrons']==31 and d['alpha']==16 and d['beta']==15 and len(d['Z'])==11
        assert len(d['frequencies_cm_1'])==27 and sum(d['frequencies_cm_1']<0)==(s=='TS')
        assert np.max(np.abs(d['frequencies_cm_1']-p['frequencies_cm_1']))<.001
        assert abs(d['energy_hartree']-r['E_hartree'])<1e-8
        assert abs(p['thermo_hartree']['G']-r['G_corrected_hartree'])<1e-8
        for k in ['E','G_corrected']:
            energy = r['E_hartree'] if k=='E' else r['G_corrected_hartree']
            ref = thermo['R']['E_hartree'] if k=='E' else thermo['R']['G_corrected_hartree']
            assert abs((energy-ref)*HARTREE_KCAL-r['relative_kcal_mol'][k])<1e-7
        stable=source(r['stability']['source'],r['stability']['source_sha256']).read_text()
        assert 'The wavefunction is stable under the perturbations considered.' in stable
        geometries[s]=d['xyz_A']
    dft=load('analysis/report-review-20261019/dft-crosscheck-audit.json')
    for model,m in dft['methods'].items():
        e={}
        for s,r in m['states'].items():
            text=source(r['output'],r['output_sha256']).read_text()
            ftext=source(r['fchk'],r['fchk_sha256']).read_text()
            assert text.count('Normal termination of Gaussian')==2
            assert 'The wavefunction is already stable.' in text
            energies=[float(v) for v in re.findall(r'SCF Done:.*?=\s*(-?\d+\.\d+)',text)]
            assert len(energies)==2 and abs(energies[0]-energies[1])<1e-8
            e[s]=fchk_field(ftext,'Total Energy'); assert abs(e[s]-r['E_hartree'])<1e-9
            assert fchk_field(ftext,'Number of basis functions')==222
            xyz=fchk_field(ftext,'Current cartesian coordinates').reshape(11,3)*BOHR_A
            assert np.max(np.abs(xyz-geometries[s]))<1e-7
            for kind in ['Mulliken','Hirshfeld']: assert populations(text,kind)==r['populations'][kind]
        for k,v in differences(e).items(): assert abs(v-m['electronic_differences_kcal_mol'][k])<1e-7
    cc=load('analysis/report-review-20261019/cc-benchmark-audit.json'); e={}
    for s,r in cc['states'].items():
        d=parse_output(source(r['output'],r['output_sha256']).read_text())
        for k in ['E_hartree','UHF_reference_S2','CCSD_linearized_S2','T1_diagnostic']:
            assert abs(d[k]-r[k])<1e-9,(s,k)
        e[s]=d['E_hartree']
    for k,v in differences(e).items(): assert abs(v-cc['CC_electronic_differences_kcal_mol'][k])<1e-7
    review=load('analysis/report-review-20261019/cc-diagnostics-review.json')
    assert review['reviewed_flags']==cc['review_flags']
    spin=load('analysis/irc-20261019/spin-audit.json')['rows']; assert len(spin)==56
    for r in spin:
        p=parse_text(source(r['source'],r['source_sha256']).read_text())
        assert p['last_Mulliken_populations']==r['Mulliken_populations']
        assert abs(sum(a['spin'] for a in p['last_Mulliken_populations'])-1)<2e-5
        assert abs(p['last_energy_hartree']-r['energy_hartree'])<1e-8
    pm7=load('analysis/report-review-20261019/pm7-complete-audit.json')
    points=pm7['joined_points']; assert len(points)==16
    original=parse_text((ROOT/'data/calculations/pm7-scan/calc.out').read_text())
    tail=parse_text((ROOT/'data/calculations/pm7-tail/calc.out').read_text())
    accepted=original['completed_optimizations']+tail['completed_optimizations']
    assert len(accepted)==16
    for a,b in zip(accepted,points): assert abs(a['energy_hartree']-b['energy_hartree'])<1e-10
    report=load(manifest['report_audit'])
    assert digest(ROOT/manifest['report_pdf'])==report['pdf_sha256']
    assert digest(ROOT/'experiment-report/radical_isomerization.tex')==report['TeX_sha256']
    assert report['pages']==16 and report['visually_reviewed_pages']==list(range(1,17))
    link=load('analysis/publication-20261019/pdf-link-audit.json')
    assert link['pdf_sha256']==report['pdf_sha256']
    assert link['cover_unchanged'] and link['scientific_content_unchanged']
    assert link['repository_link_PDF_pages']==[15]
    return {'verified_files':len(manifest['files']),'verified_outputs':len(outputs),
            'PM7_points':16,'IRC_spin_points':56,'report_pages':16,'status':'passed'}

if __name__=='__main__': print(json.dumps(verify(),indent=2))
