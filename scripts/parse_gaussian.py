"""Gaussian scalar and accepted-optimization parser for this eleven-atom task."""
import argparse
import json
from pathlib import Path
import re

NUM=r'[-+]?\d+(?:\.\d*)?(?:[DdEe][-+]?\d+)?'
def f(x): return float(x.replace('D','E').replace('d','e'))

def parse_text(text):
    lines=text.splitlines()
    geometry=None;energy=None;method=None;spin=None;points=[];freq=[];last_pop=None
    for i,line in enumerate(lines):
        if 'Input orientation:' in line or 'Standard orientation:' in line:
            j=i+1;dashes=0
            while j<len(lines) and dashes<2:
                if re.fullmatch(r'\s*-+\s*',lines[j]):dashes+=1
                j+=1
            atoms=[]
            while j<len(lines):
                fields=lines[j].split()
                if len(fields)!=6 or not fields[0].isdigit():break
                atoms.append({'index':int(fields[0]),'Z':int(fields[1]),'xyz_A':[f(v) for v in fields[3:]]});j+=1
            if atoms:geometry=atoms
        match=re.search(r'SCF Done:\s+E\(([^)]+)\)\s*=\s*('+NUM+')',line)
        if match:method=match.group(1);energy=f(match.group(2))
        match=re.search(r'\bEnergy=\s*('+NUM+r')\s+NIter=',line)
        if match:method='semiempirical';energy=f(match.group(1))
        match=re.search(r'S\*\*2 before annihilation\s+('+NUM+r'),\s+after\s+('+NUM+')',line)
        if match:spin={'before':f(match.group(1)),'after':f(match.group(2))}
        if 'Optimization completed.' in line or 'Optimization completed on the basis of negligible forces.' in line:
            if geometry is None or energy is None:raise ValueError('Optimization without geometry/energy')
            points.append({'completion_line_1based':i+1,'geometry':geometry,'energy_hartree':energy,'method':method,'S2':spin})
        if 'Frequencies --' in line:freq.extend(f(x) for x in line.split('--',1)[1].split())
        if re.match(r'\s*Mulliken charges and spin densities:',line):
            pop=[]
            for row in lines[i+1:i+20]:
                fields=row.split()
                if len(fields)==4 and fields[0].isdigit() and fields[1] in {'C','H'}:
                    pop.append({'atom':int(fields[0]),'element':fields[1],'charge':f(fields[2]),'spin':f(fields[3])})
                elif pop:break
            if pop:last_pop=pop
    thermo={}
    for label,key in [('Zero-point correction','ZPE'),('Thermal correction to Gibbs Free Energy','G_correction'),('Sum of electronic and zero-point Energies','E_ZPE'),('Sum of electronic and thermal Enthalpies','H'),('Sum of electronic and thermal Free Energies','G')]:
        values=re.findall(re.escape(label)+r'\s*=\s*('+NUM+')',text)
        if values:thermo[key]=f(values[-1])
    electrons=re.findall(r'(\d+) alpha electrons\s+(\d+) beta electrons',text)
    return {'normal_termination':'Normal termination of Gaussian' in '\n'.join(lines[-12:]),
      'error_termination':'Error termination' in text,'completed_optimizations':points,
      'last_geometry':geometry,'last_energy_hartree':energy,'last_S2':spin,
      'electron_counts':list(map(int,electrons[-1])) if electrons else None,
      'frequencies_cm_1':freq,'imaginary_frequencies_cm_1':[v for v in freq if v<0],
      'last_Mulliken_populations':last_pop,'thermo_hartree':thermo}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output');a=p.parse_args()
    print(json.dumps(parse_text(Path(a.output).read_text(errors='replace')),indent=2))
