"""Read-only scientific acceptance of the three bounded CC single points.

This script never launches or alters a calculation. It requires all three
outputs and immutable source geometries before writing a separate audit.
"""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from gaussian_data import stationary_data, HARTREE_KCAL

ROOT=Path(__file__).resolve().parents[1]
NUMBER=r'[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def get(pattern,text):
    hits=re.findall(pattern,text,re.MULTILINE)
    assert hits,pattern
    return float(hits[-1])

def parse_output(text):
    assert 'ORCA TERMINATED NORMALLY' in text[-4000:],'Output incomplete or abnormal'
    assert 'SCF CONVERGED AFTER' in text and 'SCF NOT CONVERGED' not in text
    assert '--- The CCSD iterations have converged ---' in text
    assert 'Parameters used in iterative triples' in text
    for spin in ['AAA','BBB','AAB','ABB']:
        assert f'TRIPLES CORRECTION ({spin})' in text,spin
    assert 'Frozen core treatment                      ... chemical core (8 el)' in text
    assert 'Re-Reading QRO orbitals for the actual calculation' in text
    assert get(r'Number of Electrons\s+NEL\s+\.+\s+(\d+)',text)==31
    assert get(r'Multiplicity\s+Mult\s+\.+\s+(\d+)',text)==2
    assert get(r'Basis Dimension\s+Dim\s+\.+\s+(\d+)',text)==222
    # Parse the explicitly printed final CCSD residual, not localization/SCF.
    cc=text.split('OPEN-SHELL COUPLED CLUSTER ITERATIONS')[-1].split('--- The CCSD iterations have converged ---')[0]
    iterations=re.findall(r'^\s*(\d+)\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s*$',cc,re.MULTILINE)
    assert iterations
    last=iterations[-1]
    assert int(last[0])<100 and float(last[4])<1e-5
    largest=text.split('LARGEST PNO AMPLITUDES')[-1].split('Total # unrestricted pairs')[0]
    amplitudes=[float(v) for v in re.findall(r'\s('+NUMBER+r')\s*$',largest,re.MULTILINE)]
    assert amplitudes
    energy=get(r'FINAL SINGLE POINT ENERGY\s+('+NUMBER+')',text)
    ccsd=get(r'^E\(CCSD\)\s*\.\.\.\s*('+NUMBER+')',text)
    triples=get(r'Triples Correction \(T\)\s*\.\.\.\s*('+NUMBER+')',text)
    printed_ccsdt=get(r'^E\(CCSD\(T\)\)\s*\.\.\.\s*('+NUMBER+')',text)
    assert abs(energy-printed_ccsdt)<2e-9 and abs(energy-(ccsd+triples))<3e-9
    return {'E_hartree':energy,'E_CCSD_hartree':ccsd,'iterative_triples_correction_hartree':triples,
      'UHF_reference_S2':get(r'Expectation value of <S\*\*2>\s*:\s*('+NUMBER+')',text),
      'CCSD_linearized_S2':get(r'<S\*\*2>\(linearized\)\s*\.\.\.\s*('+NUMBER+')',text),
      'T1_diagnostic':get(r'T1 diagnostic\s*\.\.\.\s*('+NUMBER+')',text),
      'largest_printed_PNO_amplitude_abs':max(abs(v) for v in amplitudes),
      'CCSD_last_iteration':int(last[0]),'CCSD_last_delta_E_hartree':float(last[3]),
      'CCSD_last_residual':float(last[4]),'chemical_core_electrons_frozen':8,
      'normal_termination':True,'iterative_triples':True,'all_four_spin_triples_blocks_present':True,'basis_functions':222,
      'electrons':31,'multiplicity':2,
      'reference':'UHF-derived quasi-restricted orbitals for open-shell DLPNO; S2 diagnostics retain their separately labeled definitions'}
