"""Accept six fixed-geometry Gaussian comparisons from immutable raw evidence."""
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from gaussian_data import BOHR_A, HARTREE_KCAL, fchk_field

ROOT = Path(__file__).resolve().parents[1]
Z = [6, 6, 6, 1, 1, 1, 1, 1, 6, 1, 1]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def populations(text, kind):
    marker = ('Mulliken charges and spin densities:' if kind == 'Mulliken' else
              'Hirshfeld charges, spin densities, dipoles, and CM5 charges using IRadAn=')
    assert marker in text
    block = text.rsplit(marker, 1)[1]
    rows = []
    started = False
    for line in block.splitlines():
        m = re.match(r'^\s*(\d+)\s+([CH])\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', line)
        if m:
            started = True
            rows.append({'atom': int(m[1]), 'element': m[2], 'charge': float(m[3]), 'spin': float(m[4])})
        elif started:
            break
    assert [r['atom'] for r in rows] == list(range(1, 12)), (kind, len(rows))
    assert [6 if r['element'] == 'C' else 1 for r in rows] == Z
    assert abs(sum(r['spin'] for r in rows) - 1) < 2e-5
    assert abs(sum(r['charge'] for r in rows)) < 2e-5
    return rows

def differences(energies):
    return {k: (energies[a] - energies[b]) * HARTREE_KCAL
            for k, a, b in [('forward', 'TS', 'R'), ('reverse', 'TS', 'P'), ('reaction', 'P', 'R')]}
