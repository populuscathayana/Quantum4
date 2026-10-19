import hashlib
import json
from pathlib import Path
import re
import unittest
from crosscheck_parser import differences,populations

ROOT=Path(__file__).resolve().parents[1]

class PublicationTests(unittest.TestCase):
    def test_manifest_paths_and_types(self):
        rows=json.loads((ROOT/'data/manifest.json').read_text())['files']
        self.assertEqual(len(rows),len({r['path'] for r in rows}))
        for r in rows:
            p=Path(r['path']);self.assertFalse(p.is_absolute());self.assertNotIn('..',p.parts)
            self.assertNotIn(p.suffix,{'.chk','.gbw','.cube','.key','.pem'})
            self.assertEqual(hashlib.sha256((ROOT/p).read_bytes()).hexdigest(),r['sha256'])
    def test_readme_local_links(self):
        for name in ['README.md','ENVIRONMENT.md','POSTPROCESS.md','THIRD_PARTY_NOTICES.md']:
            for link in re.findall(r'\]\(([^)]+)\)',(ROOT/name).read_text()):
                if '://' not in link and not link.startswith('#'):
                    self.assertTrue((ROOT/link.split('#')[0]).exists(),link)
    def test_energy_cycle(self):
        d=differences({'R':-10,'TS':-9.99,'P':-10.005})
        self.assertAlmostEqual(d['forward']-d['reverse'],d['reaction'])
    def test_reject_truncated_population(self):
        with self.assertRaises(AssertionError):
            populations('Mulliken charges and spin densities:\n 1 C 0.000000 1.000000\n','Mulliken')

if __name__=='__main__':unittest.main()
