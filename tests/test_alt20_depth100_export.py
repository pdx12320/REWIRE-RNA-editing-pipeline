import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('exporter', Path(__file__).resolve().parents[1] / 'pipeline/scripts/rna/export_alt20_depth100.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class ExportTests(unittest.TestCase):
    def setUp(self):
        self.samples = ['R1', 'R2', 'R3', 'R4']
        self.site = dict(chromosome='chr1', position='100', ref='C', alt='T', editing_rate='.2')
        self.rows = [dict(self.site, sample=s, depth='100', alt_reads='20') for s in self.samples]

    def test_inclusive_boundaries(self):
        self.assertEqual(len(m.select([self.site], self.rows, self.samples)), 1)

    def test_one_rep_fails(self):
        self.rows[0]['alt_reads'] = '19'
        self.assertEqual(m.select([self.site], self.rows, self.samples), [])

    def test_depth_fails(self):
        self.rows[0]['depth'] = '99'
        self.assertEqual(m.select([self.site], self.rows, self.samples), [])

    def test_missing_is_not_zero(self):
        with self.assertRaises(ValueError): m.select([self.site], self.rows[:-1], self.samples)

    def test_duplicate_fails(self):
        with self.assertRaises(ValueError): m.select([self.site], self.rows + self.rows[:1], self.samples)

    def test_pooled_rate_rejected(self):
        self.site['editing_rate'] = '.5'
        with self.assertRaises(ValueError): m.select([self.site], self.rows, self.samples)

if __name__ == '__main__': unittest.main()
