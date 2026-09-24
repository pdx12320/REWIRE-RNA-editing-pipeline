import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location(
    'strand_exporter', Path(__file__).resolve().parents[1] /
    'pipeline/scripts/rna/export_positive_strand.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class StrandExportTests(unittest.TestCase):
    def setUp(self):
        self.fields = list(m.KEY) + ['gene', 'sequence_101nt', 'editing_rate']
        self.rows = [dict(chromosome='chr1', position=str(i), ref=ref, alt=alt,
                          gene='GENE_A|GENE_B', sequence_101nt='A'*50+'C'+'G'*50,
                          editing_rate='0.18542718019999999')
                     for i, (ref, alt) in enumerate([('C', 'T'), ('G', 'A')], 1)]

    def test_both_strands_preserve_every_original_value(self):
        fields, rows = m.annotate(self.fields, self.rows)
        self.assertEqual(fields[:3], ['chromosome', 'position', 'strand'])
        self.assertEqual([r['strand'] for r in rows], ['+', '-'])
        for original, result in zip(self.rows, rows):
            self.assertEqual(original, {k: result[k] for k in self.fields})
        self.assertNotIn('strand', self.rows[0])

    def test_existing_correct_strand_is_preserved(self):
        fields, rows = m.annotate(self.fields, self.rows)
        self.assertEqual(m.annotate(fields, rows), (fields, rows))

    def test_conflicting_or_unknown_existing_strand_is_rejected(self):
        for strand in ['-', '.', '']:
            with self.subTest(strand=strand), self.assertRaises(ValueError):
                m.annotate(self.fields + ['strand'], [{**self.rows[0], 'strand': strand}])

    def test_not_an_arbitrary_variant_strand_predictor(self):
        for ref, alt in [('A', 'G'), ('T', 'C'), ('C', 'A'), ('c', 't')]:
            with self.subTest(ref=ref, alt=alt), self.assertRaises(ValueError):
                m.annotate(self.fields, [{**self.rows[0], 'ref': ref, 'alt': alt}])

    def test_duplicate_site(self):
        with self.assertRaises(ValueError):
            m.annotate(self.fields, self.rows + self.rows[:1])

    def test_bad_sequence(self):
        for seq in ['C'*100, 'A'*101, 'N'*50+'C'+'G'*50]:
            with self.subTest(seq=seq), self.assertRaises(ValueError):
                m.annotate(self.fields, [{**self.rows[0], 'sequence_101nt': seq}])

    def test_optional_sequence_and_empty_table(self):
        self.assertEqual(m.annotate(list(m.KEY), [dict(zip(m.KEY, ['chr1','1','G','A']))])[1][0]['strand'], '-')
        self.assertEqual(m.annotate(self.fields, [])[1], [])

    def test_bad_header_or_row(self):
        for fields, rows in [(None, []), (self.fields + ['gene'], []), (['gene'], []),
                             (self.fields, [{**self.rows[0], 'position': '0'}]),
                             (self.fields, [{**self.rows[0], 'gene': None}]),
                             (self.fields, [{**self.rows[0], None: ['extra']}])]:
            with self.subTest(fields=fields, rows=rows), self.assertRaises(ValueError):
                m.annotate(fields, rows)

    def test_file_roundtrip_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            source, out = Path(directory)/'in.tsv', Path(directory)/'out.tsv'
            with source.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=self.fields, delimiter='\t')
                writer.writeheader()
                writer.writerows(self.rows)
            original = source.read_bytes()
            self.assertEqual(m.export(source, out), {'+': 1, '-': 1})
            result = out.read_bytes()
            with self.assertRaises(FileExistsError):
                m.export(source, out)
            with self.assertRaises(FileExistsError):
                m.export(source, source)
            self.assertEqual(source.read_bytes(), original)
            self.assertEqual(out.read_bytes(), result)
            with out.open() as handle:
                self.assertEqual(len(list(csv.DictReader(handle, delimiter='\t'))), 2)

    def test_invalid_input_does_not_create_output(self):
        with tempfile.TemporaryDirectory() as directory:
            source, out = Path(directory)/'in.tsv', Path(directory)/'out.tsv'
            source.write_text('chromosome\tposition\tref\talt\nchr1\t1\tA\tG\n')
            with self.assertRaises(ValueError):
                m.export(source, out)
            self.assertFalse(out.exists())


if __name__ == '__main__':
    unittest.main()
