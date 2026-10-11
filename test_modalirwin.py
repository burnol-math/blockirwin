# -*- coding: utf-8 -*-
"""Regression tests for the GPT-6 implementation of Burnol's modal series.

Check selected results against the independent blockirwin.py method,
plus recurrence caches, diagnostics, input handling, and the CLI.
The mathematical formulas are due to J.-F. Burnol; see
README-modalirwin.md and the repository CC BY-SA 4.0 LICENSE.
"""
import io
from pathlib import Path
import re
import subprocess
import sys
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch

from blockirwin import blockirwin
from modalirwin import modalirwin
from compare_modalirwin import compare_irwin


class TestModalExpansion(unittest.TestCase):
    def coefficient_records(self, output, *, positive_count):
        """Validate and parse one-field-per-line coefficient diagnostics."""
        records = []
        constant = None
        current = None
        labels = {'n', 'budget digits', 'working digits', 'T', 'lambda', 'L', 'R', 'C', "C'"}
        value_column = len('working digits') + 3
        # Exclude the verbose summary, which follows the mode records.
        coefficient_text = output.split('Modal summary:', 1)[0]
        for line in coefficient_text.splitlines():
            if not line.strip():
                continue
            match = re.fullmatch(r"(constant|n|budget digits|working digits|lambda|L|R|C'|C|T)( +)= (\S+)", line)
            self.assertIsNotNone(match, f'Invalid coefficient line: {line!r}')
            label, _, value = match.groups()
            if label == 'constant':
                self.assertIsNone(constant)
                self.assertIsNone(current)
                constant = value
                continue
            self.assertIn(label, labels)
            self.assertEqual(line.index(value), value_column,
                             f'Unaligned field: {line!r}')
            if label == 'n':
                current = {'n': value}
                records.append(current)
            else:
                self.assertIsNotNone(current, line)
                self.assertNotIn(label, current)
                current[label] = value
        return constant, records

    def check(self, base, word, k, digits=12, **opts):
        reference = blockirwin(base, word, k, digits)
        computed = modalirwin(base, word, k, digits, **opts)
        self.assertEqual(reference, computed, (base, word, k))

    def test_one_digit_kempner(self):
        self.check(2, '0', 0)
        self.check(10, '0', 0)
        self.check(2, '1', 0)

    def test_all_zero_overlapping_kempner(self):
        self.check(2, '00', 0, 12)
        self.check(3, '00', 0, 10)

    def test_kempner_with_leading_zero(self):
        self.check(2, '01', 0, 14)
        self.check(2, '10', 0, 12)

    def test_positive_count_and_overlaps(self):
        self.check(2, '11', 1, 12)
        self.check(10, '42', 1, 14)
        self.check(10, '42', 3, 14)
        self.check(3, '12', 2, 12)

    def test_list_digits_and_explicit_ell(self):
        self.check(16, [15, 15], 2, 10, ell=2)
        self.check(2, '01', 1, 10, ell=4)

    def test_rounding_and_precision(self):
        self.check(2, '01', 0, 20)
        self.check(2, '1', 1, 14)

    def test_comparison_function_and_coimport(self):
        with redirect_stdout(io.StringIO()) as output:
            result = compare_irwin(10, '9', 1, 10)
        self.assertEqual(result['blockirwin'], result['modalirwin'])
        self.assertIn('elapsed', output.getvalue())
        self.assertIn('modalirwin', output.getvalue())
        self.assertRegex(output.getvalue(), r'(?m)^base\s+= 10$')
        self.assertRegex(output.getvalue(), r'(?m)^occurrences\s+= 1$')
        self.assertRegex(output.getvalue(), r'(?m)^blockirwin result\s+= \S+$')
        self.assertRegex(output.getvalue(), r'(?m)^modalirwin result\s+= \S+$')

    def test_longer_words_and_stability(self):
        self.check(2, '000', 0, 12)
        self.check(2, '0110', 0, 11)
        self.check(2, '0110', 0, 35)
        self.check(2, '010', 2, 12)
        self.check(3, '010', 0, 12)
        self.check(3, '121', 1, 12)

    def test_higher_precision(self):
        self.check(2, '01', 0, 45)
        self.check(2, '11', 1, 45)
        self.check(10, '42', 3, 60)

    def test_arbitrary_base_digit_sequence(self):
        self.check(37, [1, 2], 1, 12)

    def test_cached_triangular_coefficients(self):
        from mpmath import mp
        from modalirwin import _ModalEngine
        with mp.workdps(65):
            engine = _ModalEngine(2, [0, 1, 1, 0], 3, need_kempner=True)
            for which in ('u', 'v'):
                for row in (1, 2, 5, 9, 13):
                    for col in range(row + 1):
                        d, e = engine._de(row, col, which)
                        degree = row - col
                        t = mp.power(2, -(row + 1))
                        expected_e = mp.fsum(sign*t**length*integer**degree
                                             for length, integer, sign in
                                             (engine._uterms if which == 'u'
                                              else engine._vterms)[engine.b:])
                        expected_d = (int(row == col) + expected_e +
                                      mp.fsum(-t*digit**degree
                                              for digit in range(engine.b)))
                        self.assertTrue(mp.almosteq(d, expected_d))
                        self.assertTrue(mp.almosteq(e, expected_e))

    def test_truncation_and_mode_limit(self):
        self.assertEqual(modalirwin(2, '01', 0, 18, trunc=True),
                         blockirwin(2, '01', 0, 18, trunc=True))
        with self.assertRaises(ArithmeticError):
            modalirwin(2, '00', 0, 12, max_modes=2)

    def test_modalcoeffs_positive_count(self):
        from mpmath import mp
        output = io.StringIO()
        with redirect_stdout(output):
            result = modalirwin(2, '01', 2, 16, modalcoeffs=True)
        self.assertEqual(result, blockirwin(2, '01', 2, 16))
        constant, records = self.coefficient_records(output.getvalue(),
                                                     positive_count=True)
        self.assertIsNone(constant)
        self.assertTrue(records)
        terms = []
        for number, values in enumerate(records, start=1):
            self.assertEqual(values['n'], str(number))
            self.assertEqual(set(values), {'n','budget digits','working digits','lambda','L','R','C',"C'",'T'})
            self.assertTrue(values['T'].startswith(('+', '-')))
            with mp.workdps(50):
                lam, left, right = [mp.mpf(values[x]) for x in ('lambda','L','R')]
                coeff = mp.mpf(values['C'])
                coeff_prime = mp.mpf(values["C'"])
                term = mp.mpf(values['T'])
                self.assertTrue(mp.almosteq(coeff, left*right, rel_eps=mp.mpf('1e-14')))
                self.assertTrue(mp.almosteq(coeff_prime, coeff/lam, rel_eps=mp.mpf('1e-14')))
                self.assertTrue(mp.almosteq(term, coeff*lam, rel_eps=mp.mpf('1e-14')))
                terms.append(term)
        with mp.workdps(45):
            self.assertLess(abs(mp.fsum(terms) - mp.mpf(result)), mp.mpf('1e-14'))

    def test_modalcoeffs_kempner_constant_and_terms(self):
        from mpmath import mp
        for word, has_constant in [('00', True), ('01', False)]:
            with self.subTest(word=word):
                output = io.StringIO()
                with redirect_stdout(output):
                    result = modalirwin(2, word, 0, 14, modalcoeffs=True)
                self.assertEqual(result, blockirwin(2, word, 0, 14))
                constant, records = self.coefficient_records(output.getvalue(),
                                                             positive_count=False)
                self.assertEqual(constant, '1.0' if has_constant else None)
                self.assertTrue(records)
                for number, values in enumerate(records, start=1):
                    self.assertEqual(values['n'], str(number))
                    self.assertEqual(set(values), {'n', 'budget digits', 'working digits', 'T'})
                    self.assertTrue(values['T'].startswith(('+', '-')))
                with mp.workdps(40):
                    total = (mp.mpf(constant) if has_constant else mp.zero) + mp.fsum(
                        mp.mpf(values['T']) for values in records)
                    self.assertLess(abs(total - mp.mpf(result)), mp.mpf('1e-12'))

    def test_modalcoeffs_verbose_and_special_empty_sum(self):
        output = io.StringIO()
        with redirect_stdout(output):
            result = modalirwin(10, '42', 3, 10, verbose=True, modalcoeffs=True)
        self.assertEqual(result, blockirwin(10, '42', 3, 10))
        text = output.getvalue()
        constant, records = self.coefficient_records(text, positive_count=True)
        self.assertIsNone(constant)
        self.assertTrue(records)
        self.assertIn('Modal summary:', text)
        self.assertRegex(text, r'(?m)^estimated inclusive-tail bound\s+= \S+$')
        self.assertRegex(text, r'(?m)^threshold\s+= \S+$')
        self.assertIn('Tail bound is estimated, inclusive', text)
        self.assertNotIn('mode    1: term', text)
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(modalirwin(2, '1', 0, 10, modalcoeffs=True), '0')
        self.assertIn('no modal terms', output.getvalue())

    def test_modalcoeffs_cli(self):
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).with_name('modalirwin.py')),
             '2', '01', '1', '10', '--modalcoeffs'],
            capture_output=True, text=True, check=True)
        constant, records = self.coefficient_records(
            '\n'.join(completed.stdout.splitlines()[:-1]),
            positive_count=True)
        # Remove the final returned value printed by the CLI.
        self.assertIsNone(constant)
        self.assertEqual(records[0]['n'], '1')
        self.assertEqual(completed.stdout.splitlines()[-1],
                         blockirwin(2, '01', 1, 10))

    def test_decreasing_precision_and_diagnostics(self):
        output = io.StringIO()
        with redirect_stdout(output):
            result = modalirwin(2, '01', 0, 45, modalcoeffs=True,
                                verbose=True)
        self.assertEqual(result, blockirwin(2, '01', 0, 45))
        _, records = self.coefficient_records(output.getvalue(),
                                              positive_count=False)
        budgets = [int(record['budget digits']) for record in records]
        actual = [int(record['working digits']) for record in records]
        self.assertTrue(any(value < 57 for value in budgets))
        self.assertTrue(any(value < 57 for value in actual))
        self.assertTrue(all(value >= 12 for value in budgets))
        self.assertTrue(all(value >= budget for value, budget in zip(actual, budgets)))
        self.assertIn('cache rebuilds', output.getvalue())

    def test_raised_argument_signs_not_merged(self):
        from mpmath import mp
        from modalirwin import _ModalEngine
        with mp.workdps(50):
            engine = _ModalEngine(2, [1, 1], 3)
            raised = engine._raised_arguments(2, 1)
            self.assertTrue(all(pos >= 0 and neg >= 0 for pos, neg in raised.values()))
            self.assertTrue(any(pos > 0 and neg > 0 for pos, neg in raised.values()))
            term, _, positive, negative = engine.modal_term(2, 1, mp.mpf('1e-35'))
            self.assertTrue(mp.almosteq(term, positive-negative,
                                        rel_eps=mp.mpf('1e-30')))

    def test_early_stopping_and_verbose_layout(self):
        for digits in (20, 40):
            with self.subTest(digits=digits):
                output = io.StringIO()
                with redirect_stdout(output):
                    result = modalirwin(10, '314', 10, digits, verbose=True)
                self.assertEqual(result, blockirwin(10, '314', 10, digits))
                text = output.getvalue()
                modes = re.findall(r'(?m)^mode\s+= (\d+)$', text)
                self.assertEqual(modes, ['1', '2', '3'])
                self.assertRegex(text, r'(?m)^modes\s+= 3$')
                self.assertEqual(len(re.findall(r'(?m)^term\s+= ', text)), 3)
                self.assertEqual(len(re.findall(r'(?m)^\|term\|\s+= ', text)), 3)
                self.assertNotIn(';', text)

    def test_exact_binary_empty_case_does_not_construct_engine(self):
        with patch('modalirwin._ModalEngine', side_effect=AssertionError('engine called')):
            self.assertEqual(modalirwin(2, '1', 0, 20), '0')
            self.assertEqual(modalirwin(2, '1', 0, 20, modalcoeffs=True), '0')

    def test_initial_guard_and_relative_threshold(self):
        from mpmath import mp
        output = io.StringIO()
        with redirect_stdout(output):
            modalirwin(10, '314', 10, 20, verbose=True)
        fields = dict(re.findall(r'(?m)^([a-z][a-z -]+?)\s+= (\S+)$',
                                 output.getvalue()))
        self.assertEqual(fields['working decimal digits'], '32')
        self.assertEqual(fields['guard decimal digits'], '12')
        self.assertGreaterEqual(int(fields['successive decreases']), 2)
        with mp.workdps(50):
            expected = 10**3 * mp.log(10) * mp.power(10, -32)
            self.assertTrue(mp.almosteq(mp.mpf(fields['threshold']), expected,
                                         rel_eps=mp.mpf('1e-4')))

    def test_high_occurrence_count(self):
        self.check(2, '01', 6, 16)

    def test_strict_digit_input(self):
        # Nonintegral digits must never be truncated to integers.
        with self.assertRaises(TypeError):
            modalirwin(10, [1, 2.5], 0, 10)

    def test_invalid(self):
        for args in ((1, '1', 0, 10), (2,'2',0,10), (2,'1',-1,10),
                     (10, 12, 1, 10), (2, '', 0, 10)):
            with self.assertRaises((ValueError, TypeError)):
                modalirwin(*args)


if __name__ == '__main__':
    unittest.main()
