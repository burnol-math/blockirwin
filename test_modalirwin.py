# -*- coding: utf-8 -*-
"""Regression checks for Burnol's modal expansions, implemented by GPT-6.

Produced October 10, 2026.  The formulas are entirely J.-F. Burnol's;
see README-modalirwin.md and the repository CC BY-SA 4.0 LICENSE.
"""
import io
from pathlib import Path
import re
import subprocess
import sys
import unittest
from contextlib import redirect_stdout

from blockirwin import blockirwin
from modalirwin import modalirwin
from compare_modalirwin import compare_irwin


class TestModalExpansion(unittest.TestCase):
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
        lines = output.getvalue().splitlines()
        self.assertTrue(lines)
        terms = []
        for number, line in enumerate(lines, start=1):
            self.assertRegex(line, rf'^n=\s*{number}\s+lambda=')
            values = dict(re.findall(r"(lambda|L|R|C'|C|T)=([^ ]+)", line))
            self.assertEqual(set(values), {'lambda','L','R','C',"C'",'T'})
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
                lines = output.getvalue().splitlines()
                if has_constant:
                    self.assertEqual(lines.pop(0), 'constant = 1.0')
                else:
                    self.assertNotIn('constant', output.getvalue())
                for number, line in enumerate(lines, start=1):
                    self.assertRegex(line, rf'^n=\s*{number}\s+T=[+-]')
                    self.assertNotIn('lambda=', line)
                    self.assertNotIn('L=', line)
                    self.assertNotIn('R=', line)
                with mp.workdps(40):
                    total = (mp.one if has_constant else mp.zero) + mp.fsum(
                        mp.mpf(line.split('T=')[1]) for line in lines)
                    self.assertLess(abs(total - mp.mpf(result)), mp.mpf('1e-12'))

    def test_modalcoeffs_verbose_and_special_empty_sum(self):
        output = io.StringIO()
        with redirect_stdout(output):
            result = modalirwin(10, '42', 3, 10, verbose=True, modalcoeffs=True)
        self.assertEqual(result, blockirwin(10, '42', 3, 10))
        self.assertIn('First-omitted modal term estimate', output.getvalue())
        self.assertIn("C'=", output.getvalue())
        self.assertNotIn('mode    1: term', output.getvalue())
        output = io.StringIO()
        with redirect_stdout(output):
            self.assertEqual(modalirwin(2, '1', 0, 10, modalcoeffs=True), '0')
        self.assertIn('no modal terms', output.getvalue())

    def test_modalcoeffs_cli(self):
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).with_name('modalirwin.py')),
             '2', '01', '1', '10', '--modalcoeffs'],
            capture_output=True, text=True, check=True)
        self.assertRegex(completed.stdout, r'^n=\s*1\s+lambda=')
        self.assertEqual(completed.stdout.splitlines()[-1],
                         blockirwin(2, '01', 1, 10))

    def test_numerical_stability_retry_and_clean_output(self):
        # This slowly convergent bordered ternary Kempner sum previously
        # exhausted the default guard digits.  A failed attempt must not
        # leak a partial, duplicated list of modal coefficients.
        output = io.StringIO()
        with redirect_stdout(output):
            result = modalirwin(3, '101', 0, 10,
                                modalcoeffs=True, verbose=True)
        self.assertEqual(result, blockirwin(3, '101', 0, 10))
        text = output.getvalue()
        self.assertIn('Precision safeguard: restarted', text)
        mode_numbers = [int(num) for num in
                        re.findall(r'^n=\s*(\d+)\s+T=', text, re.M)]
        self.assertGreater(len(mode_numbers), 10)
        self.assertEqual(mode_numbers, list(range(1, len(mode_numbers) + 1)))

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
