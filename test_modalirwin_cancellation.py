"""Focused regression checks for mode-level precision management.

These tests inspect signed contributions, the absence of wrapper-level
restarts, stopping diagnostics, and observed precision reduction. They
do not certify a mathematical error bound for every computed mode.
"""
import io
import re
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch
from mpmath import mp

import modalirwin as module
from blockirwin import blockirwin


class CancellationPolicyTests(unittest.TestCase):
    def test_no_whole_calculation_guard_restart(self):
        attempts = []
        def fake_worker(*args, **kwargs):
            attempts.append(kwargs.get('_guard', 12))
            return 'complete'
        with patch.object(module, '_modalirwin_once', side_effect=fake_worker):
            self.assertEqual(module.modalirwin(10, '314', 10, 20), 'complete')
        self.assertEqual(attempts, [12])

    def test_separate_positive_and_negative_contributions(self):
        for b, word, k, j in ((10, [3,1,4], 10, 2), (2, [0,1], 0, 2)):
            with self.subTest(b=b, k=k), mp.workdps(55):
                engine = module._ModalEngine(b, word, 3 if b == 2 else 2,
                                              need_kempner=(k == 0))
                engine._base_dps = 55
                term, factors, positive, negative = engine.modal_term(
                    j, k, mp.mpf('1e-38'), modalcoeffs=True)
                self.assertGreaterEqual(positive, 0)
                self.assertGreaterEqual(negative, 0)
                self.assertTrue(mp.almosteq(term, positive-negative,
                                            rel_eps=mp.mpf('1e-35'),
                                            abs_eps=mp.mpf('1e-70')))

    def test_three_mode_stop_and_reference(self):
        for digits in (20, 40):
            with self.subTest(digits=digits):
                out = io.StringIO()
                with redirect_stdout(out):
                    value = module.modalirwin(10, '314', 10, digits,
                                              verbose=True, max_modes=10)
                self.assertEqual(value, blockirwin(10, '314', 10, digits))
                self.assertRegex(out.getvalue(), r'(?m)^modes\s+= 3$')
                self.assertRegex(out.getvalue(), r'(?m)^guard decimal digits\s+= 12$')
                self.assertIn('maximum modal cancellation digits', out.getvalue())

    def test_later_modes_use_less_precision_at_40_digits(self):
        """Check that the 40-digit example uses fewer than 52 digits."""
        report = io.StringIO()
        with redirect_stdout(report):
            result = module.modalirwin(2, '01', 0, 40, verbose=True)
        self.assertEqual(result, blockirwin(2, '01', 0, 40))
        self.assertRegex(report.getvalue(),
                         r'(?m)^working decimal digits\s+= 52$')
        found = re.search(r'(?m)^minimum modal decimal digits\s+= (\d+)$',
                          report.getvalue())
        self.assertIsNotNone(found)
        self.assertLess(int(found.group(1)), 52)

    def test_binary_empty_sum_is_exact_without_engine(self):
        with patch.object(module, '_ModalEngine', side_effect=AssertionError):
            self.assertEqual(module.modalirwin(2, '1', 0, 12), '0')


if __name__ == '__main__':
    unittest.main()
