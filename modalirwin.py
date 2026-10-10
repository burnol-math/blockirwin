# -*- coding: utf-8 -*-
"""Block-Irwin and Kempner sums by J.-F. Burnol's modal expansions.

The mathematics is entirely due to Jean-François Burnol:
* "Block-count constrained harmonic sums: spectral expansion and
  block-directed Euler--Maclaurin" (2026), Part 2, for k >= 1;
* "The modal expansion of Kempner sums" (2026), for k = 0.

Independent Python + mpmath implementation produced by OpenAI GPT-6 on
October 10, 2026.  This module coexists with the separate implementation
in blockirwin.py, produced by Claude Opus 5.5 and Claude Sonnet 5.

Modal stopping and cancellation checks are numerical heuristics, not certified bounds.
For usage, source references, tests and license see README-modalirwin.md.
Licensed under CC BY-SA 4.0; see the repository LICENSE file.
"""

from __future__ import annotations

from collections import defaultdict
from contextlib import redirect_stdout
from io import StringIO
from math import comb
from operator import index
import sys
import time
from mpmath import mp

__version__ = '0.3.0'


class _ModalNumericalInstability(ArithmeticError):
    """A modal term needs more working precision."""


def _print_fields(fields):
    """Print one labeled diagnostic value per line, aligned by label width."""
    width = max(len(label) for label, _ in fields)
    for label, value in fields:
        print(f'{label:<{width}} = {value}')


def _digits(base, word):
    if isinstance(word, int):
        raise TypeError('w must be a string or digit sequence (leading zeros matter)')
    if isinstance(word, str):
        if base > 36:
            raise ValueError('for bases > 36, use a sequence of digits')
        digits = [int(char, 36) for char in word.lower()]
    else:
        # Reject nonintegral entries instead of silently truncating floats.
        try:
            digits = [index(x) for x in word]
        except TypeError as exc:
            raise TypeError('digit sequences must contain integers') from exc
    if not digits or any(x < 0 or x >= base for x in digits):
        raise ValueError('w must be a nonempty word of base-b digits')
    return digits


def _word_value(b, seq):
    n = 0
    for d in seq:
        n = b*n + d
    return n


class _ModalEngine:
    """Triangular recurrences for the modal moments and coefficients.

    Both base moment recurrences have diagonal D(b**(-m-1)); their
    off-diagonal coefficients use different word-overlap orientations.
    All caches belong to one engine at one working precision.
    """

    def __init__(self, b, word, ell, need_kempner=False):
        self.b, self.w, self.p, self.ell = b, word, len(word), ell
        self.need_kempner = need_kempner
        self.nw = _word_value(b, word)
        self.overlaps = [h for h in range(1, self.p)
                         if word[:self.p-h] == word[h:]]
        self.right = [(h, _word_value(b, word[self.p-h:])) for h in self.overlaps]
        self.left = [(h, _word_value(b, word[:h])) for h in self.overlaps]
        self.C1 = mp.fsum([mp.one] + [mp.power(b, -h) for h in self.overlaps])
        self.mass_u = b**self.p * self.C1
        self.mass_v = mp.mpf(b**self.p)
        # Signed finite shift lists for the eta_0 and mu_0 recurrences.
        self._vterms = self._terms(self.right)
        self._uterms = self._terms(self.left, left=True)
        self.v = [self.mass_v]
        self.u = [self.mass_u]
        self.M = {}        # M[j] contains all computed moments of nu_j
        self.R = [self.mass_v]  # Right coefficients: integrals of eigenpolynomials against eta_0.
        self.U = [self.mass_u]  # Kempner counterparts against mu_0.
        self._basis_gross = [(abs(self.mass_v), abs(self.mass_u))]
        self._lams = {}
        self._diagonal = {}
        self._coefficient_rows = {}
        self._stats = defaultdict(int)

    def _terms(self, overlaps, left=False):
        base = [(1, a, -1) for a in range(self.b)]
        tails = []
        for h, suffix in overlaps:
            tails.append((h, suffix, 1))
            for a in range(self.b):
                tails.append((h+1, suffix*self.b+a if left else a*self.b**h+suffix, -1))
        tails.append((self.p, self.nw, 1))
        return base + tails

    def _de(self, row, col, which='v'):
        """Return a recurrence coefficient, caching each row's shift powers."""
        degree = row-col
        key = (which, row)
        cached = self._coefficient_rows.get(key)
        if cached is None:
            t = mp.power(self.b, -row-1)
            terms = self._vterms if which == 'v' else self._uterms
            # The evaluation point depends on the row, not on the degree.
            # Cache its weights, then extend the degree table on demand.
            weighted = [(integer, sign*t**length) for length, integer, sign in terms]
            initial_e = mp.fsum(weight for _, weight in weighted[self.b:])
            initial_d = mp.one + initial_e + mp.fsum(weight for _, weight in weighted[:self.b])
            cached = ([initial_d], [initial_e], weighted)
            self._coefficient_rows[key] = cached
        d_values, e_values, weighted = cached
        while len(d_values) <= degree:
            d = len(d_values)
            a_part = mp.zero
            e_part = mp.zero
            for index, (integer, weight) in enumerate(weighted):
                if integer == 0:
                    continue
                contribution = weight * integer**d
                if index < self.b:
                    a_part += contribution
                else:
                    e_part += contribution
            d_values.append(a_part + e_part)
            e_values.append(e_part)
        return d_values[degree], e_values[degree]

    def _D(self, exponent):
        cached = self._diagonal.get(exponent)
        if cached is None:
            t = mp.power(self.b, -exponent)
            C = mp.one + mp.fsum(t**h for h in self.overlaps)
            cached = (1-self.b*t)*C+t**self.p
            self._diagonal[exponent] = cached
        return cached

    def lam(self, j):
        if j not in self._lams:
            t = mp.power(self.b, -j)
            C = mp.one+mp.fsum(t**h for h in self.overlaps)
            D = (1-self.b*t)*C+t**self.p
            self._lams[j] = ((1-self.b*t)*(C-1)+t**self.p)/D
        return self._lams[j]

    def _extend_basic(self, order):
        """Extend eta_0 moments and, only for k=0, autonomous mu_0 moments."""
        for m in range(len(self.v), order+1):
            d0 = self._D(m+1)
            vtotal = mp.zero
            utotal = mp.zero
            for i in range(m):
                factor = comb(m,i)
                dv, _ = self._de(m,i,'v')
                vtotal += factor*dv*self.v[i]
                if self.need_kempner:
                    du, _ = self._de(m,i,'u')
                    utotal += factor*du*self.u[i]
            self.v.append(-vtotal/d0)
            if self.need_kempner:
                forcing = mp.fsum(mp.power(self.b, -(m+1)*h)*(value**m)
                                  for h,value in self.left)
                self.u.append((forcing-utotal)/d0)
            self._stats['basic_moments'] += 1

    def _mom(self, j, q):
        """Moment of the tilted probability measure nu_j (Part 2)."""
        moments = self.M.setdefault(j, [mp.one])
        lamj = self.lam(j)
        while len(moments) <= q:
            r = len(moments)
            row = j+r-1
            numerator = mp.zero
            for i in range(r):
                col = j+i-1
                d, e = self._de(row,col)
                numerator += comb(r,i)*(e-lamj*d)*moments[i]
            denominator = self._D(j+r)*(lamj-self.lam(j+r))
            moments.append(numerator/denominator)
            self._stats['nu_moments'] += 1
        return moments[q]

    def _basis_coefficients(self, m):
        """Obtain modal right coefficients from monomial base moments."""
        self._extend_basic(m)
        while len(self.R) <= m:
            index = len(self.R)
            total_R = mp.zero
            total_U = mp.zero
            gross_R = abs(self.v[index])
            gross_U = abs(self.u[index]) if self.need_kempner else mp.zero
            for i in range(index):
                c = comb(index,i)*self._mom(i+1,index-i)
                summand_R = c*self.R[i]
                total_R += summand_R
                gross_R += abs(summand_R)
                if self.need_kempner:
                    summand_U = c*self.U[i]
                    total_U += summand_U
                    gross_U += abs(summand_U)
            # Record gross coefficient magnitudes for cancellation checks.
            self._basis_gross.append((gross_R, gross_U))
            self.R.append(self.v[index]-total_R)
            if self.need_kempner:
                self.U.append(self.u[index]-total_U)
        return self.R[m], self.U[m] if self.need_kempner else None

    def _raised_arguments(self, j, initial):
        """Raise auxiliary A_j arguments to at least b**(ell-1).

        Applies the finite identity (Ajnlevelraising) from Part 2,
        merging repeated arguments and preserving signed weights.
        """
        counts = defaultdict(lambda: mp.zero)
        counts[int(initial)] = mp.one
        cutoff = self.b**(self.ell-1)
        lam = self.lam(j)
        c = (1-lam)/lam
        while True:
            small = [(n, coeff) for n, coeff in counts.items()
                     if n < cutoff and coeff]
            if not small:
                break
            n = min(n for n, _ in small)
            coeff = counts.pop(n)
            for a in range(self.b):
                counts[self.b*n+a] += coeff
            for h, suffix in self.right:
                counts[self.b**h*n+suffix] += coeff*c
                for a in range(self.b):
                    counts[self.b**(h+1)*n+a*self.b**h+suffix] -= coeff*c
            counts[self.b**self.p*n+self.nw] += coeff*c
        return {n: coeff for n, coeff in counts.items() if coeff}

    def linear_A(self, j, arguments, tolerance):
        """Evaluate an A_j combination and its positive/negative totals.

        All signed Newton contributions are classified *before* their
        cancellation.  This permits a per-modal-term cancellation check.
        The omitted-row tolerance remains heuristic.
        """
        raised_plus = defaultdict(lambda: mp.zero)
        raised_minus = defaultdict(lambda: mp.zero)
        for n, coeff in arguments.items():
            if not coeff:
                continue
            for N, rcoeff in self._raised_arguments(j,n).items():
                contribution = coeff*rcoeff
                if contribution >= 0:
                    raised_plus[N] += contribution
                else:
                    raised_minus[N] -= contribution
        # Do not cancel these contributions before measuring their sizes.
        # The net coefficients are needed only for the Newton cutoff test.
        raised = {N: raised_plus[N]-raised_minus[N]
                  for N in raised_plus.keys() | raised_minus.keys()}
        if not raised:
            return mp.zero, mp.zero, mp.zero
        powers = {N: mp.power(N,-j) for N in raised}
        positives, negatives = [], []
        q = 0
        while True:
            binomial = comb(j+q-1,q)
            moment = self._mom(j,q)
            signed_factor = (-1 if q%2 else 1)*binomial*moment
            # Evaluate signed contributions at the engine precision.
            for N in raised:
                positive_part = signed_factor*raised_plus[N]*powers[N]
                negative_part = -signed_factor*raised_minus[N]*powers[N]
                for contribution in (positive_part, negative_part):
                    if contribution >= 0:
                        positives.append(contribution)
                    else:
                        negatives.append(-contribution)
            self._stats['newton_rows'] += 1
            for N in powers:
                powers[N] /= N
            q += 1
            omitted = comb(j+q-1,q)*mp.fsum(abs(raised[N])*abs(powers[N])
                                                for N in raised)
            if omitted < tolerance:
                break
            if q > 20000:
                raise ArithmeticError('Newton series did not reach requested tolerance')
        positive = mp.fsum(positives)
        negative = mp.fsum(negatives)
        return positive-negative, positive, negative

    def modal_term(self, j, k, tolerance, modalcoeffs=False):
        """Return one mode and the positive/negative contributions to it.

        For k>=1 the right modal factor scales the signed Newton sums.
        For k=0 the left/right coefficients enter the Newton combination.
        """
        m = j-1
        R, U = self._basis_coefficients(m)
        if j == 1:
            if k >= 1:
                left = mp.log(self.b)
                term = self.mass_v * left
                return term, (mp.one, left, self.mass_v), abs(term), mp.zero
            if all(x == 0 for x in self.w):
                positive = self.mass_u*mp.log(self.b)
                negative = self.mass_v*mp.fsum(
                    mp.log1p(mp.one/(a*self.b**self.p))
                    for a in range(1,self.b))
            else:
                positive = self.mass_u*mp.log(self.b)
                negative = (self.mass_v*mp.log1p(mp.one/self.nw)
                            if self.w[0] else mp.zero)
            return positive-negative, None, positive, negative
        if k >= 1:
            lam = self.lam(j)
            multiplier = R*lam**(k-1)
            tolerance = tolerance/max(mp.mpf('1e-60'), abs(multiplier))
            a = defaultdict(lambda: mp.zero)
            for d in range(1,self.b):
                a[d] += lam
            for h,suffix in self.right:
                for d in range(1,self.b):
                    a[d*self.b**h+suffix] += lam-1
            if self.w[0]:
                a[self.nw] -= lam-1
            left, positive, negative = self.linear_A(j,a,tolerance)
            if m%2:
                left, positive, negative = -left, negative, positive
            term = left*multiplier
            if multiplier < 0:
                positive, negative = negative, positive
            magnitude = abs(multiplier)
            return term, (lam, left, R), positive*magnitude, negative*magnitude
        a = defaultdict(lambda: mp.zero)
        if all(x == 0 for x in self.w):
            for n in range(self.b,self.b*self.b):
                a[n] += U
            for d in range(1,self.b):
                a[d*self.b**self.p] -= R
        else:
            for d in range(1,self.b):
                a[d] += U
            if self.w[0]:
                a[self.nw] -= R
        term, positive, negative = self.linear_A(j,a,tolerance)
        if m%2:
            term, positive, negative = -term, negative, positive
        return term, None, positive, negative


def _modalirwin_once(b: int, w, k: int, nboffigures: int = 100,
                     ell=None, verbose=False, trunc=False, *, max_modes=1000,
                     modalcoeffs=False, _guard=12) -> str:
    """Run one calculation with the given decimal guard precision."""
    if int(b) != b or b < 2:
        raise ValueError('b must be an integer >= 2')
    b = int(b)
    digits = _digits(b,w)
    if int(k) != k or k < 0:
        raise ValueError('k must be a nonnegative integer')
    k = int(k)
    if int(nboffigures) != nboffigures or nboffigures < 6:
        raise ValueError('nboffigures must be >= 6')
    nboffigures = int(nboffigures)
    if ell is None:
        ell = 3 if b == 2 else 2
    if int(ell) != ell or ell < 2:
        raise ValueError('ell must be an integer >= 2')
    ell = int(ell)
    if int(max_modes) != max_modes or max_modes < 2:
        raise ValueError('max_modes must be an integer >= 2')
    if k == 0 and b == 2 and digits == [1]:
        if modalcoeffs:
            print('No admissible positive integers; no modal terms (I=0).')
        return '0'

    started = time.perf_counter()
    # Each attempt uses its supplied decimal guard precision.
    full_dps = nboffigures + _guard
    terms = []
    decreasing_run = 0
    cancellation_max = 0
    with mp.workdps(full_dps):
        # A modal term is undetectable relative to the natural scale
        # b**len(w) * log(b) at requested figures plus twelve digits.
        scale = mp.power(b, len(digits)) * mp.log(b)
        target = scale * mp.power(10, -(nboffigures+12))
        engine = _ModalEngine(b,digits,ell,need_kempner=(k == 0))
        answer = mp.zero
        if k == 0 and all(x == 0 for x in digits):
            constant = mp.fsum(mp.one/a for a in range(1,b))
            answer += constant
            if modalcoeffs:
                _print_fields([('constant', mp.nstr(constant, nboffigures))])
                print()
        # Recurrences and Newton contributions use one working precision.
        for j in range(1,max_modes+1):
            # Scale the Newton-series cutoff to this modal contribution.
            term, factors, positive, negative = engine.modal_term(
                j, k, target/(100 if j == 1 else 10),
                modalcoeffs=modalcoeffs)
            magnitude = abs(term)
            # Estimate lost decimal digits from positive and negative
            # components of THIS mode only, not from other modal terms.
            gross = positive+negative
            # The right/left factors are themselves built by subtractive
            # recurrences.  Include cancellation in these factors when
            # estimating the precision needed for this one modal term.
            if j > 1:
                R, U = engine.R[j-1], (engine.U[j-1] if k == 0 else None)
                gross_R, gross_U = engine._basis_gross[j-1]
                if k >= 1:
                    if R:
                        gross *= max(mp.one, gross_R/abs(R))
                    elif gross_R:
                        gross = mp.inf
                else:
                    coeff_amplification = mp.one
                    if U and gross_U:
                        coeff_amplification = max(coeff_amplification,
                                                  gross_U/abs(U))
                    elif gross_U:
                        gross = mp.inf
                    if R and gross_R:
                        coeff_amplification = max(coeff_amplification,
                                                  gross_R/abs(R))
                    elif gross_R:
                        gross = mp.inf
                    gross *= coeff_amplification
            loss = (max(mp.zero, mp.log10(gross/magnitude))
                    if magnitude and gross else
                    (mp.inf if gross else mp.zero))
            cancellation_max = max(cancellation_max, loss)
            # With d working decimal digits, unresolved rounding is on
            # the scale of gross*10**(-d), even if term rounds to zero.
            # This is a numerical heuristic, not an interval error bound.
            rounding_proxy = gross*mp.power(10, -full_dps)
            if (gross and rounding_proxy >= target and
                    (magnitude < target or loss > _guard-2)):
                raise _ModalNumericalInstability(
                    f'mode {j}: positive/negative cancellation '
                    f'requires more than {_guard} guard digits')
            if terms:
                decreasing_run = (decreasing_run+1 if
                                  magnitude < abs(terms[-1]) else 0)
                # An isolated upward movement is not by itself evidence of
                # instability; repeated nondecrease without meeting the
                # stopping rule is handled by a bounded mode limit.
            # For an unresolved cancellation smaller than the threshold,
            # rounding_proxy < target is essential before treating it as
            # negligible.  Do not mistake an inaccurate zero for evidence.
            answer += term
            terms.append(term)
            # Heuristic proxy for the entire omitted tail: twice the
            # most recent term's magnitude.  Not a rigorous error bound.
            omitted_tail_estimate = 2 * abs(term)
            if modalcoeffs:
                if j > 1:
                    print()
                signed_term = ('+' if term >= 0 else '') + mp.nstr(term, nboffigures)
                fields = [('n', j)]
                if k >= 1:
                    lam, left, right = factors
                    product = left*right
                    fields.extend([
                        ('lambda', mp.nstr(lam, nboffigures)),
                        ('L', mp.nstr(left, nboffigures)),
                        ('R', mp.nstr(right, nboffigures)),
                        ('C', mp.nstr(product, nboffigures)),
                        ("C'", mp.nstr(product/lam, nboffigures)),
                    ])
                fields.append(('T', signed_term))
                _print_fields(fields)
            elif verbose:
                if j > 1:
                    print()
                _print_fields([
                    ('mode', j),
                    ('term', mp.nstr(term, 12)),
                    ('|term|', mp.nstr(abs(term), 5)),
                ])
            if (j >= 3 and decreasing_run >= 2 and
                    2*(magnitude+rounding_proxy) < target):
                break
        else:
            raise ArithmeticError(f'modal expansion exceeded max_modes={max_modes}')
        result = mp.nstr(answer, nboffigures+7 if trunc else nboffigures,
                        strip_zeros=False, min_fixed=-mp.inf, max_fixed=mp.inf)
        if trunc:
            result = result[:-7]
        if verbose:
            print('\nModal summary:')
            _print_fields([
                ('modes', len(terms)),
                ('last modal term', mp.nstr(terms[-1], 6)),
                ('tilted moments', engine._stats['nu_moments']),
                ('base moments', engine._stats['basic_moments']),
                ('cached coefficient rows', len(engine._coefficient_rows)),
                ('Newton rows', engine._stats['newton_rows']),
                ('maximum modal cancellation digits', mp.nstr(cancellation_max, 5)),
                ('elapsed time', f'{time.perf_counter()-started:.3f}s'),
                ('estimated omitted-tail magnitude',
                 mp.nstr(omitted_tail_estimate, 5)),
                ('working decimal digits', full_dps),
                ('guard decimal digits', _guard),
                ('successive decreases', decreasing_run),
                ('reference scale', mp.nstr(scale, 6)),
                ('threshold', mp.nstr(target, 5)),
            ])
            print('Stopping criterion: heuristic (NOT a certified error bound)')
    return result



def modalirwin(b: int, w, k: int, nboffigures: int = 100,
               ell=None, verbose=False, trunc=False, *, max_modes=1000,
               modalcoeffs=False) -> str:
    """Compute I(b,w,k) using the published modal expansion.

    The positional interface matches ``blockirwin``; return a decimal
    string of ``nboffigures`` significant figures.  ``ell`` controls
    level raising of auxiliary modal Stieltjes arguments, not the level
    setting of the other method.  ``modalcoeffs`` prints signed mode
    terms and, for k >= 1, lambda, L, R, C and C' (mode n starts at 1).
    ``max_modes`` limits the number of terms; reaching it raises
    ArithmeticError.

    Stop after three consecutive term magnitudes in strict descent and once
    the tail proxy is below b**len(w)*log(b)*10**(-nboffigures-12).
    Each term has separate positive and negative Newton contributions;
    substantial unresolved cancellation triggers fresh-cache retries at
    12, 24, 36, ... decimal guard digits. All error indicators and the
    modal tail criterion remain heuristic, not certified bounds.
    """
    # Each attempt owns an independent recurrence cache. Suppress output from
    # unsuccessful attempts, so displayed modes always belong to one run.
    guard = 12
    max_guard = 240
    while guard <= max_guard:
        output = StringIO()
        try:
            with redirect_stdout(output):
                result = _modalirwin_once(
                    b, w, k, nboffigures, ell, verbose, trunc,
                    max_modes=max_modes, modalcoeffs=modalcoeffs, _guard=guard)
        except _ModalNumericalInstability as exc:
            if guard == max_guard:
                raise ArithmeticError(
                    f'modal calculation remained unstable at '
                    f'{max_guard} guard digits') from exc
            guard += 12
            continue
        if verbose and guard > 12:
            print('Precision safeguard: restarted calculation')
            _print_fields([('guard decimal digits', guard)])
        sys.stdout.write(output.getvalue())
        return result
    raise AssertionError('unreachable precision-retry state')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('b', type=int, help='integer base (at least 2)')
    parser.add_argument('w', help='digit word; comma-separated digits for b>36')
    parser.add_argument('k', type=int, help='number of word occurrences')
    parser.add_argument('nboffigures', type=int, nargs='?', default=100)
    parser.add_argument('ell', type=int, nargs='?', default=None)
    parser.add_argument('--verbose', action='store_true',
                        help='print modal diagnostics and timing')
    parser.add_argument('--modalcoeffs', action='store_true',
                        help='print signed modal terms and coefficient factors')
    parser.add_argument('--trunc', action='store_true',
                        help='truncate the decimal result instead of rounding')
    parser.add_argument('--max-modes', type=int, default=1000,
                        help='maximum number of modal terms')
    ns = parser.parse_args()
    word = ns.w if ns.b <= 36 else [int(s) for s in ns.w.split(',')]
    print(modalirwin(ns.b, word, ns.k, ns.nboffigures, ns.ell,
                     ns.verbose, ns.trunc, max_modes=ns.max_modes,
                     modalcoeffs=ns.modalcoeffs))
