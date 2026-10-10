# -*- coding: utf-8 -*-
"""Block-Irwin and Kempner sums by J.-F. Burnol's modal expansions.

The mathematics is entirely due to Jean-François Burnol:
* "Block-count constrained harmonic sums: spectral expansion and
  block-directed Euler--Maclaurin" (2026), Part 2, for k >= 1;
* "The modal expansion of Kempner sums" (2026), for k = 0.

Independent Python + mpmath implementation produced by OpenAI GPT-6 on
October 10, 2026.  This module coexists with the separate implementation
in blockirwin.py, produced by Claude Opus 5.5 and Claude Sonnet 5.

Truncation and variable precision are heuristic, not certified bounds.
For usage, source references, tests and license see README-modalirwin.md.
Licensed under CC BY-SA 4.0; see the repository LICENSE file.
"""

from __future__ import annotations

from collections import defaultdict
from contextlib import redirect_stdout
from io import StringIO
from math import comb, ceil
from operator import index
import sys
import time
from mpmath import mp

__version__ = '0.2.0'


class _ModalNumericalInstability(ArithmeticError):
    """A coefficient recurrence became unreliable at this precision."""



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
            for i in range(index):
                c = comb(index,i)*self._mom(i+1,index-i)
                total_R += c*self.R[i]
                if self.need_kempner:
                    total_U += c*self.U[i]
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
        """Evaluate a finite linear combination of A_j(n) by Newton series.

        Raise each argument before evaluation.  The first omitted row
        with the tilted moment replaced by 1 estimates its magnitude;
        only the moment recurrences retain full engine precision.
        """
        raised = defaultdict(lambda: mp.zero)
        for n, coeff in arguments.items():
            if not coeff:
                continue
            for N, rcoeff in self._raised_arguments(j,n).items():
                raised[N] += coeff*rcoeff
        raised = {N: weight for N, weight in raised.items() if weight}
        if not raised:
            return mp.zero
        # Grouped Newton rows share the same M_{j;q}.
        total = mp.zero
        powers = {N: mp.power(N,-j) for N in raised}
        coeff_guard = max(0, int(mp.mag(max(abs(x) for x in raised.values()))*0.30103)+2)
        q = 0
        # Omitted-row estimate uses absolute argument weights; this avoids
        # cancellation but is not an estimate of the entire series tail.
        while True:
            binomial = comb(j+q-1,q)
            bound = binomial*mp.fsum(abs(raised[N])*abs(powers[N])
                                    for N in raised)
            # Obtain any new nu_j moments at the engine's *full* precision.
            moment = self._mom(j,q)
            # Remaining precision budget decreases with the absolute size
            # of the Newton row.  This is a heuristic, not a certification.
            vanished_digits = max(0, int(-mp.mag(bound)*0.30103)) if bound else 9999
            row_dps = min(self._base_dps, max(25, self._base_dps
                          - max(0, vanished_digits-coeff_guard-6)))
            with mp.workdps(row_dps):
                weights = mp.fsum(raised[N]*powers[N] for N in raised)
                next_term = (-1 if q%2 else 1)*binomial*weights*moment
            total += next_term
            self._stats['newton_rows'] += 1
            self._stats['reduced_rows'] += int(row_dps < self._base_dps)
            for N in powers:
                powers[N] /= N
            q += 1
            # Use ABSOLUTE weights (avoid cancellation of the a priori bound).
            omitted = comb(j+q-1,q)*mp.fsum(abs(raised[N])*abs(powers[N])
                                               for N in raised)
            if omitted < tolerance:
                break
            if q > 20000:
                raise ArithmeticError('Newton series did not reach requested tolerance')
        return total

    def modal_term(self, j, k, tolerance, modalcoeffs=False):
        m = j-1
        R, U = self._basis_coefficients(m)
        if j == 1:
            if k >= 1:
                left = mp.log(self.b)
                term = self.mass_v * left
                return (term, (mp.one, left, self.mass_v)) if modalcoeffs else term
            if all(x == 0 for x in self.w):
                term = (self.mass_u*mp.log(self.b) - self.mass_v*
                        mp.fsum(mp.log1p(mp.one/(a*self.b**self.p))
                                for a in range(1,self.b)))
            else:
                term = (self.mass_u*mp.log(self.b) -
                        (self.mass_v*mp.log1p(mp.one/self.nw)
                         if self.w[0] else mp.zero))
            return (term, None) if modalcoeffs else term
        if k >= 1:
            lam = self.lam(j)
            tolerance = tolerance/max(mp.mpf('1e-60'),abs(R*lam**(k-1)))
            a = defaultdict(lambda: mp.zero)
            for d in range(1,self.b):
                a[d] += lam
            for h,suffix in self.right:
                for d in range(1,self.b):
                    a[d*self.b**h+suffix] += lam-1
            if self.w[0]:
                a[self.nw] -= lam-1
            left = (-1 if m%2 else 1)*self.linear_A(j,a,tolerance)
            term = left*lam**(k-1)*R
            return (term, (lam, left, R)) if modalcoeffs else term
        a = defaultdict(lambda: mp.zero)
        if all(x == 0 for x in self.w):
            # Level-2 Kempner series: finite initial harmonic sum external.
            for n in range(self.b,self.b*self.b):
                a[n] += U
            for d in range(1,self.b):
                a[d*self.b**self.p] -= R
        else:
            for d in range(1,self.b):
                a[d] += U
            if self.w[0]:
                a[self.nw] -= R
        term = (-1 if m%2 else 1)*self.linear_A(j,a,tolerance)
        return (term, None) if modalcoeffs else term


def _modalirwin_once(b: int, w, k: int, nboffigures: int = 100,
                     ell=None, verbose=False, trunc=False, *, max_modes=1000,
                     modalcoeffs=False, _extra_guard=0) -> str:
    """Run one modal calculation at a fixed full recurrence precision.

    This worker is private.  The public entry point handles safe restarts
    if the numerical coefficient-growth safeguard is triggered.
    """
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
    # Extra digits protect moment/eigenpolynomial cancellation.  The public
    # wrapper can retry this worker with a larger guard if it is insufficient.
    guard = max(45, ceil(0.8*nboffigures)+30,
                (len(digits)-1)*nboffigures+30 if b == 2 else 0)
    full_dps = nboffigures + guard + _extra_guard
    # Three additional target digits provide a pragmatic stopping margin.
    target = mp.power(10, -(nboffigures+3))
    engine = None
    terms = []
    small_count = 0
    with mp.workdps(full_dps):
        engine = _ModalEngine(b,digits,ell,need_kempner=(k == 0))
        engine._base_dps = full_dps
        answer = mp.zero
        if k == 0 and all(x == 0 for x in digits):
            constant = mp.fsum(mp.one/a for a in range(1,b))
            answer += constant
            if modalcoeffs:
                print(f'constant = {mp.nstr(constant, nboffigures)}')
        # Full precision for cached recurrences, reduced precision for small
        # Newton rows.  The outer modal stopping estimate is empirical.
        for j in range(1,max_modes+1):
            # Per-mode A_j accuracy is scaled to the modal contribution.
            # This does not affect the precision of stored M, u, v or R.
            mode_result = engine.modal_term(j, k, target/(100 if j == 1 else 10),
                                            modalcoeffs=modalcoeffs)
            if modalcoeffs:
                term, factors = mode_result
            else:
                term = mode_result
            if (j >= 10 and terms and abs(term) > 0 and
                    mp.mag(abs(term)) >
                    max((mp.mag(abs(x)) for x in terms[-5:] if x),
                        default=-mp.inf) + 20):
                raise _ModalNumericalInstability(
                    'modal coefficient recurrence lost numerical stability')
            answer += term
            terms.append(term)
            # Proxy for the first omitted term: twice the most recent
            # magnitude.  Require three consecutive small modes to stop;
            # neither this proxy nor the tail is rigorously bounded.
            first_omitted_estimate = 2*abs(term)
            if j >= 5 and first_omitted_estimate < target:
                small_count += 1
            else:
                small_count = 0
            if modalcoeffs:
                signed_term = ('+' if term >= 0 else '') + mp.nstr(term, nboffigures)
                if k >= 1:
                    lam, left, right = factors
                    product = left*right
                    print(f'n={j:4d}  lambda={mp.nstr(lam, nboffigures)}  '
                          f'L={mp.nstr(left, nboffigures)}  '
                          f'R={mp.nstr(right, nboffigures)}  '
                          f'C={mp.nstr(product, nboffigures)}  '
                          f"C'={mp.nstr(product/lam, nboffigures)}  "
                          f'T={signed_term}')
                else:
                    print(f'n={j:4d}  T={signed_term}')
            elif verbose:
                print(f'mode {j:4d}: term {mp.nstr(term, 12)}, '
                      f'|term|={mp.nstr(abs(term), 5)}')
            if small_count >= 3:
                break
        else:
            raise ArithmeticError(f'modal expansion exceeded max_modes={max_modes}')
        result = mp.nstr(answer, nboffigures+7 if trunc else nboffigures,
                        strip_zeros=False, min_fixed=-mp.inf, max_fixed=mp.inf)
        if trunc:
            result = result[:-7]
        if verbose:
            print(f'{len(terms)} modes, last = {mp.nstr(terms[-1],6)}; '
                  f'{engine._stats["nu_moments"]} tilted moments, '
                  f'{engine._stats["basic_moments"]} base moments, '
                  f'{len(engine._coefficient_rows)} cached coefficient rows, '
                  f'{engine._stats["newton_rows"]} Newton rows '
                  f'({engine._stats["reduced_rows"]} reduced precision); '
                  f'{time.perf_counter()-started:.3f}s')
            print(f'First-omitted modal term estimate: '
                  f'{mp.nstr(first_omitted_estimate, 5)}; '
                  f'threshold {mp.nstr(target, 5)} '
                  '(heuristic; NOT a certified error bound)')
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
    ``max_modes`` limits runtime; reaching it raises ArithmeticError.

    Precision reductions, modal stopping and automatic stability retries
    are numerical heuristics, not rigorous accuracy certificates.
    """
    # Most calls succeed without a restart.  A new worker must be created
    # on each retry: cached recurrence values cannot be promoted in place.
    # Redirect output until the successful attempt, so a failed attempt
    # cannot print a misleading partial coefficient list.
    retries = (0, 80, 240)
    for attempt, extra_guard in enumerate(retries):
        output = StringIO()
        try:
            with redirect_stdout(output):
                result = _modalirwin_once(
                    b, w, k, nboffigures, ell, verbose, trunc,
                    max_modes=max_modes, modalcoeffs=modalcoeffs,
                    _extra_guard=extra_guard)
        except _ModalNumericalInstability as exc:
            if attempt + 1 == len(retries):
                raise ArithmeticError(
                    'modal calculation remained numerically unstable '
                    'after increasing working precision') from exc
            continue
        if verbose and attempt:
            print(f'Precision safeguard: restarted with {extra_guard} '
                  'additional guard digits.')
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
