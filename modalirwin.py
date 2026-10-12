# -*- coding: utf-8 -*-
"""Block-Irwin and Kempner sums by J.-F. Burnol's modal expansions.

The mathematics is entirely due to Jean-François Burnol:
* "Block-count constrained harmonic sums: spectral expansion and
  block-directed Euler--Maclaurin" (2026), Part 2, for k >= 1;
* "The modal expansion of Kempner sums" (2026), for k = 0.

Independent Python + mpmath implementation produced by OpenAI GPT-6 on
October 10, 2026.  Level raising and the Newton loop rewritten by
Claude Opus 5.5 on October 12, 2026.  This module coexists with the
separate implementation in blockirwin.py, produced by Claude Opus 5.5
and Claude Sonnet 5.

The stopping test uses an estimated tail bound, not a certified bound.
For usage, source references, tests and license see README-modalirwin.md.
Licensed under CC BY-SA 4.0; see the repository LICENSE file.
"""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from heapq import heappush, heappop
from math import comb, ceil
from operator import index
import time
from mpmath import mp

__version__ = '0.6.0'


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


def _horner(coefficients, x):
    """Value at x of the polynomial with the given nonnegative integer
    coefficients (constant term first); zero for an empty tuple."""
    value = mp.zero
    for coefficient in reversed(coefficients):
        value = value*x + coefficient
    return value


@lru_cache(maxsize=64)
def _raised_table(b, w, ell, initial):
    """Level raising with exact integer coefficients, valid for all modes.

    initial is a tuple of pairs (n, poly), poly a tuple of integers: the
    coefficients of 1, c, c**2, ... where c = (1-lambda_j)/lambda_j.  An
    argument n < b**(ell-1) is replaced by its children
        b*n + a                      with factor  1   (0 <= a < b),
        b**h*n + s                   with factor  c,
        b**(h+1)*n + a*b**h + s      with factor -c   (0 <= a < b),
        b**p*n + n(w)                with factor  c,
    for each overlap period h of w, s being the value of the last h
    digits of w.  The children of n are at least b*n, so taking the
    arguments in increasing order (a heap) raises each argument once,
    after all its contributions have arrived.  Contributions are merged
    exactly.  Returns a tuple of triples (N, pos, neg), N increasing,
    where pos and neg are the positive coefficients and the absolute
    values of the negative ones; the structure depends on (b, w, ell)
    and the initial arguments only, not on the mode j.
    """
    p = len(w)
    nw = _word_value(b, w)
    right = [(h, _word_value(b, w[p-h:])) for h in range(1, p)
             if w[:p-h] == w[h:]]
    cutoff = b**(ell-1)
    table = {}
    pending = []

    def add(n, poly, shift, sign):
        current = table.get(n)
        if current is None:
            current = table[n] = []
            if n < cutoff:
                heappush(pending, n)
        need = len(poly)+shift
        if len(current) < need:
            current.extend([0]*(need-len(current)))
        for i, coefficient in enumerate(poly):
            current[i+shift] += sign*coefficient

    for n, poly in initial:
        add(n, poly, 0, 1)
    while pending:
        n = heappop(pending)
        poly = table.pop(n)
        if not any(poly):
            continue
        for a in range(b):
            add(b*n+a, poly, 0, 1)
        for h, suffix in right:
            add(b**h*n+suffix, poly, 1, 1)
            for a in range(b):
                add(b**(h+1)*n+a*b**h+suffix, poly, 1, -1)
        add(b**p*n+nw, poly, 1, 1)
    result = []
    for N in sorted(table):
        poly = table[N]
        if any(poly):
            result.append((N, tuple(max(x, 0) for x in poly),
                           tuple(max(-x, 0) for x in poly)))
    return tuple(result)


class _ModalEngine:
    """Triangular recurrences for the modal moments and coefficients.

    Both base moment recurrences have diagonal D(b**(-m-1)); their
    off-diagonal coefficients use different word-overlap orientations.
    Recurrence values cached at an earlier working precision are reused
    as later modes lower their precision.  A precision retry constructs
    a new engine rather than upgrading the cached numbers in place.
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
        # Signed shifts in the base-measure moment recurrences.
        self._vterms = self._terms(self.right)
        self._uterms = self._terms(self.left, left=True)
        self.v = [self.mass_v]
        self.u = [self.mass_u]
        self.M = {}        # M[j] stores moments of the tilted measure nu_j.
        self.R = [self.mass_v]  # Right modal coefficients for eta_0.
        self.U = [self.mass_u]  # Kempner modal coefficients for mu_0.
        self._lams = {}
        self._diagonal = {}
        self._coefficient_rows = {}
        self._tables = {}  # raised integer tables, see _table
        self._stats = defaultdict(int)
        self._internal_loss = mp.zero

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
        """Return a shift coefficient; extend its row cache as needed."""
        degree = row-col
        key = (which, row)
        cached = self._coefficient_rows.get(key)
        if cached is None:
            t = mp.power(self.b, -row-1)
            terms = self._vterms if which == 'v' else self._uterms
            # The row fixes the shift weights; only the degree varies.
            weighted = [(integer, sign*t**length) for length, integer, sign in terms]
            initial_e = mp.fsum(weight for _, weight in weighted[self.b:])
            initial_d = mp.one + initial_e + mp.fsum(weight for _, weight in weighted[:self.b])
            cached = ([initial_d], [initial_e], weighted)
            self._coefficient_rows[key] = cached
        d_values, e_values, weighted = cached
        while len(d_values) <= degree:
            d = len(d_values)
            a_parts, e_parts = [], []
            for index, (integer, weight) in enumerate(weighted):
                if integer:
                    contribution = weight*integer**d
                    (a_parts if index < self.b else e_parts).append(contribution)
            d_values.append(mp.fsum(a_parts + e_parts))
            e_values.append(mp.fsum(e_parts))
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
            overlap = mp.fsum(t**h for h in self.overlaps)
            C = mp.one+overlap
            D = (1-self.b*t)*C+t**self.p
            self._lams[j] = ((1-self.b*t)*overlap+t**self.p)/D
        return self._lams[j]

    def _signed_sum(self, values):
        """Sum by sign and track the largest local cancellation indicator."""
        positives, negatives = [], []
        for x in values:
            if x > 0:
                positives.append(x)
            elif x < 0:
                negatives.append(-x)
        plus, minus = mp.fsum(positives), mp.fsum(negatives)
        result = plus-minus
        if plus and minus:
            loss = mp.log10((plus+minus)/abs(result)) if result else mp.inf
            self._internal_loss = max(self._internal_loss, loss)
        return result

    def _extend_basic(self, order):
        """Extend eta_0 moments and, only for k=0, autonomous mu_0 moments."""
        for m in range(len(self.v), order+1):
            d0 = self._D(m+1)
            vs, us = [], []
            for i in range(m):
                factor = comb(m,i)
                dv, _ = self._de(m,i,'v')
                vs.append(factor*dv*self.v[i])
                if self.need_kempner:
                    du, _ = self._de(m,i,'u')
                    us.append(factor*du*self.u[i])
            self.v.append(-self._signed_sum(vs)/d0)
            if self.need_kempner:
                forcing = mp.fsum(mp.power(self.b, -(m+1)*h)*(value**m)
                                  for h,value in self.left)
                self.u.append(self._signed_sum([forcing] + [-x for x in us])/d0)
            self._stats['basic_moments'] += 1

    def _mom(self, j, q):
        """Moment of the tilted probability measure nu_j (Part 2)."""
        moments = self.M.setdefault(j, [mp.one])
        lamj = self.lam(j)
        while len(moments) <= q:
            r = len(moments)
            row = j+r-1
            summands = []
            for i in range(r):
                col = j+i-1
                d, e = self._de(row,col)
                summands.append(comb(r,i)*(e-lamj*d)*moments[i])
            denominator = self._D(j+r)*(lamj-self.lam(j+r))
            moments.append(self._signed_sum(summands)/denominator)
            self._stats['nu_moments'] += 1
        return moments[q]

    def _basis_coefficients(self, m):
        """Build right coefficients and, for Kempner sums, left ones."""
        self._extend_basic(m)
        while len(self.R) <= m:
            index = len(self.R)
            r_parts = [self.v[index]]
            u_parts = [self.u[index]] if self.need_kempner else []
            for i in range(index):
                c = comb(index,i)*self._mom(i+1,index-i)
                r_parts.append(-c*self.R[i])
                if self.need_kempner:
                    u_parts.append(-c*self.U[i])
            self.R.append(self._signed_sum(r_parts))
            if self.need_kempner:
                self.U.append(self._signed_sum(u_parts))
        return self.R[m], self.U[m] if self.need_kempner else None

    def _table(self, name):
        """Raised table of one group of initial arguments (see _raised_table).

        'k1': the arguments of A_j for k >= 1, with the factor lambda_j
        taken out (lambda_j - 1 = -lambda_j c);  'U' and 'R': the arguments
        multiplied by the Kempner coefficients U and R for k = 0.
        """
        if name not in self._tables:
            b, w, nw = self.b, self.w, self.nw
            initial = defaultdict(list)
            if name == 'k1':
                for d in range(1, b):
                    initial[d].append((1,))
                for h, suffix in self.right:
                    for d in range(1, b):
                        initial[d*b**h+suffix].append((0, -1))
                if w[0]:
                    initial[nw].append((0, 1))
            elif all(x == 0 for x in w):
                if name == 'U':
                    for n in range(b, b*b):
                        initial[n].append((1,))
                else:
                    for d in range(1, b):
                        initial[d*b**self.p].append((-1,))
            elif name == 'U':
                for d in range(1, b):
                    initial[d].append((1,))
            elif w[0]:
                initial[nw].append((-1,))
            merged = []
            for n in sorted(initial):
                poly = [0]*max(len(x) for x in initial[n])
                for x in initial[n]:
                    for i, coefficient in enumerate(x):
                        poly[i] += coefficient
                merged.append((n, tuple(poly)))
            self._tables[name] = _raised_table(b, tuple(w), self.ell,
                                               tuple(merged))
        return self._tables[name]

    def linear_A(self, j, groups, tolerance):
        """Evaluate the Newton sums of mode j over the raised arguments.

        groups is a list of (scalar, table): the weight of a raised
        argument N is the sum of scalar*P(c) over the groups, P the
        integer polynomial of N in the table and c = (1-lambda_j)/lambda_j.
        As c > 0, the positive and negative coefficients of P give the
        positive and negative parts of the weight; their cancellation is
        recorded by _signed_sum.  Newton-series contributions are then
        accumulated by sign, one row at a time.
        """
        lam = self.lam(j)
        c = (1-lam)/lam
        plus = defaultdict(lambda: mp.zero)
        minus = defaultdict(lambda: mp.zero)
        for scalar, table in groups:
            if not scalar:
                continue
            size = abs(scalar)
            for N, pos, neg in table:
                P, M = _horner(pos, c), _horner(neg, c)
                if scalar < 0:
                    P, M = M, P
                if P:
                    plus[N] += size*P
                if M:
                    minus[N] += size*M
        positive_args, negative_args = [], []
        for N in sorted(plus.keys() | minus.keys()):
            weight = self._signed_sum([plus[N], -minus[N]])
            if weight > 0:
                positive_args.append((N, weight))
            elif weight < 0:
                negative_args.append((N, -weight))
        if not positive_args and not negative_args:
            return mp.zero, mp.zero, mp.zero
        Np = [N for N, _ in positive_args]
        Nm = [N for N, _ in negative_args]
        wp = [x for _, x in positive_args]
        wm = [x for _, x in negative_args]
        pp = [mp.power(N, -j) for N in Np]
        pm = [mp.power(N, -j) for N in Nm]
        positives, negatives = [], []
        q = 0
        while True:
            Sp = mp.fdot(wp, pp) if wp else mp.zero
            Sm = mp.fdot(wm, pm) if wm else mp.zero
            # Omitted rows: binomial*sum |w_N| N^(-j-q), moments being <= 1.
            if q and comb(j+q-1, q)*(Sp+Sm) < tolerance:
                break
            if q > 20000:
                raise ArithmeticError('Newton series did not reach requested tolerance')
            signed_factor = ((-1 if q % 2 else 1)*comb(j+q-1, q)
                             * self._mom(j, q))
            if signed_factor >= 0:
                positives.append(signed_factor*Sp)
                negatives.append(signed_factor*Sm)
            else:
                positives.append(-signed_factor*Sm)
                negatives.append(-signed_factor*Sp)
            self._stats['newton_rows'] += 1
            pp = [x/N for x, N in zip(pp, Np)]
            pm = [x/N for x, N in zip(pm, Nm)]
            q += 1
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
            left, positive, negative = self.linear_A(
                j, [(lam, self._table('k1'))], tolerance)
            if m%2:
                left, positive, negative = -left, negative, positive
            term = left*multiplier
            if multiplier < 0:
                positive, negative = negative, positive
            magnitude = abs(multiplier)
            return term, (lam, left, R), positive*magnitude, negative*magnitude
        term, positive, negative = self.linear_A(
            j, [(U, self._table('U')), (R, self._table('R'))], tolerance)
        if m%2:
            term, positive, negative = -term, negative, positive
        return term, None, positive, negative


def _precision_for_term(magnitude, scale, full_dps):
    """Next-mode precision budget before cancellation adjustments."""
    if not magnitude:
        return 12
    return max(12, min(full_dps, int(ceil(full_dps +
                 min(mp.zero, mp.log10(magnitude/scale))))))


def _modalirwin_once(b: int, w, k: int, nboffigures=100,
                     ell=None, verbose=False, trunc=False, *, max_modes=1000,
                     modalcoeffs=False, _guard=12) -> str:
    """Evaluate modes with size-based budgets and selective retries.

    The stopping test uses a numerical tail indicator, not a certified
    remainder bound; max_modes and retry limits can raise ArithmeticError.
    """
    if int(b) != b or b < 2:
        raise ValueError('b must be an integer >= 2')
    b = int(b)
    digits = _digits(b, w)
    if int(k) != k or k < 0:
        raise ValueError('k must be a nonnegative integer')
    k = int(k)
    if int(nboffigures) != nboffigures or nboffigures < 6:
        raise ValueError('nboffigures must be >= 6')
    nboffigures = int(nboffigures)
    if ell is None:
        # Fastest levels at 40 digits for b = 2, 3, 4, measured on
        # October 12, 2026; other bases keep the earlier default 2.
        ell = {2: 5, 3: 4, 4: 3}.get(b, 2)
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
    full_dps = nboffigures + _guard
    decreasing_precision = nboffigures > 30
    terms = []
    decreasing_run = 0
    cancellation_max = mp.zero
    retries = 0
    minimum_dps = full_dps
    maximum_dps = full_dps
    previous_internal_loss = mp.zero
    previous_work_dps = full_dps
    full_cache_required = False
    cache_rebuilds = 0
    suffix_precision_floor = 0
    with mp.workdps(full_dps):
        scale = mp.power(b, len(digits))*mp.log(b)
        target = scale*mp.power(10, -full_dps)
        engine = _ModalEngine(b, digits, ell, need_kempner=(k == 0))
        answer = mp.zero
        if k == 0 and all(x == 0 for x in digits):
            constant = mp.fsum(mp.one/a for a in range(1, b))
            answer += constant
            if modalcoeffs:
                _print_fields([('constant', mp.nstr(constant, nboffigures))])
                print()
        for j in range(1, max_modes+1):
            budget = (_precision_for_term(abs(terms[-1]), scale, full_dps)
                      if terms and decreasing_precision else full_dps)
            work_dps = max(budget, suffix_precision_floor)
            if decreasing_precision and terms:
                if mp.isfinite(previous_internal_loss):
                    work_dps = min(previous_work_dps, full_dps, budget +
                                   max(0, int(ceil(previous_internal_loss))) + 4)
                else:
                    work_dps = full_dps
            if full_cache_required:
                work_dps = max(full_dps, suffix_precision_floor)
            elif work_dps > previous_work_dps:
                # A precision increase cannot improve cached low-precision
                # dependencies. Rebuild at full precision and keep it.
                full_cache_required = True
                work_dps = full_dps
                with mp.workdps(work_dps):
                    engine = _ModalEngine(b, digits, ell,
                                          need_kempner=(k == 0))
                cache_rebuilds += 1
            local_retries = 0
            growth_rechecked = False
            while True:
                # Cache values have the precision of the attempt that
                # originally computed them; a retry uses a fresh engine.
                engine._internal_loss = mp.zero
                with mp.workdps(work_dps):
                    term, factors, positive, negative = engine.modal_term(
                        j, k, target/(100 if j == 1 else 10),
                        modalcoeffs=modalcoeffs)
                    gross = positive+negative
                    magnitude = abs(term)
                    loss = (mp.log10(gross/magnitude) if gross and magnitude
                            else (mp.inf if gross else mp.zero))
                    loss = max(loss, engine._internal_loss)
                    # Estimate cancellation in this term and in its
                    # recurrences before deciding whether to retry.
                    if mp.isfinite(loss) and magnitude:
                        # A smaller contribution admits a lower absolute
                        # accuracy budget, before accounting for cancellation.
                        needed = max(budget, int(ceil(
                            nboffigures + _guard + mp.log10(magnitude/scale) + loss)))
                    elif gross:
                        needed = max(work_dps, full_dps)
                    else:
                        needed = budget
                    if terms and magnitude > abs(terms[-1]) and not growth_rechecked:
                        # Recheck an unexpectedly larger mode with fresh
                        # dependencies: its growth could be numerical.
                        extra_growth = (max(1, int(ceil(loss)))
                                        if mp.isfinite(loss) else 12)
                        needed = max(needed, work_dps+extra_growth)
                        growth_rechecked = True
                    # Guard against a small net term arising from large
                    # opposite-sign contributions at insufficient precision.
                    rounding_proxy = gross*mp.power(10, -work_dps)
                    if rounding_proxy >= target and magnitude < target:
                        needed = max(needed, work_dps + max(1, int(ceil(
                            mp.log10(rounding_proxy/target)))+2))
                if (decreasing_precision and not full_cache_required and
                        (not mp.isfinite(engine._internal_loss) or
                         work_dps-engine._internal_loss < 12)):
                    # Unresolved internal cancellation at reduced precision
                    # requires rebuilding from the original precision.
                    full_cache_required = True
                    needed = max(needed, full_dps)
                if needed <= work_dps:
                    break
                if local_retries >= 4 or needed > full_dps+120:
                    raise ArithmeticError(
                        f'mode {j}: cancellation requires excessive precision')
                work_dps = max(work_dps+1, needed)
                retries += 1
                local_retries += 1
                with mp.workdps(work_dps):
                    engine = _ModalEngine(b, digits, ell, need_kempner=(k == 0))
                cache_rebuilds += 1
            if work_dps > full_dps:
                # Keep the higher precision for subsequent cached work.
                suffix_precision_floor = max(suffix_precision_floor, work_dps)
            previous_internal_loss = engine._internal_loss
            previous_work_dps = work_dps
            minimum_dps = min(minimum_dps, work_dps)
            maximum_dps = max(maximum_dps, work_dps)
            cancellation_max = max(cancellation_max, loss)
            if terms:
                decreasing_run = (decreasing_run+1 if
                                  magnitude < abs(terms[-1]) else 0)
            # Add the accepted term at the full accumulator precision.
            answer += term
            terms.append(term)
            tail_bound = 2*magnitude
            if modalcoeffs:
                if j > 1:
                    print()
                shown = min(nboffigures, max(2, work_dps-12))
                signed_term = ('+' if term >= 0 else '') + mp.nstr(term, shown)
                fields = [('n', j), ('budget digits', budget),
                          ('working digits', work_dps)]
                if k >= 1:
                    lam, left, right = factors
                    with mp.workdps(work_dps):
                        product = left*right
                        prime = product/lam
                    fields.extend([
                        ('lambda', mp.nstr(lam, shown)),
                        ('L', mp.nstr(left, shown)),
                        ('R', mp.nstr(right, shown)),
                        ('C', mp.nstr(product, shown)),
                        ("C'", mp.nstr(prime, shown)),
                    ])
                fields.append(('T', signed_term))
                _print_fields(fields)
            elif verbose:
                if j > 1:
                    print()
                _print_fields([
                    ('mode', j),
                    ('term', mp.nstr(term, 12)),
                    ('|term|', mp.nstr(magnitude, 5)),
                    ('working digits', work_dps),
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
                ('local precision retries', retries),
                ('cache rebuilds', cache_rebuilds),
                ('elapsed time', f'{time.perf_counter()-started:.3f}s'),
                ('estimated inclusive-tail bound', mp.nstr(tail_bound, 5)),
                ('working decimal digits', full_dps),
                ('minimum modal decimal digits', minimum_dps),
                ('maximum modal decimal digits', maximum_dps),
                ('guard decimal digits', _guard),
                ('successive decreases', decreasing_run),
                ('reference scale', mp.nstr(scale, 6)),
                ('threshold', mp.nstr(target, 5)),
            ])
            print('Tail bound is estimated, inclusive of the last term.')
    return result


def modalirwin(b: int, w, k: int, nboffigures: int = 100,
               ell=None, verbose=False, trunc=False, *, max_modes=1000,
               modalcoeffs=False) -> str:
    """Compute I(b,w,k) using Burnol's published modal expansion.

    Return ``nboffigures`` decimal significant figures. ``ell`` sets the
    auxiliary level-raising threshold.  ``modalcoeffs`` prints modal terms
    and, for k >= 1, their coefficient factors and per-mode precisions.

    Above 30 requested digits, the preceding mode sets a size-based
    precision budget for the next mode, relative to b**len(w)*log(b).
    Cancellation checks can override that budget and rebuild the cache.
    The stopping threshold includes twelve decimal guard digits.
    Cancellation and tail indicators estimate accuracy but do not certify
    twelve correct guard digits or provide a rigorous remainder bound.
    The estimated tail bound includes the final computed term.
    """
    return _modalirwin_once(
        b, w, k, nboffigures, ell, verbose, trunc,
        max_modes=max_modes, modalcoeffs=modalcoeffs)


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
