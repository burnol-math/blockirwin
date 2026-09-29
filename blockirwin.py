# -*- mode: python ; coding: utf-8; -*-
"""blockirwin.py -- block-Irwin sums I(b, w, k) to arbitrary precision

THIS FILE WAS PRODUCED ON MONDAY SEP 28, 2026 BY CLAUDE OPUS 5.5
AND CLAUDE SONNET 5
FROM THE FORMULAS OF THE PAPER BY JEAN-FRANÇOIS BURNOL CITED BELOW.
The mathematics is fully authored by Jean-François Burnol.

This code is licensed under the
Creative Commons Attribution-ShareAlike 4.0 International License.
Full license text: https://creativecommons.org/licenses/by-sa/4.0/

For b > 1 an integer, w a non-empty word of base-b digits (leading
zeros allowed) and k >= 0, the block-Irwin sum is

    I(b, w, k) = sum of 1/n over the integers n > 0 whose minimal
                 base-b representation contains exactly k occurrences
                 of w (occurrences counted with overlaps).

Reference:

    J.-F. Burnol, Block-count constrained harmonic sums: spectral
    expansion and block-directed Euler--Maclaurin, arXiv:2609.17045,
    Part 2: Proposition 11, Proposition 12 (level raising), and
    Section 30 (finite triangular moment recurrences).

    J.-F. Burnol, The modal expansion of Kempner sums, arXiv:2609.31986,
    continuation of the above for the modal expansion when k = 0 (the
    modal expansion is not used here).

Method:

1) Proposition 11 expresses I(b, w, k) as a finite integer combination
   of the Stieltjes values V_j(N) = int dη_j(x)/(N+x), j in {k, k-1},
   of the measures η_j on [0,1) (the "Burnol measures" of the verbose
   output), plus nothing else.

2) Proposition 12 (level raising) replaces V_j(N) by V_j and V_{j-1}
   values at arguments >= bN, plus the exact term 1/N if j = 0.  It is
   applied to every argument having fewer than ell digits in base b, so
   that in the end all arguments N satisfy N >= b^(ell-1).  Identical
   (j, N) pairs are merged, and so are the exact terms.

3) Each remaining V_j(N) is the alternating series
       V_j(N) = sum_{m>=0} (-1)^m v_{j;m} / N^(m+1)
   where v_{j;m} is the m-th moment of η_j.  As η_j has total mass b^p
   (p = |w|) and support in [0,1), 0 <= v_{j;m} <= b^p and the terms
   decrease, so each series is truncated a priori, with a rigorous
   bound, after a number of terms depending on N.  The moments are
   computed by the finite triangular recurrences of Section 30
   (Proposition 13), all moments of index m for 0 <= j <= k at once.

4) Precision management, as in kempner.py: the moments of index m only
   enter through terms of size at most b^p/N^(m+1), so the working
   precision decreases by log2(min N) bits at each step m.  Numbers
   stored at an earlier step (moments, power sums, constants) are
   rounded down to the running precision whenever it has decreased by
   10%, so that every multiplication involves mantissas of about the
   running precision only.

About mpmath and pure powers of 2: mpmath 1.4.0 and 1.4.1 with
their pure-Python backend (i.e. when gmpy2 is not installed) convert
Python integers which have many trailing zero bits (e.g. 2**10000)
to mpf numbers extremely slowly (quadratic time; see
https://github.com/mpmath/mpmath/issues/1073, fixed in the
development version, unreleased as of 2026-09-28; mpmath 1.3.0 is not
affected).  Such integers arise naturally here, in particular for
bases b = 2 and b = 16.  This script never lets mpmath convert such
an integer: every conversion strips trailing zero bits first (see
_mpf_from_int).  Installing gmpy2 is recommended anyhow for speed.

No attempt at parallelization is made.

Usage (interactive):

    >>> from blockirwin import blockirwin
    >>> blockirwin(10, "42", 1, 50)
    >>> blockirwin(2, "11", 0, 100, verbose=True)
    >>> blockirwin(16, "ff", 3, 1000, ell=3)

Usage (command line):

    python blockirwin.py b w k [nboffigures [ell]]
"""

__version__ = "1.2.0"
__date__ = "2026-09-29"

from typing import List, Optional
from math import ceil, log, log2
from heapq import heappush, heappop
import time
from mpmath import mp, mpf, nstr, inf

_FIXED = dict(strip_zeros=False, min_fixed=-inf, max_fixed=inf)  # never
# switch to scientific notation, so a fixed number of trailing characters
# always corresponds to a fixed number of extra decimal digits (needed for
# the trunc=True slicing below to be safe at any magnitude).

_DIGITS = "0123456789abcdefghijklmnopqrstuvwxyz"


# ----------------------------------------------------------------------
# small utilities
# ----------------------------------------------------------------------

def _mpf_from_int(n: int) -> mpf:
    """mpf(n) rounded to current precision, without handing to mpmath an
    integer with trailing zero bits (see mpmath issue 1073)."""
    if n == 0:
        return mpf(0)
    t = (n & -n).bit_length() - 1
    if t:
        return mp.ldexp(mpf(n >> t), t)
    return mpf(n)


def _parse_block(b: int, w) -> List[int]:
    if isinstance(w, str):
        assert b <= 36, "a string w requires b <= 36, else use a list of ints"
        digits = [int(ch, 36) for ch in w]
    elif isinstance(w, int):
        raise TypeError("w must be a string of base-b digits or a list of "
                        "digits (an int would lose leading zeros)")
    else:
        digits = [int(d) for d in w]
    assert len(digits) > 0, "w must be non-empty"
    for d in digits:
        assert 0 <= d < b, f"{d} is not a base {b} digit"
    return digits


def _tostr(b: int, digits: List[int]) -> str:
    if b <= 36:
        return "".join(_DIGITS[d] for d in digits)
    return "[" + ",".join(map(str, digits)) + "]"


def _num2str(b: int, n: int) -> str:
    L = []
    while n:
        n, r = divmod(n, b)
        L.append(r)
    return _tostr(b, L[::-1])


def _val(b: int, digits: List[int]) -> int:
    n = 0
    for d in digits:
        n = n * b + d
    return n


def _overlap_periods(w: List[int]) -> List[int]:
    """The overlap periods of w: the integers h with 0 < h < p such that
    w[:p-h] == w[h:], i.e. w has a border of length p-h.  The suffix c_h
    of w of length h is then the corresponding element of C_w^+: two
    occurrences of w at distance h overlap, and w c_h contains w twice."""
    p = len(w)
    return [h for h in range(1, p) if w[:p - h] == w[h:]]


def default_ell(b: int, w, k: int, nboffigures: int) -> int:
    """Rule of thumb for the level ell (see README.md).

    Measured on September 28, 2026 with pure-Python mpmath 1.4.1, for
    b in {2, 3, 4, 5, 10, 16}, 100 to 800 digits, 0 <= k <= 10 and blocks
    of length 1 to 5: the time at the chosen level is within 10% of the
    time at the best level in 251 of 273 measured cases, and within 29%
    in all of them.  Above 800 digits the values are not measured.  Other
    bases use a formula which fits the measured bases 3 to 16 within one
    level, but which has not been measured for them.
    """
    p = len(w)
    d = nboffigures
    if b == 2:
        if p == 1:
            if k == 0:
                e = round(0.4 * d + 10)
            elif k == 1:
                e = round(15 + 0.05 * d)
            elif k == 2:
                e = round(13 + 0.028 * d)
            else:
                e = 10 if d <= 150 else 12 if d <= 450 else 13
        elif k == 0:
            e = 11 if d <= 150 else 12 if d <= 300 else 14
        else:
            e = 10 if d <= 150 else 11 if d <= 300 else 12
    elif b == 3:
        e = (7 if d <= 150 else 8) if p == 1 else (6 if d <= 300 else 7)
    elif b == 4:
        e = 5 if d <= 300 else 6
    elif b == 5:
        e = 4 if d <= 150 else 5
    elif b == 10:
        if d < 400:
            e = 3
        elif k == 0:
            e = 3 if d <= 600 else 4
        elif d < 600:
            e = 4 if k >= 2 else 3
        else:
            e = 4
    elif b == 16:
        e = 3
    else:
        e = 1 + round(log(d + 100) / log(b))
    return max(2, e)


# ----------------------------------------------------------------------
# Propositions 11 and 12: the table of Stieltjes values
# ----------------------------------------------------------------------

def _build_table(b: int, w: List[int], k: int, ell: int):
    """Propositions 11 and 12.

    Returns (table, const): table maps (j, N) to the integer coefficient
    of V_j(N), const maps n to the integer coefficient of 1/n, so that
    I(b, w, k) = sum const[n]/n + sum table[j, N] V_j(N).  All N in table
    have at least ell digits in base b.  The arguments are examined in
    increasing order (a heap), which is safe since level raising only
    creates arguments at least b times larger."""
    p = len(w)
    nw = _val(b, w)
    # (|c|, n(c)) for c in C_w^+, the complements of the borders of w
    SC = [(h, _val(b, w[p - h:])) for h in _overlap_periods(w)]

    bpow = [b**i for i in range(p + 1)]     # b^0, ..., b^p

    table = {}                  # (j, N) -> integer coefficient of V_j(N)
    const = {}                  # n -> integer coefficient of 1/n
    pending = []                # heap of the arguments not yet examined
    seen = set()

    def add(j, N, c):
        if j < 0:           # V_{-1} = 0
            return
        key = (j, N)
        table[key] = table.get(key, 0) + c
        if N not in seen:
            seen.add(N)
            heappush(pending, N)

    def addDelta(j, N, c):  # c * (V_j - V_{j-1})(N)
        add(j, N, c)
        add(j - 1, N, -c)

    # Proposition 11: sum over the leading digits a, then over the words
    # ac with c in C_w^+, minus the word w itself if it has no leading 0.
    for a in range(1, b):
        add(k, a, 1)
    for (h, nc) in SC:
        for a in range(1, b):
            addDelta(k, a * bpow[h] + nc, 1)
    if w[0] != 0:
        addDelta(k, nw, -1)

    # Proposition 12, applied to all arguments with fewer than ell digits.
    # New arguments are >= b N, so the heap pops them in a safe order.
    threshold = b**(ell - 1)
    while pending and pending[0] < threshold:
        n = heappop(pending)
        for j in range(k, -1, -1):
            c = table.pop((j, n), 0)
            if c == 0:
                continue
            if j == 0:
                const[n] = const.get(n, 0) + c
            bn = b * n
            for a in range(b):
                add(j, bn + a, c)
            for (h, nc) in SC:
                addDelta(j, bpow[h] * n + nc, -c)
                base = bpow[h] * bn + nc          # n followed by a then c
                for a in range(b):
                    addDelta(j, base + a * bpow[h], c)
            addDelta(j, bpow[p] * n + nw, -c)

    table = {key: c for key, c in table.items() if c != 0}
    const = {n: c for n, c in const.items() if c != 0}
    return table, const


# ----------------------------------------------------------------------
# main function
# ----------------------------------------------------------------------

def blockirwin(b: int,
               w,
               k: int,
               nboffigures: int = 100,
               ell: Optional[int] = None,
               verbose=False,
               trunc: bool = False) -> str:
    """Computes the block-Irwin sum I(b, w, k)

    Sum of 1/n over positive integers n whose base-b representation has
    exactly k occurrences (overlaps allowed) of the block w.  Uses
    Propositions 11 and 12 and the moment recurrences of Section 30 of
    arXiv:2609.17045 (see also its continuation arXiv:2609.31986, on the
    modal expansion when k = 0).

    :param b int: the base, at least 2
    :param w: the block, either a string of base-b digits (0-9a-z, needs
              b <= 36) or a list of ints; leading zeros are significant
    :param k int: the number of occurrences, at least 0
    :param nboffigures int: (optional, defaults to 100) number of
                            asked-for significant figures
    :param ell int: (optional) the level, at least 2: all Stieltjes
                    arguments are raised to at least ell base-b digits.
                    Defaults to default_ell(b, w, k, nboffigures).
    :param verbose: (optional, defaults to False) if True, print details
                    on the computation; the list of arguments is
                    abridged if longer than 40 lines, unless verbose=2
    :param trunc bool: (optional, defaults to False) whether to truncate
                       rather than round
    :rtype: str
    :return: nboffigures significant figures of I(b, w, k)
    """
    t_start = time.time()
    assert int(b) == b and b > 1, f"base {b} must be an integer > 1"
    b = int(b)
    w = _parse_block(b, w)
    p = len(w)
    assert int(k) == k and k >= 0, f"k = {k} must be an integer >= 0"
    k = int(k)
    assert int(nboffigures) == nboffigures, "nboffigures must be an integer"
    nboffigures = int(nboffigures)
    assert nboffigures > 5, "ask for at least 6 significant figures"
    if ell is None:
        ell = default_ell(b, w, k, nboffigures)
    assert int(ell) == ell and ell >= 2, f"level ell = {ell} must be >= 2"
    ell = int(ell)

    wstr = _tostr(b, w)
    nw = _val(b, w)

    def say(*args):
        if verbose:
            print(*args)

    say(f"I(b, w, k) with b = {b}, w = \"{wstr}\" (n(w) = {nw}),"
        f" p = |w| = {p}, k = {k}")

    # Every positive integer starts with digit 1 in base 2.
    if b == 2 and w == [1] and k == 0:
        say("empty sum: every positive integer has a leading binary 1")
        return "0"

    H = _overlap_periods(w)
    if verbose:
        if H:
            for h in H:
                print(f"  border \"{_tostr(b, w[:p - h])}\""
                      f" with complement \"{_tostr(b, w[p - h:])}\""
                      f" (overlap period {h})")
        else:
            print("  w has no border")

    # ---------------- precision targets
    # P_target bits give nboffigures + 1 decimal digits.  The truncations
    # of the series are chosen so that their total error is at most
    # 2^(-P_eps) times (an estimate of) the result.  P_max adds guard bits
    # for the rounding errors of the arithmetic, which is done at P_max
    # bits for the terms of index 0 and at decreasing precision after.
    P_target = ceil((nboffigures + 1) * log2(10))
    P_eps = P_target + 3
    P_max = P_eps + 20 + ceil(log2(k + 1))
    mp.prec = P_max

    # ---------------- Propositions 11 and 12
    table, const = _build_table(b, w, k, ell)

    # One list of the Stieltjes arguments: (cut, N, [(j, c), ...]).
    byN = {}
    for (j, N), c in table.items():
        byN.setdefault(N, []).append((j, c))
    del table
    say(f"level ell = {ell}: {sum(len(cs) for cs in byN.values())} Stieltjes"
        f" values V_j(N) at {len(byN)} distinct arguments N,"
        f" plus {len(const)} exact terms 1/n")

    S = mp.fsum(_mpf_from_int(c) / _mpf_from_int(n)
                for n, c in const.items())
    if not byN:
        say("no Stieltjes values remain")
        return nstr(S, nboffigures, **_FIXED)

    bp = b**p                     # total mass of every η_j
    # Rough size of the result, from the terms of index m = 0 of the
    # series (v_{j;0} = b^p for every j): it scales the error target.
    S0 = S + _mpf_from_int(bp) * mp.fsum(_mpf_from_int(c) / _mpf_from_int(N)
                                         for N, cs in byN.items()
                                         for _, c in cs)
    log2eps = float(mp.log(abs(S0), 2)) - P_eps
    log2tau = log2eps - log2(len(byN))
    # For an argument N the terms have absolute value at most A b^p/N^(m+1)
    # with A the sum of the |c|, and the tail of an alternating series with
    # decreasing terms is bounded by its first omitted term.  cut = number
    # of kept terms.  Sorted by decreasing cut: the arguments still active
    # at step m are a prefix of the list.
    args = []
    for N, cs in byN.items():
        A = sum(abs(c) for _, c in cs)
        cut = max(0, ceil((log2(A * bp) - log2tau) / log2(N) - 1 + 1e-9))
        args.append((cut, N, cs))
    del byN
    args.sort(key=lambda t: (-t[0], t[1]))
    n_act = sum(1 for t in args if t[0] > 0)
    Mmax = args[0][0]
    # used[j]: number of moments of η_j entering the series
    # need[j]: number of moments of η_j to compute (η_{j+1} needs η_j)
    used = [0] * (k + 1)
    for cut, _, cs in args:
        for j, _ in cs:
            if cut > used[j]:
                used[j] = cut
    need = used[:]
    for j in range(k - 1, -1, -1):
        need[j] = max(need[j], need[j + 1])

    # ---------------- moments of the measures eta_j (Section 30)
    # E_w = C_w - b^{-1} A C_w + b^{-p} P_w and D_w = I - b^{-1} A + E_w.
    # The words s of E_w carry signs: +c (c in C_w^+), -ac (a a digit),
    # +w.  A word acts through the point x(s) = n(s)/b^|s| in [0,1).
    # Pairs with equal (|s|, n(s)) and opposite signs cancel exactly
    # (+c_{h+1} and -ac_h with a the digit of w before c_h, and +w with
    # -w_1 c_{p-1}), so signs are accumulated and zeros dropped.
    signs = {}
    for h in H:
        nc = _val(b, w[p - h:])
        signs[(h, nc)] = signs.get((h, nc), 0) + 1
        for a in range(b):
            key = (h + 1, a * b**h + nc)
            signs[key] = signs.get(key, 0) - 1
    signs[(p, nw)] = signs.get((p, nw), 0) + 1
    words = {}
    for (L, n), sg in signs.items():
        if sg:
            words.setdefault(L, []).append((sg, n))
    del signs

    # Notation.  For the words s of E_w, grouped by length L:
    #   y_{L,r} = sum_s sign(s) x(s)^r,
    # and for the digits a of Sigma_b (the -b^{-1} A part of D_w):
    #   g_r = sum_{a=0}^{b-1} (a/b)^r        (with 0^0 = 1).
    # Proposition 13 then reads, for m >= 0 and 0 <= j <= k, with
    # u = m+1-r, delta = 1 if j = m = 0 and delta = 0 otherwise,
    # and e_r(m) = sum_L y_{L,r} b^{-L u}:
    #   v_{j;m} + sum_{r=0}^m binom(m,r) (e_r(m) - g_r b^{-u}) v_{j;m-r}
    #           - sum_{r=0}^m binom(m,r) e_r(m) v_{j-1;m-r}  =  delta.
    # The r = 0 term of the first sum is d_{m,m} - 1 times v_{j;m}; the
    # r = 0 term of the second sum involves v_{j-1;m}, which is computed
    # first since j increases.  Kept across steps: only y_{L,r}, g_r and
    # the moments v_{j;i}, which the convolutions over r need; binom(m,r)
    # and the powers b^{-L u} are produced on the fly during the pass
    # over r.  Everything is an mpf at the running precision.
    Ls1 = sorted(set(words) | {1})       # word lengths, and 1 for g
    nL = len(Ls1)
    i1 = Ls1.index(1)                    # position of L = 1 in Ls1
    wl = [t for t, L in enumerate(Ls1) if L in words]
    nq = len(wl)                         # number of word lengths
    # E_w is not zero (its word w, or else the words ac_{p-1} with a not
    # equal to w_1, survive the cancellations), so nq >= 1.
    assert nq >= 1
    bL = [_mpf_from_int(b**L) for L in Ls1]                   # b^L
    binvL = [1 / x for x in bL]                               # b^{-L}
    Bm1 = [mpf(1)] * nL                  # b^{-L(m+1)} at step m
    sg_of = [[sg for sg, _ in words[Ls1[t]]] for t in wl]
    x_of = [[_mpf_from_int(n) / bL[t] for _, n in words[Ls1[t]]] for t in wl]
    xpow = [[mpf(1)] * len(words[Ls1[t]]) for t in wl]      # x(s)^m
    y = [[] for _ in wl]                                      # y[q][r]
    da = [mpf(a) / b for a in range(b)]                       # a/b
    dpow = [mpf(1)] * b                                       # (a/b)^m
    g = []                                                    # g[r]
    v = [[] for _ in range(k + 1)]                            # v[j][i]

    # Series.  For the arguments still active (a prefix of args), the
    # integer ips[i] = floor(2^F / N_i^(m+1)) holds N_i^{-m-1} in fixed
    # point; F leaves room for the truncation errors of all the terms.
    # T[j] = sum_i c_{i,j} ips[i] is then multiplied by v_{j;m}.
    Nlist = [N for _, N, _ in args[:n_act]]
    Clist = [cs for _, _, cs in args[:n_act]]
    Amax = max((sum(abs(c) for _, c in cs) for cs in Clist), default=1)
    F = P_max + ceil(log2(max(1, n_act) * (Mmax + 1) * Amax * bp)) + 2
    ips = [(1 << F) // N for N in Nlist]
    prec_delta = log2(min(Nlist)) if n_act else 1.0
    current_prec_float = float(P_max)
    prec_stored = P_max                  # precision of the stored numbers
    t_moments = 0.0
    t_series = 0.0
    jtop = k                             # moments needed for j <= jtop

    def reround(lists):
        # unary plus rounds an mpf to the current precision
        for lst in lists:
            lst[:] = [+x for x in lst]

    for m in range(Mmax):
        prec = max(30, ceil(current_prec_float))
        mp.prec = prec
        if prec < 0.9 * prec_stored:
            reround(v + y + x_of + xpow + [g, da, dpow, binvL, Bm1])
            prec_stored = prec
        t0 = time.time()
        while jtop >= 0 and need[jtop] <= m:
            jtop -= 1

        # ---- power sums of index m, and b^{-L(m+1)}
        if m:
            for q in range(nq):
                xp, xs = xpow[q], x_of[q]
                for i in range(len(xp)):
                    xp[i] *= xs[i]
            for a in range(1, b):        # dpow[0] = 0^m, 1 then 0
                dpow[a] *= da[a]
            dpow[0] = mpf(0)
        for q in range(nq):
            z = mpf(0)
            for sg, xm in zip(sg_of[q], xpow[q]):
                z += xm if sg == 1 else sg * xm
            y[q].append(z)
        z = mpf(0)
        for xm in dpow:
            z += xm
        g.append(z)
        for t in range(nL):
            Bm1[t] *= binvL[t]

        # ---- single pass over r = 0..m for all j <= jtop at once.
        # C = binom(m,r) by one multiplication and one division, and
        # pw[t] = b^{-L u} with u = m+1-r by one multiplication each.  The
        # coefficient of index r is formed once and used at once.
        acc = [mpf(0)] * (jtop + 1)
        if m == 0:
            acc[0] = mpf(1)
        C = mpf(1)
        pw = Bm1[:]
        for r in range(m + 1):
            if r:
                C = C * (m - r + 1) / r
                for t in range(nL):
                    pw[t] *= bL[t]
            e = y[0][r] * pw[wl[0]]
            for q in range(1, nq):
                e += y[q][r] * pw[wl[q]]
            cD = C * (e - g[r] * pw[i1])
            if r == 0:
                D0 = 1 + cD              # d_{m,m}
                cE0 = e                  # C = 1 here
            else:
                i = m - r
                for j in range(jtop + 1):
                    acc[j] -= cD * v[j][i]
                if jtop:
                    cE = C * e
                    for j in range(1, jtop + 1):
                        acc[j] += cE * v[j - 1][i]
        for j in range(jtop + 1):
            z = acc[j]
            if j:
                z += cE0 * v[j - 1][m]
            v[j].append(z / D0)
        t1 = time.time()
        t_moments += t1 - t0

        # ---- series terms of index m, for the arguments still active
        while n_act and args[n_act - 1][0] <= m:
            n_act -= 1
        del ips[n_act:]
        T = [0] * (k + 1)
        for i in range(n_act):
            xi = ips[i]
            for j, c in Clist[i]:
                if c == 1:
                    T[j] += xi
                elif c == -1:
                    T[j] -= xi
                else:
                    T[j] += c * xi
            ips[i] = xi // Nlist[i]
        rowsum = mpf(0)
        for j in range(k + 1):
            if T[j]:
                rowsum += v[j][m] * _mpf_from_int(T[j])
        t_series += time.time() - t1

        mp.prec = P_max
        rowsum = mp.ldexp(rowsum, -F)
        if m % 2:
            S -= rowsum
        else:
            S += rowsum
        current_prec_float -= prec_delta

    mp.prec = P_max

    if verbose:
        print(f"target: {nboffigures} significant figures;"
              f" working precision {P_max} bits, decreasing by"
              f" {prec_delta:.2f} bits per moment")
        _print_arguments(b, args, const, verbose)
        print("moments of the Burnol measures eta_j:")
        for j in range(k, -1, -1):
            if need[j] == 0:
                print(f"  eta_{j}: no moments needed"
                      " (contribution already negligible)")
            else:
                print(f"  eta_{j}: {need[j]} moments computed"
                      f" (m = 0..{need[j] - 1}), {used[j]} used in the series")
        print(f"times: moments {t_moments:.3f}s, series {t_series:.3f}s,"
              f" total {time.time() - t_start:.3f}s")
        S_string = nstr(S, nboffigures + 3, **_FIXED)
        print(f"Last digits are {S_string[-12:-3]}({S_string[-3:]})")
        print(f"I({b}, \"{wstr}\", {k}) {'truncated' if trunc else 'rounded'}"
              f" to {nboffigures} significant figures is:")

    if trunc:
        return nstr(S, nboffigures + 6, **_FIXED)[:-6]
    return nstr(S, nboffigures, **_FIXED)


def _print_arguments(b, args, const, verbose):
    lines = []
    for cut, N, cs in sorted(args, key=lambda t: t[1]):
        coeffs = ", ".join(f"{c:+d} V_{j}" for j, c in sorted(cs, reverse=True))
        lines.append(f"  \"{_num2str(b, N)}\" = {N}: {coeffs}; {cut} terms")
    print(f"Stieltjes arguments N (base {b} = decimal), with the coefficients"
          f" of V_j(N) and number of series terms kept:")
    if verbose != 2 and len(lines) > 40:
        for x in lines[:20]:
            print(x)
        print(f"  ... ({len(lines) - 40} lines skipped, use verbose=2) ...")
        for x in lines[-20:]:
            print(x)
    else:
        for x in lines:
            print(x)
    if const:
        items = [f"{const[n]:+d}/\"{_num2str(b, n)}\"(={n})"
                 for n in sorted(const)]
        print(f"exact terms c/n ({len(items)}):")
        if verbose != 2 and len(items) > 40:
            items = items[:20] + [f"... ({len(items) - 40} skipped) ..."] \
                + items[-20:]
        for i in range(0, len(items), 6):
            print("  " + "  ".join(items[i:i + 6]))

if __name__ == "__main__":
    import sys
    args = sys.argv[1:]
    if 3 <= len(args) <= 5:
        b = int(args[0])
        w = args[1] if b <= 36 else [int(x) for x in args[1].split(",")]
        k = int(args[2])
        nb = int(args[3]) if len(args) > 3 else 100
        ell = int(args[4]) if len(args) > 4 else None
        print(blockirwin(b, w, k, nb, ell))
    else:
        print("""
blockirwin(b, w, k[, nboffigures, ell, verbose, trunc])

Computes I(b, w, k), the sum of 1/n over the positive integers n whose
base-b representation contains exactly k occurrences (overlaps allowed)
of the block w (a string of base-b digits such as "0110" or "ff", or a
list of ints; leading zeros are significant).

The last parameters are optional:
- nboffigures (default 100): number of significant figures
- ell (default from default_ell()): the level parameter
- verbose (False): print details (True, or 2 for complete listings)
- trunc (False): whether to truncate or round.

For use on command line: python blockirwin.py b w k [nboffigures [ell]]
(for b > 36, give w as comma separated digits).
""")
