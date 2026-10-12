# Modal expansions of block-Irwin and Kempner sums — `modalirwin.py`

## Scope and mathematical sources

`modalirwin.py` evaluates the block-Irwin sum `I(b, w, k)` with the
**modal-expansion formulas developed by Jean-François Burnol**. The
occurrences of the digit word `w` in positive integers written in base `b`
are counted with overlaps. Leading zeros in the *word* are significant.

The entire mathematical method comes from these two original papers:

1. Jean-François Burnol, *Block-count constrained harmonic sums: spectral
   expansion and block-directed Euler–Maclaurin*,
   [arXiv:2609.17045](https://arxiv.org/abs/2609.17045), **Part 2**.
   Its Theorem 3 gives the modal expansion for **`k >= 1`**; subsequent
   sections give the right and left modal coefficients, tilted-measure
   moments, and auxiliary Stieltjes level raising.
2. Jean-François Burnol, *The modal expansion of Kempner sums*,
   [arXiv:2609.31986](https://arxiv.org/abs/2609.31986), especially its
   initial sections on **`k = 0`**, the occurrence-free measure, and the
   exceptional forbidden word `0^p`. The autonomous occurrence-free
   moment recurrence appears later in this paper.

**No new mathematical formulas or theoretical error estimates are claimed.**
The code is a numerical realization of Professor Burnol's results.
It is independent of `blockirwin.py`, which implements a different method
based on Propositions 11–13 of the first paper. Both modules can be imported
at the same time; the existing module is not modified.

## Acknowledgement and authorship

**THIS IMPLEMENTATION, `modalirwin.py`, ITS TEST SUITE `test_modalirwin.py`,
AND THE COMPARISON SCRIPT `compare_modalirwin.py` WERE PRODUCED ON
OCTOBER 10, 2026, BY OPENAI CHATGPT (GPT-6), FROM THE TWO PAPERS BY
JEAN-FRANÇOIS BURNOL CITED ABOVE.**

**The mathematics is entirely authored by Professor Jean-François Burnol.**
The numerical implementation was created independently of the method used
in `blockirwin.py`. The latter, its original test suite, and
`README-method1.md` were produced by Claude Opus 5.5 and Claude Sonnet 5
(Anthropic) on September 28–29, 2026, using Professor Burnol's work.
The implementations are complementary, not replacements for each other.

On October 12, 2026, Claude Opus 5.5 (Anthropic) revised `modalirwin.py`
into version 0.6.0.  The auxiliary level raising, previously redone for
every mode and every initial argument with a scan of all arguments at each
step, is now built once per call as polynomials with integer coefficients
in $c=(1-\lambda_j)/\lambda_j$, evaluated for each mode.  Exact
cancellations, which were counted as an infinite loss of digits and forced
full precision in all modes, no longer occur.  The Newton stopping test
reuses the per-row sums, and the default levels were changed (see the API
table).  The test `test_raised_argument_signs_not_merged` of
`test_modalirwin.py`, which checked the former design, was replaced by
`test_raised_table_matches_direct_raising`.

This contribution and the repository are distributed under the existing
**Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)**
license. The full legal text is in [`LICENSE`](LICENSE); see also
[the CC BY-SA 4.0 license page](https://creativecommons.org/licenses/by-sa/4.0/).
Retain attribution and follow the license's share-alike requirements when
redistributing adaptations.

## Requirements and public API

Python 3 and [`mpmath`](https://mpmath.org/), plus Python's standard library.
No imports from `blockirwin.py` are needed inside the modal module.

**Environment used for the reported validations:** CPython **3.13.5**,
`mpmath` **1.3.0**, without `gmpy2`.  The tests have not established
compatibility or performance under `mpmath` 1.4.x.

```python
from modalirwin import modalirwin

modalirwin(b, w, k, nboffigures=100, ell=None,
           verbose=False, trunc=False, *, max_modes=1000,
           modalcoeffs=False)
```

The result is a **decimal string**, compatible with the return format of
`blockirwin`. The parameters are:

| Parameter | Meaning |
| --- | --- |
| `b` | Integer base, at least 2 |
| `w` | Nonempty word: a digit string for bases up to 36, or a sequence of integer digits for any base |
| `k` | Number of occurrences, an integer at least 0 |
| `nboffigures` | Requested number of significant figures, at least 6 (default 100) |
| `ell` | Minimum number of base-`b` digits in the **auxiliary modal Stieltjes arguments** after level raising; default 5 in base 2, 4 in base 3, 3 in base 4 and 2 otherwise |
| `verbose` | Print per-mode magnitudes, recurrence counts, timings, and the heuristic stopping indicator |
| `trunc` | Truncate the decimal result rather than round it |
| `max_modes` | Maximum modal indices to evaluate; raise `ArithmeticError` if no stopping decision is obtained |
| `modalcoeffs` | Print signed modal contributions and, for `k >= 1`, their modal factors |

`ell` applies to the modal method's auxiliary functions and **does not
have the same computational meaning** as the homonymous parameter of
`blockirwin.py`. When comparing the two methods, leave it unspecified
unless you deliberately want to override each method's default.

Examples:

```python
from blockirwin import blockirwin
from modalirwin import modalirwin

print(blockirwin(2, "01", 0, 20))
print(modalirwin(2, "01", 0, 20))
print(modalirwin(10, "42", 3, 30))
print(modalirwin(37, [1, 2], 1, 20))
```

## Modal coefficients and signed terms

Request `modalcoeffs=True` to see modes in **increasing order `n=1,2,...`**.
This index is the same as the paper's `j` and equals `m+1` where the
paper indexes a series from `m=0`.

For `k >= 1`, each mode is printed as a separate, aligned block of fields,
one numerical value per line: `n`, `lambda`, `L`, `R`,
`C=L*R`, `C'=L*R/lambda`, and the actual signed contribution
`T=C*lambda**(k-1)` (equivalently `T=C'*lambda**k`).
The first mode has `lambda=1`, `L=log(b)`, `R=b**len(w)`.
**Any sign appearing in `L` is already included in `C`, `C'`, and `T`;**
there is no additional alternating-sign multiplier on the display.
A leading `+` is shown for positive terms.

For `k = 0`, only the numbered signed contributions `n` and `T` are printed
(one numerical value per line).
When `w=0^p`, the separate finite constant `sum(1/a, a=1,...,b-1)`
is printed before the modes. The degenerate binary forbidden word `w='1'`
with `k=0` has no admissible positive integers, so its sum is zero.

```python
modalirwin(2, "01", 2, 20, modalcoeffs=True)
modalirwin(2, "00", 0, 20, modalcoeffs=True)
```

The display uses values from the actual computation; it does not recalculate
the expensive modal integrals. Coefficients, per-mode diagnostics, and the
verbose summary each put one numerical value on each line, with aligned
values and spaces around the equal signs. `verbose` and `modalcoeffs` may
be combined. Each invocation computes one `k`; the returned value remains
a decimal string.

## Computational outline

For `k >= 1`, the code combines the right modal coefficients with the
left ones (computed from auxiliary Stieltjes functions) and eigenvalue
powers. It uses the eigenpolynomial change of basis and the triangular
recurrences for the moments of `eta_0` and the tilted probabilities `nu_j`.

For `k = 0`, it uses the modal Kempner expansion and additionally computes
the occurrence-free `mu_0` moments. The word `0^p` receives the further
level-raising transformation required in the Kempner paper.

Recurrence coefficients and moments are cached at the precision at
which they were computed. Working precision can decrease between modes,
while the accumulated sum is maintained at the requested digits plus
12 guard digits. When increased accuracy is needed, a term may be
recomputed using a fresh recurrence engine rather than attempting to
recover lost digits from a lower-precision cache.

**Accuracy and precision policy.** For requests of 30 significant digits or
fewer, the initial working precision is fixed at the requested digits plus
12 guard digits.  Above 30 digits, the budget for a new modal contribution
is reduced according to the magnitude of the preceding contribution
relative to `b**len(w) * log(b)`.  Positive and negative contributions are
accumulated separately in the Newton series and monitored for cancellation;
the moment and coefficient recurrences also contribute to the cancellation
indicator.  When the estimated loss of precision is too large, the current
mode is retried with increased precision and a fresh recurrence cache.
Higher precision may also be retained for subsequent modes when required.
The diagnostics `modalcoeffs=True` and `verbose=True` expose the precision
used and the number of retries.

The stopping threshold is `b**len(w) * log(b) * 10**(-(nboffigures+12))`.
The twelve digits here are a **working precision and stopping margin**,
not a mathematical guarantee of twelve correct guard digits after all
intermediate operations.
After two consecutive decreases in modal-term magnitude, the calculation
stops when twice the magnitude of the last term, together with its
cancellation uncertainty indicator, is below the threshold.  Twice the
last term is an **estimated bound for the tail starting with that term,
inclusive**, not an estimate of the remaining tail after that term.
The stopping indicator is not a mathematically certified error bound.
Reaching `max_modes` or the allowed precision-retry limit raises an
exception rather than returning an unchecked value. The exceptional
sum `I(2, "1", 0)` is exactly zero.

### mpmath version caveat: large powers of two

As described in [`README-method1.md`](README-method1.md#about-mpmath-versions),
the **pure-Python** backend in `mpmath` **1.4.0 and 1.4.1** is affected
by a slow conversion of Python integers with many trailing zero bits to
`mpf` ([mpmath issue 1073](https://github.com/mpmath/mpmath/issues/1073)).
That README reports the issue absent in 1.3.0 and when using `gmpy2`.
`blockirwin.py` has an explicit trailing-zero stripping workaround.

**`modalirwin.py` is potentially affected too.** In particular, its
recurrence coefficients involve `integer**degree` multiplied by
`mpmath` values (in `_ModalEngine._de()`), and level raising and moment
computations also involve large Python integers. For binary words and
powers-of-two radices, these operands can have many trailing zero bits.
The modal implementation has **no corresponding conversion workaround**;
we have not tested its performance under `mpmath` 1.4.0/1.4.1. For a
known unaffected configuration, use `mpmath` 1.3.0, or consult the
issue and `README-method1.md` when using a newer version or `gmpy2`.
This caveat concerns potential performance degradation, not evidence
of incorrect numerical results.

### Numerically expensive cases

The modal series are geometrically convergent, as established in Burnol's
papers.  **Geometric convergence does not imply that every instance is
fast at high requested precision.**  Some small-base examples, notably
Kempner sums (`k=0`) for binary words, require many modes before their
contributions become small enough.  Moreover, the work for each mode can
increase with the modal index: evaluating its coefficients and tilted
moments extends triangular recurrences and Newton series, and cancellation
can require additional precision or a fresh cache.  Thus two separate
costs matter: the **number of modes** and the **cost per mode**.

Heuristically, the coefficient of mode $m$ uses the moments of orders
$1,\dots,m$ of the tilted measures $\nu_m,\dots,\nu_1$, and a moment of
order $r$ comes from a recurrence over all lower orders.  With $M$ modes
this is about $M^2/2$ moments and $M^3/6$ operations, whatever the level
`ell`; in version 0.6.0 this took about two thirds of the time for binary
`11`, `k=1`, at 40 digits.  For $k\ge1$ mode $j$ carries the factor
$\lambda_j^{k-1}$, so fewer modes are needed as $k$ grows: for binary `01`
and `11` at 40 digits, $M(k+1)$ stayed within 3% for each word, for
$k=1,3,5,10$.
Hence the modal expansion is costly for small `k`, and fast for large
`k`, when a few modes beyond $b^p\log b$ suffice.

For example, the following are single-run measurements of a recent
working revision of `modalirwin.py`, in base 2, at **20 requested
significant digits** (CPython 3.13.5, `mpmath` 1.3.0). Times measure the function call,
excluding Python startup and imports; they are not performance guarantees.

| Word | `k` | Modal terms | Modal elapsed | `blockirwin.py` elapsed |
| --- | ---: | ---: | ---: | ---: |
| `01` | 0 | 54 | 0.58 s | not measured |
| `00` | 0 | 98 | 2.61 s | 0.003 s |
| `11` | 0 | 107 | 4.31 s | 0.003 s |
| `000` | 0 | 96 | 2.94 s | 0.012 s |
| `010` | 0 | 80 | 2.28 s | 0.006 s |
| `0110` | 0 | 94 | 4.09 s | not measured |
| `11` | 1 | 54 | 0.55 s | 0.004 s |
| `111` | 1 | 36 | 0.50 s | not measured |

All 20 combinations comprising binary `00`, `11`, and the eight three-digit
binary words, each for `k=0` and `k=1`, were also compared at 20 digits;
all completed with matching requested decimal strings.  The two methods
have very different computational costs, and the alternative method can be
much faster.  This does **not** establish that every case agrees at 40 or
60 digits.

In the final revision, binary `01` with `k=0` at 40 digits took
2.61 s and 87 modes (minimum 47 working digits, versus 52 initially); binary `00` and `11` with `k=0` each exceeded a
12-second per-case limit.  Those two cases were **not completed** at
40 digits in this measurement.

Additional independently checked measurements from earlier revisions
(same environment; use them as cost illustrations, not final-version benchmarks):

| Word | `k` | Digits | Modal terms | Modal time | `blockirwin.py` time | Agreement |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `001` | 0 | 40 | 59 | 1.08 s | 0.006 s | yes |
| `01` | 0 | 60 | 121 | 8.17 s | 0.003 s | yes |
| `11` | 1 | 60 | 120 | 5.72 s | 0.008 s | yes |
| `111` | 1 | 60 | 80 | 4.81 s | 0.021 s | yes |

A final-version check at 60 digits gave 7.31 s for binary `01`,
`k=0` (minimum 63 working digits versus 72 initially), and
4.75 s for binary `11`, `k=1` (72 working digits throughout
because its cancellation-sensitive recurrences prevented reduction).
Both matched `blockirwin.py`. These checks demonstrate that reducing
precision does not automatically reduce total running time.

More generally, these figures illustrate
implementation cost, not a universal ordering of convergence rates or a
claim that all Kempner cases are expensive.

For reproducible comparisons and independent numerical cross-checks, use
`compare_modalirwin.py` and report requested digits, number of computed
modes, both elapsed times, and whether the returned decimal strings agree.
The alternative `blockirwin.py` algorithm can be preferable for costly
modal cases.  A completed comparison at 20 digits does not establish
correctness at 40 or 60 digits.

### Decimal-radix performance and large occurrence counts

The dependence on `k` differs markedly between the two algorithms.
`blockirwin.py` does **not** use modal terms; its moment recurrences
can become more expensive as the occurrence count grows. The modal
expansion can instead converge after fewer modes for large `k`.
For instance, with `b=10`, `w="314"`, at 40 digits, earlier timings
were 0.0124 s and 0.0234 s for method 1 (`k=1` and `k=10`), versus
0.277 s and 0.0045 s for `modalirwin.py` 0.5.0. These measurements
illustrate contrasting trends, not a general monotonicity theorem.

With 0.6.0 and on another, slower, hardware (best of three runs,
October 12, 2026, CPython 3.13.16, `mpmath` 1.4.1), the same example
gave 0.026 s and 0.039 s for method 1 (`k=1` and `k=10`), versus
0.305 s and 0.0071 s for `modalirwin.py`; both methods returned the
same decimal strings.  On this hardware `blockirwin.py` ran about
1.7 to 2 times slower than in the earlier timings.

For ordinary decimal examples at **100 significant digits**, previous
single-run comparisons of `modalirwin.py` 0.3.0 and 0.5.0 gave:

| Word | `k` | 0.3.0 time | 0.5.0 time | Gain with 0.5.0 |
| --- | ---: | ---: | ---: | ---: |
| `42` | 3 | 1.133 s | 0.884 s | 22% |
| `123` | 1 | 1.914 s | 1.601 s | 16% |
| `7` | 2 | 3.505 s | 2.752 s | 21% |

All three agreed with `blockirwin.py` to the requested precision.
The modified solver's minimum working precision fell from its initial
112 digits to 24, 47, and 31 digits respectively. Savings are
**case-dependent**: for binary `01`, `k=0`, at 60 digits, the revised
code was about 21% slower than 0.3.0, despite reducing precision.
These are historical measurements, **not new benchmarks**; neither
algorithm is universally faster.

Version 0.6.0, single runs on October 12, 2026 (CPython 3.13.16,
`mpmath` 1.4.1, another machine than the tables above), all agreeing
with `blockirwin.py`; times in seconds, modes in parentheses:

| Base, word | Digits | `k=1` | `k=5` | `k=10` | `k=20` |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2, `11` | 40 | 3.59 (87) | 0.266 (29) | 0.163 (16) | 0.092 (9) |
| `blockirwin.py` | | 0.010 | 0.023 | 0.030 | 0.063 |
| 10, `42` | 40 | 0.498 (26) | 0.077 (6) | 0.016 (3) | 0.011 (3) |
| `blockirwin.py` | | 0.038 | 0.035 | 0.047 | 0.072 |
| 10, `42` | 100 | 3.06 (56) | 0.733 (12) | 0.367 (6) | 0.140 (3) |
| `blockirwin.py` | | 0.086 | 0.133 | 0.172 | 0.278 |

## Command line, comparisons, and tests

From this repository's root:

```sh
python modalirwin.py 2 01 0 20
python modalirwin.py 2 01 2 20 --modalcoeffs
python modalirwin.py 2 00 0 20 --modalcoeffs --verbose
python modalirwin.py 37 1,2 1 20
python compare_modalirwin.py 2 01 0 20
python -m unittest test_modalirwin.py test_modalirwin_cancellation.py
python -m pytest -q test_modalirwin.py test_modalirwin_cancellation.py test_blockirwin.py
python benchmark_modalirwin_binary.py --digits 20,40,60 --total-seconds 170
```

`benchmark_modalirwin_binary.py` tests binary `00`, `11`, and every
three-digit binary word for both `k=0` and `k=1`, at the selected
precisions.  Each independent method has an individual timeout; the driver
also has a total time limit and writes its CSV output after every case.
A timeout is recorded as *unverified*, not as success.

`compare_modalirwin.py` imports the two modules simultaneously, prints
each result, the elapsed wall-clock times, and the difference of their
printed values. The methods are numerically independent and use distinct
working algorithms. Timings should be compared at equivalent requested
precision; differences depend on hardware and Python/mpmath versions.

`test_modalirwin.py` checks the independent results against
`blockirwin.py` for positive `k`, Kempner `k=0`, one-digit and overlapping
blocks, words with leading zeros, the all-zero exception, explicit `ell`,
large bases, higher requested precisions, and the numeric truncation flag.
It also checks the recurrence row cache, coefficient-output identities,
mode numbering and signs, co-importing, the command-line program,
input validation, and cancellation-triggered precision retries.
`test_modalirwin_cancellation.py` checks signed-term accumulation,
12-digit guard settings, per-mode diagnostics, precision reduction on
an example, and the absence of wrapper-level global restarts.
The independent `test_blockirwin.py` suite does not use `modalirwin.py`.

For the other method, see [`README-method1.md`](README-method1.md).
