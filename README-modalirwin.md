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

This contribution and the repository are distributed under the existing
**Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)**
license. The full legal text is in [`LICENSE`](LICENSE); see also
[the CC BY-SA 4.0 license page](https://creativecommons.org/licenses/by-sa/4.0/).
Retain attribution and follow the license's share-alike requirements when
redistributing adaptations.

## Requirements and public API

Python 3 and [`mpmath`](https://mpmath.org/), plus Python's standard library.
No imports from `blockirwin.py` are needed inside the modal module.

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
| `ell` | Minimum number of base-`b` digits in the **auxiliary modal Stieltjes arguments** after level raising; default 3 in base 2 and 2 otherwise |
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

The expensive coefficient rows of the recurrences are cached per engine;
the cache and moment computations are confined to one working precision.
The Newton-series contributions use the same working precision,
with positive and negative contributions accumulated separately.

**Accuracy policy.** Every *individual modal term* is evaluated by
separately accumulating the positive and negative contributions of its
Newton series, including contributions that cancel when raised arguments
coincide. Cancellation inside the right and left modal coefficient
recurrences is also included in the precision indicator. Comparing the
gross contribution with the resulting difference estimates the number of
decimal digits lost by subtraction. A term which appears negligible is
not accepted on that basis unless the estimated uncertainty from this
cancellation is also negligible on the requested scale.

The first attempt has **12 decimal guard digits** beyond `nboffigures`.
If a term is not adequately resolved at that precision, the calculation
restarts with fresh caches and **12 additional decimal guard digits**:
12, 24, 36, 48, ... (maximum 240). The attempt is restarted because all
moment and modal-coefficient caches must share one working precision.
This check is applied separately to **each term**, including `k=0`.
The working-precision/cancellation estimate is heuristic, not a rigorous
interval bound on rounding errors.

The natural scale is `b**len(w) * log(b)`. The fixed significance
threshold is this scale multiplied by `10**(-(nboffigures+12))`,
regardless of the guard precision reached after retries. The series
stops when three successive modal terms decrease strictly in absolute
value (two decreases), provided twice the sum of the latest term's
magnitude and its estimated cancellation uncertainty is below the
threshold. This tail indicator is heuristic, not a rigorous error
bound. Auxiliary Newton series use a heuristic first-omitted-row
criterion. Reaching `max_modes` or the
precision ceiling raises an exception instead of reporting a result.
The exceptional exact sum `I(2, "1", 0)` returns `"0"` immediately.

The chosen method and precision policy are experimental; timings and the
number of modes depend on the block, occurrence count, base, and hardware.
For high-stakes precision requirements, independently cross-check results.
Some parameters, including longer-word Kempner cases, can be
computationally expensive. The current full regression suite did not
finish within a three-minute limit; it reached
`test_longer_words_and_stability` without reporting a failure. This is
not evidence that all remaining tests pass. `blockirwin.py` remains
available for independent cross-checks.

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
```

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
input validation, and the automatic working-precision safeguard.
`test_modalirwin_cancellation.py` provides short tests of signed-term
accumulation and linear guard increments. The independent
`test_blockirwin.py` suite is not modified.

For the other method, see [`README-method1.md`](README-method1.md).
