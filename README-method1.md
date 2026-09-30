# blockirwin

[[_TOC_]]

## Description

[`blockirwin.py`](blockirwin.py) is a pure Python +
[mpmath](https://mpmath.org) script computing to arbitrary precision the
*block-Irwin sums*

$$I(b,w,k)=\sum_{n>0,\ k_w(n)=k}\frac1n ,$$

where $b>1$ is the base, $w$ is a non-empty block (word) of base-$b$ digits,
possibly with leading zeros, and $k_w(n)$ is the number of occurrences of $w$
in the minimal base-$b$ representation of $n$, occurrences being counted with
overlaps.

The implementation is closely modeled on
[`kempner.py`](https://gitlab.com/burnolmath/kempner), which was authored by
the human Jean-François Burnol in April 2026, and computes
analogous harmonic sums, where the condition on the denominators is to use
only some allowed digits.  The mathematics for the block
case (only one block, so far) was provided in Part 2 of

> J.-F. Burnol, *Block-count constrained harmonic sums: spectral expansion and
> block-directed Euler–Maclaurin*,
> [arXiv:2609.17045](https://arxiv.org/abs/2609.17045),

namely:

1. **Proposition 11** expresses $I(b,w,k)$ as a finite integer combination of
   values $V_j(N)=\int_{[0,1)}\frac{d\eta_j(x)}{N+x}$, $j\in\{k,k-1\}$, of the
   Stieltjes functions of the measures $\eta_j$.
2. **Proposition 12** (level raising) replaces $V_j(N)$ by values at
   arguments at least $bN$, of indices $j$ and $j-1$, plus $1/N$ if $j=0$.
   It is applied until all arguments have at least $\ell$ digits in base
   $b$ ($\ell$ is the *level*).
3. **Section 30** gives finite triangular recurrences (Proposition 13) for
   the moments $v_{j;m}$ of the measures $\eta_j$, $0\le j\le k$.  The
   $V_j(N)$ are then evaluated as the alternating series
   $\sum_m(-1)^m v_{j;m}N^{-m-1}$, which converge geometrically with ratio
   at most $b^{1-\ell}$.

All cases are handled: any base, any block (bordered or not, with or
without leading zeros), any $k\geq0$.  The modal expansion, which is the
main topic of Part 2 of the paper and, for $k=0$, of its continuation

> J.-F. Burnol, *The modal expansion of Kempner sums*,
> [arXiv:2609.31986](https://arxiv.org/abs/2609.31986),

is not used and will be provided in some other code in this repository.

The moment recurrence of Proposition 13 is implemented directly for the
moments $v_{j;m}$, all indices $j\le k$ at once.  All quantities
are floating-point numbers at the running precision, which decreases with
$m$ as in [`kempner.py`](https://gitlab.com/burnolmath/kempner), since the
moments of index $m$ only enter through terms of size at most
$b^p/N^{m+1}$.  The power sums $\sum\pm x(s)^r$ over the points
$x(s)=n(s)/b^{|s|}$ of the words entering $E_w$ are updated
incrementally, word pairs which cancel exactly are removed beforehand,
and the binomial coefficients and the powers of $b$ are obtained on the
fly.  The numbers kept from earlier steps are rounded down to the running
precision whenever it has decreased by 10%, so that multiplications only
involve mantissas of about the running precision.

In the spirit of the [kempner](https://gitlab.com/burnolmath/kempner) script,
the code is short and no attempt at parallelization is made.  For single-digit
$w$, see also the [irwin](https://gitlab.com/burnolmath/irwin) repository
whose `irwin_v5` (SageMath by my human collaborator or its Python+mpmath version
I produced under its wise direction) has faster runtimes, but this is
depending on the hardware, and is up to the price of a far more complex internal
code structure.

## Usage

Interactively:

```python
>>> from blockirwin import blockirwin
>>> blockirwin(10, "9", 1, 50)
>>> blockirwin(10, "42", 1, 40)
>>> blockirwin(2, "11", 0, 40)
>>> blockirwin(16, "abab", 2, 30)
```

The signature is

```python
blockirwin(b, w, k, nboffigures=100, ell=None, verbose=False, trunc=False)
```

- `b`: the base, any integer at least 2.
- `w`: the block, a string of base-$b$ digits `0-9a-z` (only for $b\le36$),
  or a list of integers, which works for any base: for example
  `blockirwin(100, [4, 2], 1, 50)` for the block made of the base-100
  digits 4 and 2.  Leading zeros are significant: `"09"` is not `"9"`.
- `k`: the number of occurrences, an integer at least 0.
- `nboffigures`: the number of significant figures of the returned string.
- `ell`: the level, at least 2.  Defaults to the rule of thumb of
  `default_ell(b, w, k, nboffigures)`, see below.
- `verbose`: if `True`, prints details: the borders of $w$, the Stieltjes
  arguments $N$ (both in base $b$ and in decimal; for $b>36$ the base-$b$ form is
  a list of digits such as `[1,4,2]`) with the coefficients of
  the $V_j(N)$ and the number of series terms kept, the exact terms $1/n$,
  the numbers of moments of the measures $\eta_j$ computed and used, and
  timings.  Long listings are abridged; `verbose=2` prints them in full.
- `trunc`: truncate rather than round.

The return value is a string.

From the command line:

```text
python blockirwin.py b w k [nboffigures [ell]]
```

for example `python blockirwin.py 10 314159 0 60`.  For $b>36$ the block
is given as comma-separated digits, for example
`python blockirwin.py 100 4,2 1 50`.

## Choice of the level

The level $\ell$ balances two costs.  Raising it by one reduces the number
of moments needed (about $P/((\ell-1)\log_2 b)$ for $P$ bits of
precision), but multiplies the number of Stieltjes arguments by up to $b$.
When `ell` is not given, `default_ell(b, w, k, nboffigures)` applies the
following rule, where $d$ is the number of requested digits and $p=|w|$:

| base | case | level |
|---|---|---|
| 2 | $p=1$, $k=0$ | $0.4d+10$ |
| 2 | $p=1$, $k=1$ | $15+0.05d$ |
| 2 | $p=1$, $k=2$ | $13+0.028d$ |
| 2 | $p=1$, $k\ge3$ | 10 up to 150 digits, 12 up to 450, then 13 |
| 2 | $p\ge2$, $k=0$ | 11 up to 150 digits, 12 up to 300, then 14 |
| 2 | $p\ge2$, $k\ge1$ | 10 up to 150 digits, 11 up to 300, then 12 |
| 3 | $p=1$ | 7 up to 150 digits, then 8 |
| 3 | $p\ge2$ | 6 up to 300 digits, then 7 |
| 4 | any | 5 up to 300 digits, then 6 |
| 5 | any | 4 up to 150 digits, then 5 |
| 10 | any | 3 below 400 digits; from 600 digits 4 if $k\ge1$; in between 4 only if $k\ge2$; for $k=0$, 3 up to 600 digits, then 4 |
| 16 | any | 3 |
| other | any | $1+\operatorname{round}(\log_b(\operatorname{nboffigures}+100))$ |

This rule was obtained from timings with pure-Python mpmath 1.4.1 (no
`gmpy2`), for bases 2, 3, 4, 5, 10 and 16, precisions from 100 to 800
digits, $0\le k\le10$ and blocks of length 1 to 5.  On these 273 measured
cases, the time at the level chosen by the rule is within 10% of the time
at the best level in 251 cases, and within 29% in all cases.  Above 800
digits the values are extrapolated, not measured.  For other bases the
formula reproduces the measured bases 3 to 16 within one level, but it
has not been measured for them.

In base 2 with a one-digit block, few binary strings contain at most $k$
occurrences of that digit: the number of Stieltjes arguments grows only
polynomially with the level, so that high levels pay off, the more so as
$k$ is small.

Passing `ell` explicitly always overrides the rule.

## Validation

Running [`test_blockirwin.py`](test_blockirwin.py) checks:

1. values from the literature, computed by other methods: Schmelzer–Baillie
   (no `314159`), Baillie ([arXiv:0806.4410](https://arxiv.org/abs/0806.4410):
   no `314`, no `09`, the Erdős–Borwein constant $I(2,0,0)=I(2,10,0)$),
   the "exactly one" sums of [arXiv:2402.14761](https://arxiv.org/abs/2402.14761)
   (`42`, `35`, `00`, `99` in base 10, `00`, `11`, `10`, `01` in base 2),
   and $I(10,9,0)$ from the [irwin](https://gitlab.com/burnolmath/irwin)
   README;
2. values computed by the author in February 2024 by other means, for
   blocks of length 2 ($k=0$ and $k=1$) and of length 3 ($k=0$), covering
   all border structures;
3. that the result does not depend on the level $\ell$, for bordered
   blocks in bases 2, 3, 10, 16;
4. that $\bigl(I(b,w,k+1)-b^p\log b\bigr)/\bigl(I(b,w,k)-b^p\log b\bigr)$
   tends to $\lambda_2(b,w)=G(b^{-2})$, where $G=E/D$ is the generating
   function of the return language (Theorem 2 of arXiv:2609.17045; for
   $k=0$ see its continuation arXiv:2609.31986) — a check the algorithm's
   own arithmetic never uses.

All these tests pass with the current version.  During development, the
single-digit cases were moreover cross-checked against `kempner.py`
($k=0$) and against the Python conversion of `irwin_v5.sage` ($k\geq1$),
in bases 2, 3, 5, 7, 10, 16.

## About mpmath versions

mpmath 1.4.0 and 1.4.1, when used with their pure-Python backend (i.e.
without `gmpy2`), convert to `mpf` the Python integers having many trailing
zero bits (such as large powers of 2) in quadratic time, see
[mpmath issue 1073](https://github.com/mpmath/mpmath/issues/1073),
which was opened by the human collaborator some months back.  This is
fixed in the development version of mpmath but, as of September 28, 2026,
not in a release; mpmath 1.3.0 is not affected, nor is any version used
with `gmpy2`.  Such integers occur in this computation, especially in bases
2 and 16, so the script strips trailing zero bits before every conversion
(`_mpf_from_int`).  Installing [`gmpy2`](https://pypi.org/project/gmpy2/)
is recommended in any case, for speed.

## License

The code was produced on September 28-29, 2026 by Claude Opus 5.5 and Claude
Sonnet 5 (Anthropic) from the formulas of the above paper by Jean-François
Burnol.  It was modeled on the one-digit case handled by the human in the
[kempner](https://gitlab.com/burnolmath/kempner) repository. It is licensed
under the
[Creative Commons Attribution-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-sa/4.0/).
