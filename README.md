# block-Irwin sums

This repository will gather several methods for computing to arbitrary
precision the *block-Irwin sums*

$$I(b,w,k)=\sum_{n>0,\ k_w(n)=k}\frac1n ,$$

where $b>1$ is the base, $w$ a block of base-$b$ digits and $k_w(n)$ the
number of occurrences of $w$ in the base-$b$ representation of $n$,
counted with overlaps.

References:

- J.-F. Burnol, *Block-count constrained harmonic sums: spectral expansion
  and block-directed Euler–Maclaurin*,
  [arXiv:2609.17045](https://arxiv.org/abs/2609.17045).
- J.-F. Burnol, *The modal expansion of Kempner sums*,
  [arXiv:2609.31986](https://arxiv.org/abs/2609.31986).

Methods:

1. Script [`blockirwin.py`](blockirwin.py) — full documentation: [README-method1.md](README-method1.md).
   Proposition 11 of
   arXiv:2609.17045 expresses $I(b,w,k)$ through Stieltjes values of the
   measures $\eta_j$, Proposition 12 (level raising) pushes their arguments
   to at least $\ell$ digits, and these values are computed as alternating
   series in the moments of the $\eta_j$, obtained from the finite
   triangular recurrences of Proposition 13.  Any base, any block, any
   $k\ge0$.
2. Script [`modalirwin.py`](modalirwin.py) - full documentation: [README-modalirwin.md](README-modalirwin.md).
   Implements the spectral/modal series of
   arXiv:2609.17045, Part 2 ($k\ge1$), and the modal Kempner series of
   arXiv:2609.31986 ($k=0$).
   It uses the published triangular recurrences for the moments of
   $\eta_0$, $\mu_0$, and the tilted measures $\nu_j$, and the
   corresponding modal coefficients, with level-raised Newton series
   for the auxiliary Stieltjes integrals.
3. *Uncertain, may never be done:* modal expansions applied once a
   certain level has been reached.


## Usage and validation

```python
from blockirwin import blockirwin
from modalirwin import modalirwin

print(blockirwin(2, "01", 0, 20))
print(modalirwin(2, "01", 0, 20))
modalirwin(2, "01", 2, 20, modalcoeffs=True)
```

```sh
python compare_modalirwin.py 2 01 0 20
python -m pytest -q test_modalirwin.py test_blockirwin.py
```

The modules only use Python + `mpmath`

Regarding `modalirwin.py`, the numeric stopping policy, including a crude
estimate for the first omitted modal term and reduced precision for small
Newton rows, is heuristic, not a rigorous accuracy certificate.  The modal
method may be slower on some parameters.  `max_modes` raises an error rather
than quietly claiming accuracy when convergence is too slow. A detected
numerical instability triggers fresh-cache retries at higher precision.


## AI assistance

The initial `blockirwin.py`, `test_blockirwin.py`, `README-method1.md`,
and this `README.md`'s original version were produced on September 28–29, 2026
by Claude Opus 5.5 and Claude Sonnet 5 (Anthropic), from the mathematics
and earlier Maple/SageMath/mpmath code of Jean-François Burnol.

The complementary `modalirwin.py`, `test_modalirwin.py` and
`compare_modalirwin.py` were produced with OpenAI ChatGPT (GPT-6) on October
10, 2026, from the two Burnol preprints.  This `README.md` was extended to
cover this, and then edited by a human, to prune some too lengthy ChatGPT
output.

The 
mathematics is fully authored by Jean-François Burnol. It was parsed
by ChatGPT 5.6 Sol to search for errors (and improve the English language use).

The repository is licensed
under the [Creative Commons Attribution-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-sa/4.0/).

