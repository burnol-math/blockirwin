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

1. **Available — full documentation: [README-method1.md](README-method1.md).**
   Script [`blockirwin.py`](blockirwin.py).  Proposition 11 of
   arXiv:2609.17045 expresses $I(b,w,k)$ through Stieltjes values of the
   measures $\eta_j$, Proposition 12 (level raising) pushes their arguments
   to at least $\ell$ digits, and these values are computed as alternating
   series in the moments of the $\eta_j$, obtained from the finite
   triangular recurrences of Proposition 13.  Any base, any block, any
   $k\ge0$.
2. **Planned:** the modal expansion, as in arXiv:2609.17045 ($k\ge1$) and
   arXiv:2609.31986 ($k=0$).
3. **Uncertain, may never be done:** modal expansions applied once a
   certain level has been reached.

The code is licensed under the
[Creative Commons Attribution-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-sa/4.0/).
