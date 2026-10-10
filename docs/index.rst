=====================
burnolmath/blockirwin
=====================

.. note::

   This entire page was authored by OpenAI ChatGPT (GPT-6) on
   October 10, 2026.  See `AI assistance`_ at the bottom
   for details.

This page is a short mathematical introduction to the Python + ``mpmath``
software at https://gitlab.com/burnolmath/blockirwin.  It concerns the
*block-Irwin sums* studied by Jean-François Burnol.  The repository supplies
two complementary numerical implementations; the mathematical developments
summarized here are Burnol's work.  Consult the papers linked below
for proofs and complete statements.

.. contents::

The sums and their background
=============================

Fix an integer base $b\geq 2$, a nonempty word $w$ of $p$ base-$b$ digits,
and an integer $k\geq0$.  Write $k_w(n)$ for the number of occurrences of
$w$ in the usual base-$b$ representation of a positive integer $n$.
Occurrences may overlap.  The block-Irwin sum is

.. math::
   :label: blockirwin_def

   I(b,w,k)=\sum_{n\geq1,\;k_w(n)=k}\frac{1}{n}.

The case $k=0$ is a *block-Kempner sum*: the word $w$ is forbidden.
When $w$ consists of one digit, these sums go back to the classical
questions of Kempner (1914) and Irwin (1916).

For multi-digit words and $k=0$, Schmelzer and Baillie (2008)
provided a numerical algorithm for these Kempner sums; see also
Baillie's 2008 arXiv preprint documenting further Mathematica code.
These works did not treat $k>0$ for
multi-digit words.  Some questions left open in
the 2008 Schmelzer--Baillie work were solved recently by Schmelzer
(private communication to Burnol, Sep. 17, 2026).

Allouche, Hu and Morin proved
that, for a word of length $p$, the sums tend to $b^p\log b$ as
$k\to\infty$.  Independently, Burnol's work in *Acta Mathematica
Hungarica* introduced the associated measures $\mu_k$ and proved their
weak convergence to $b^p$ times Lebesgue measure, also yielding this limit.
In an earlier paper in *Advances in Applied Mathematics* (2025), Burnol
had already obtained numerical summation formulas for two-digit words,
including the cases $k=0$ (already covered by Schmelzer--Baillie) and $k=1$.

The 2026 work develops two distinct but intertwined themes.

Numerical evaluation based on triangular recurrences of moments of measures
===========================================================================

Burnol extends the numerical framework of his earlier Kempner and Irwin
papers to every nonempty word $w$ and every $k\geq0$, including
self-overlapping words.  The method expresses the harmonic sums through
auxiliary measures whose moments are obtained from *triangular
recurrences*.  A technique called *level raising* leads to convergent
numerical series that can be made rapidly geometrically convergent.

This method is implemented in ``blockirwin.py``.  It is applicable to both
$k=0$ and $k\geq1$.

Modal expansions
================

Burnol also discovered how these sums
depend on the occurrence count.  Explicit computations of moments as $k$
varies suggested a decomposition into modes; the resulting theory uses
return operators, eigenvalues and eigenpolynomials.  For $k\geq1$, it gives
an absolutely convergent expansion

.. math::
   :label: blockirwin_modal

   I(b,w,k)=b^p\log b
      +\sum_{j=2}^{\infty} C'_j(b,w)\,\lambda_j(b,w)^k,
   \qquad 1>\lambda_2>\lambda_3>\cdots>0.

The notation $C'_j$ matches the ``C'`` coefficient label in the Python
scripts' user interface; it is not notation used in Burnol's papers.

Here the eigenvalues and the signed coefficients $C'_j$ depend on
$b$ and $w$, but *not* on $k$.  Thus, once the modal coefficients are
computed, evaluating the sum for other positive occurrence counts requires
only the powers of the eigenvalues and their weighted sum.  The leading
term $b^p\log b$ also makes the earlier limit transparent.

The case $k=0$ has a separate modal expansion, developed in Burnol's
companion 2026 paper on Kempner sums.  Words made entirely of zeros require
a special treatment, including an additional level-raising step.  Such
level raising can also be applied in general and the modal expansion
can start at a higher level (this is briefly hinted at in
Burnol's papers).

A binary example: the word ``11``
=================================

In this section we examine an example with $b=2$.
The word ``11`` overlaps with itself.  The binary integer
``111`` (that is, $7$) contains two occurrences of ``11``: in its
first two positions and in its last two positions.  Hence $1/7$
contributes to $I(2,11,2)$, not to $I(2,11,1)$.

This example is not new as an individual two-digit problem: Burnol's
2025 paper already covered ``11`` for $k=0$ and $k=1$.  What the 2026
papers add is the general numerical treatment of *all* $k$, and a modal
expansion whose principal term is $4\log 2$ and whose remaining terms
involve powers of the eigenvalues.  The coefficients obtained for one
occurrence count can be reused for other positive counts.

The following values were computed to 40 significant decimal digits.
The two implementations give identical displayed results for these cases.
The values approach $4\log 2$ as $k$ increases.

.. list-table:: Binary word ``11``: values of $I(2,11,k)$
   :header-rows: 1
   :widths: 8 92

   * - $k$
     - $I(2,11,k)$ (40 significant digits)
   * - 0
     - ``3.704711752910469457886531055976801955909``
   * - 1
     - ``2.938413407650189451517403882901748505185``
   * - 2
     - ``2.811456728628971939664929206845353589612``
   * - 3
     - ``2.782542836166987025295658250418892939816``
   * - 4
     - ``2.775229126541100652913024137035591439708``
   * - 5
     - ``2.773299809564120216836783887022996209558``
   * - 6
     - ``2.772781532882320027829114765871723260160``
   * - 7
     - ``2.772641165720707274198827749911715332619``
   * - 8
     - ``2.772603007130006482705845856202270154514``
   * - 9
     - ``2.772592615849791558660113071680734422499``
   * - 10
     - ``2.772589783845005247372727781720571666760``
   * - $4\log 2$ (limit)
     - ``2.772588722239781237668928485832706272302``

The modal data make the dependence on $k$ especially explicit.  For the
binary word ``11``, the eigenvalues in the positive-$k$ expansion are

.. math::

   \lambda_n(2,11)=\frac{2^n-1}{2^{2n}-2^n-1},\qquad n\geq1.

The first eight eigenvalues and their corresponding coefficients $C'_n$
are listed below.  The coefficients are signed quantities (although all
eight displayed values happen to be positive); the index $n$ starts at
$1$, with $\lambda_1=1$ and $C'_1=4\log 2$.  The decimals in this table
are rounded to 25 significant digits.

.. list-table:: First eight modes for $I(2,11,k)$, $k\geq1$
   :header-rows: 1
   :widths: 7 19 74

   * - $n$
     - $\lambda_n$ (exact)
     - $C'_n$ (25 significant digits)
   * - 1
     - $1$
     - ``2.772588722239781237668928``
   * - 2
     - $3/11$
     - ``0.4662020590105378139526108``
   * - 3
     - $7/55$
     - ``0.2260007777178610723648151``
   * - 4
     - $15/239$
     - ``0.1175876573062431455082605``
   * - 5
     - $31/991$
     - ``0.06051075912618450960831403``
   * - 6
     - $63/4031$
     - ``0.03074152760268022734100137``
   * - 7
     - $127/16255$
     - ``0.01549748960785583823063004``
   * - 8
     - $255/65279$
     - ``0.007780754329337076506497514``

For instance, at $k=10$, the contributions from modes $n=2$ and $n=3$
are about $1.06\times10^{-6}$ and $2.52\times10^{-10}$, respectively.
At $k=100$, the contribution from mode $n=2$ is only about
$1.74\times10^{-57}$.  Once the eigenvalues and coefficients have been
calculated to suitable precision, evaluating the sum for larger positive
$k$ involves no new moment recurrences or auxiliary series computations:
one simply forms the powers $\lambda_n^k$ and sums the contributions.
Usually fewer modes are numerically significant as $k$ increases.
The eight rows above illustrate the modal data; they are not the full
truncation used for the 40-digit values of $I(2,11,k)$ in the preceding
table.  The $k=0$ value uses the separate Kempner modal expansion.

*All decimal values in the two binary tables are rounded to the stated
number of significant digits; the rational eigenvalues are exact.*

A decimal example: the word ``314``
===================================

The word ``314`` in base $10$ has no proper self-overlap.  Its first five
positive-$k$ modal eigenvalues and coefficients $C'_n$ are as follows.
The first mode has $\lambda_1=1$ and $C'_1=1000\log 10$.  In this
example the nonprincipal eigenvalues decrease very rapidly, although
the coefficients $C'_n$ themselves are of comparable size.

.. list-table:: First five modes for $I(10,314,k)$, $k\geq1$
   :header-rows: 1
   :widths: 7 29 64

   * - $n$
     - $\lambda_n$ (exact)
     - $C'_n$ (25 significant digits)
   * - 1
     - $1$
     - ``2302.585092994045684017991``
   * - 2
     - $1/900001$
     - ``3.534460625252272365389384``
   * - 3
     - $1/990000001$
     - ``3.239359263987618323927269``
   * - 4
     - $1/999000000001$
     - ``3.193498415158547982752342``
   * - 5
     - $1/999900000000001$
     - ``3.184205350809411216123445``

The following values were computed independently with ``blockirwin.py``
and ``modalirwin.py``.  Both implementations agree at all displayed digits.
As in the binary example, we use 40 *significant decimal digits* (not 40
places after the decimal point).  Since these values have four digits
before the decimal point, 36 decimal places are displayed.

.. list-table:: Decimal word ``314``: values of $I(10,314,k)$
   :header-rows: 1
   :widths: 8 92

   * - $k$
     - $I(10,314,k)$ (40 significant digits)
   * - 0
     - ``2299.829782767518338453586353611974367846``
   * - 1
     - ``2302.585096924495072937876813819397855219``
   * - 2
     - ``2302.585092994050047543235982364566287362``
   * - 3
     - ``2302.585092994045684022839807011081374343``
   * - 4
     - ``2302.585092994045684017991460071416359078``
   * - 5
     - ``2302.585092994045684017991454684370193208``
   * - 6
     - ``2302.585092994045684017991454684364207608``
   * - 7
     - ``2302.585092994045684017991454684364207601``
   * - 8
     - ``2302.585092994045684017991454684364207601``
   * - 9
     - ``2302.585092994045684017991454684364207601``
   * - 10
     - ``2302.585092994045684017991454684364207601``
   * - $1000\log 10$ (limit)
     - ``2302.585092994045684017991454684364207601``

The entries for $k\geq7$ coincide with the limit at this displayed
precision: their differences from $1000\log 10$ are smaller than the
rounding threshold.  In this example the rapid decay of the
nonprincipal eigenvalues makes the reuse of modal data for larger $k$
particularly effective.

*All decimal values in these two decimal tables are rounded to the stated
number of significant digits (25 for the modal coefficients, 40 for the
sums and their limit); the rational eigenvalues are exact.*

Software
========

The two Python + ``mpmath`` scripts are independent and can be imported
together.  For example::

    from blockirwin import blockirwin
    from modalirwin import modalirwin

    print(blockirwin(2, "11", 2, 20))
    print(modalirwin(2, "11", 2, 20))

``blockirwin.py`` implements the numerical summation method based on
triangular recurrences.  ``modalirwin.py`` implements the modal expansions;
its optional ``modalcoeffs=True`` output shows the individual eigenvalues,
modal coefficients and signed contributions for $k\geq1$ (and the signed
terms for $k=0$).

See the repository's ``README-method1.md`` and ``README-modalirwin.md`` for
instructions, testing and numerical limitations.  In particular, the
current modal implementation uses a *heuristic*, not certified, numerical
stopping criterion.  The repository carries the `Creative Commons
Attribution-ShareAlike 4.0 International license
<https://creativecommons.org/licenses/by-sa/4.0/>`__ (CC BY-SA 4.0).

References
==========

The 2026 papers

- Jean-François Burnol, *Block-count constrained harmonic sums: spectral
  expansion and block-directed Euler--Maclaurin*,
  `arXiv:2609.17045 <https://arxiv.org/abs/2609.17045>`__.
  In particular, Part 2 presents the general modal theory for $k\geq1$
  as well as the triangular-recurrence numerical method.

- Jean-François Burnol, *The modal expansion of Kempner sums*,
  `arXiv:2609.31986 <https://arxiv.org/abs/2609.31986>`__.
  This develops the $k=0$ case.

Partial list of earlier works

- Thomas Schmelzer and Robert Baillie, *Summing a curious, slowly
  convergent series*, *American Mathematical Monthly* 115 (6)
  (2008), 525--540.
  `DOI <https://doi.org/10.1080/00029890.2008.11920559>`__.

- Robert Baillie, *Summing the curious series of Kempner and Irwin*
  (2008), `arXiv:0806.4410 <https://arxiv.org/abs/0806.4410>`__.

- Thomas Schmelzer, *manuscript*, 19 pages, September 17, 2026
  (communicated privately to Jean-François Burnol).

- J.-P. Allouche, Y. Hu and C. Morin, *Ellipsephic harmonic series
  revisited*, *Acta Mathematica Hungarica* 173 (2024), 461--470.
  `DOI <https://doi.org/10.1007/s10474-024-01448-5>`__.

- Jean-François Burnol, *Summing the "exactly one 42" and similar subsums
  of the harmonic series*, *Advances in Applied Mathematics* 162
  (January 2025), 102791. `DOI <https://doi.org/10.1016/j.aam.2024.102791>`__.

- Jean-François Burnol, *Measures associated with certain ellipsephic
  harmonic series and the Allouche-Hu-Morin limit theorem*,
  *Acta Mathematica Hungarica* (2025).
  `DOI <https://doi.org/10.1007/s10474-025-01525-3>`__.


AI assistance
=============

All the mathematics in Burnol's papers, including the
numerical schemes and modal expansions, is entirely his work.

The Python implementations are AI-assisted realizations:
``blockirwin.py`` was produced with Anthropic's Claude on
September 28--29, 2026, while ``modalirwin.py`` was produced with
OpenAI ChatGPT (GPT-6) on October 10, 2026.
Anthropic's Claude code took as starting point the earlier
Burnol authored implementations for single digit case available at
https://gitlab.com/burnolmath/kempner and
at https://gitlab.com/burnolmath/irwin.

This entire introductory
web page was authored by OpenAI ChatGPT (GPT-6) on October 10, 2026.
It is modeled after https://burnolmath.gitlab.io/irwin, and was
minimally edited by a human after production.
