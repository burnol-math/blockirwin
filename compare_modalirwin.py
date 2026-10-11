# -*- coding: utf-8 -*-
"""Compare outputs and call times of two independent numerical methods.

Produced by OpenAI GPT-6 on October 10, 2026, from J.-F. Burnol's
mathematics; see README-modalirwin.md and the CC BY-SA 4.0 LICENSE.
The blockirwin.py method was separately produced by Claude Opus 5.5
and Claude Sonnet 5.  Neither implementation imports the other.

Usage: python compare_modalirwin.py b w k [nboffigures [ell]]
       python compare_modalirwin.py 2 01 0 20
"""
from __future__ import annotations

import time
from mpmath import mp
from blockirwin import blockirwin
from modalirwin import modalirwin


def compare_irwin(b, w, k, nboffigures=30, ell=None, *, verbose=False):
    """Print both results and their separately measured call durations.

    Return decimal strings, elapsed seconds, and their absolute difference
    as a formatted string. An explicit ``ell`` is passed to both methods,
    although they use it for different auxiliary level-raising operations.
    Timings exclude imports and are not repeated benchmark medians.
    """
    started = time.perf_counter()
    classical = blockirwin(b, w, k, nboffigures, ell, verbose=verbose)
    middle = time.perf_counter()
    modal = modalirwin(b, w, k, nboffigures, ell, verbose=verbose)
    finished = time.perf_counter()
    with mp.workdps(nboffigures+10):
        difference = abs(mp.mpf(classical)-mp.mpf(modal))
        difference_text = mp.nstr(difference, 6)
    print('Method comparison:')
    fields = [
        ('base', b),
        ('word', repr(w)),
        ('occurrences', k),
        ('significant figures', nboffigures),
        ('blockirwin result', classical),
        ('blockirwin elapsed', f'{middle-started:.6f}s'),
        ('modalirwin result', modal),
        ('modalirwin elapsed', f'{finished-middle:.6f}s'),
        ('absolute difference', difference_text),
    ]
    width = max(len(label) for label, _ in fields)
    for label, value in fields:
        print(f'{label:<{width}} = {value}')
    print('Modal truncation is heuristic, not a certified error bound.')
    return {'blockirwin': classical, 'modalirwin': modal,
            'blockirwin_seconds': middle-started,
            'modalirwin_seconds': finished-middle,
            'absolute_difference': difference_text}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('b', type=int)
    parser.add_argument('w')
    parser.add_argument('k', type=int)
    parser.add_argument('nboffigures', type=int, nargs='?', default=30)
    parser.add_argument('ell', type=int, nargs='?', default=None)
    parser.add_argument('--verbose', action='store_true')
    ns = parser.parse_args()
    w = ns.w if ns.b <= 36 else [int(x) for x in ns.w.split(',')]
    compare_irwin(ns.b, w, ns.k, ns.nboffigures, ns.ell, verbose=ns.verbose)
