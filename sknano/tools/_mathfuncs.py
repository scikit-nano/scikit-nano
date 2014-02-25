# -*- coding: utf-8 -*-
"""
=======================================================
Math helper functions (:mod:`sknano.tools._mathfuncs`)
=======================================================

.. currentmodule:: sknano.tools._mathfuncs

"""
from __future__ import division, print_function, absolute_import

from fractions import gcd

__all__ = ['totient_func']


def totient_func(n=int):
    """Compute the totatives of :math:`n`.

    Parameters
    ----------
    n : int

    Returns
    -------
    int : the number of totatives of :math:`n`

    """
    if not isinstance(n, int) and int(n) != n:
        raise ValueError('n must be an integer')
    n = int(n)

    tots = 0
    for k in xrange(1, n + 1):
        if gcd(n, k) == 1:
            tots += 1

    return tots
