# -*- coding: utf-8 -*-
"""
=======================================================
Math helper functions (:mod:`sknano.core._extras`)
=======================================================

.. currentmodule:: sknano.core._extras

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext en'

from fractions import gcd
import operator

from ._vector import Vector

__all__ = ['e1', 'e2', 'e3', 'xhat', 'yhat', 'zhat',
           'comparison_symbol_operator_mappings',
           'math_symbol_operator_mappings',
           'symbol_operator_mappings',
           'totient_func']

e1 = xhat = Vector([1, 0, 0])
e2 = yhat = Vector([0, 1, 0])
e3 = zhat = Vector([0, 0, 1])

symbol_operator_mappings = {}

comparison_symbol_operator_mappings = \
    {'==': operator.eq,
     '!=': operator.ne,
     '<': operator.lt,
     '<=': operator.le,
     '>': operator.gt,
     '>=': operator.ge,
     'is': operator.is_,
     'is not': operator.is_not}
symbol_operator_mappings.update(comparison_symbol_operator_mappings)

math_symbol_operator_mappings = \
    {'+': operator.add,
     '-': operator.sub,
     '*': operator.mul,
     '/': operator.truediv}
symbol_operator_mappings.update(math_symbol_operator_mappings)


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
    for k in range(1, n + 1):
        if gcd(n, k) == 1:
            tots += 1

    return tots
