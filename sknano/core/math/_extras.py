# -*- coding: utf-8 -*-
"""
=======================================================
Math helper functions (:mod:`sknano.core._extras`)
=======================================================

.. currentmodule:: sknano.core._extras

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from fractions import gcd
import operator

__all__ = ['comparison_symbol_operator_mappings',
           'math_symbol_operator_mappings',
           'symbol_operator_mappings',
           'math_operators', 'convert_condition_str',
           'totient_func']

math_operators = symbol_operator_mappings = {}

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
     '/': operator.truediv,
     '//': operator.floordiv}
symbol_operator_mappings.update(math_symbol_operator_mappings)


def convert_condition_str(obj, condition):
    try:
        attr, comparison, value = condition.split()
    except ValueError:
        print('the condition string must be of the form:\n'
              'ATTRIBUTE COMPARISON_OPERATOR VALUE')
    condition = \
        math_operators[comparison](getattr(obj, attr), float(value))

    return condition


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
