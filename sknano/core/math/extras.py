# -*- coding: utf-8 -*-
"""
=================================================================
Math helper functions (:mod:`sknano.core.math.extras`)
=================================================================

.. currentmodule:: sknano.core.math.extras

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from fractions import gcd
import operator

import numpy as np

__all__ = ['abs_cap', 'math_functions', 'function_map',
           'comparison_symbol_operator_mappings',
           'math_symbol_operator_mappings',
           'symbol_operator_mappings', 'operator_map',
           'math_operators', 'convert_condition_str',
           'totient_func']

math_functions = function_map = {}
function_map.update(dict(
    abs=np.abs,
    sqrt=np.sqrt,
    floor=np.floor,
    ceil=np.ceil,
    sin=np.sin,
    cos=np.cos,
    tan=np.tan,
    atan=np.arctan,
    asin=np.arcsin,
    acos=np.arccos,
    sinh=np.sinh,
    cosh=np.cosh,
    tanh=np.tanh,
    exp=np.exp,
    log=np.log,
    log10=np.log10))

math_operators = symbol_operator_mappings = operator_map = {}

comparison_symbol_operator_mappings = \
    {'==': operator.eq,
     '!=': operator.ne,
     '<': operator.lt,
     '<=': operator.le,
     '>': operator.gt,
     '>=': operator.ge,
     'is': operator.is_,
     'is not': operator.is_not}
operator_map.update(comparison_symbol_operator_mappings)

math_symbol_operator_mappings = \
    {'+': operator.add,
     '-': operator.sub,
     '*': operator.mul,
     '/': operator.truediv,
     '//': operator.floordiv}
operator_map.update(math_symbol_operator_mappings)


def abs_cap(val, max_abs_val=1):
    """Returns the value with its absolute value capped at max_abs_val.

    Modified implementation of
    :func:`pymatgen:pymatgen.util.num_utils.abs_cap`.

    Particularly useful in passing values to trignometric functions where
    numerical errors may result in an argument > 1 being passed in.

    Parameters
    ----------
    val : :class:`~python:float`
    max_abs_val : :class:`~python:float`
        The maximum absolute value for val. Default=1.

    Returns
    -------
    val if abs(val) < `max_abs_val` else sign of val * max_abs_val.
    """
    return max(min(val, max_abs_val), -max_abs_val)


def convert_condition_str(obj, condition):
    """Evaluates a boolean condtion of an object attribute.

    Parameters
    ----------
    obj : :class:`~python:object`
    condtion : :class:`~python:str`

    Returns
    -------
    :class:`~python:bool`

    Raises
    ------
    :class:`~python:ValueError`

    """
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
