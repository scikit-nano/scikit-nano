# -*- coding: utf-8 -*-
"""
===================================================
Core functions (:mod:`sknano.tools._corefuncs`)
===================================================

.. currentmodule:: sknano.tools._corefuncs

"""
from __future__ import division, print_function, absolute_import

__all__ = ['check_type']


def check_type(self, value, valid_type):
    """Check value type against `valid_type` tuple.

    Parameters
    ----------
    value : `object`
    valid_type : class or type or tuple of classes and/or types

    Raises
    ------
    `TypeError`
        If `value` fails `isinstance` check against `valid_type`.

    """
    if not isinstance(value, valid_type):
        raise TypeError('{} not valid type.\n'.format(value) +
                        '(Valid Type: {})'.format(valid_type))
