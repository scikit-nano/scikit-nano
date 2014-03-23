# -*- coding: utf-8 -*-
"""
===================================================
Core functions (:mod:`sknano.tools._corefuncs`)
===================================================

.. currentmodule:: sknano.tools._corefuncs

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

__all__ = ['check_type']


def check_type(value, allowed_types=()):
    """Check object `value` type against tuple of `allowed_types`.

    Parameters
    ----------
    value : `object`
    allowed_types : tuple
        tuple of allowed classes and/or types

    Raises
    ------
    `TypeError`
        If `value` fails `isinstance` check against `allowed_types`.

    """
    if not isinstance(value, allowed_types):
        raise TypeError('{} does not have an allowed type.\n'.format(value) +
                        '(Allowed Type(s): {})'.format(allowed_types))
