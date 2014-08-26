# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class with extended feature set (:mod:`sknano.core.atoms._bonds`)
===============================================================================

An "eXtended" `Atom` class for structure analysis.

.. currentmodule:: sknano.core.atoms._bonds

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numpy as np
from sknano.core.math import Vector

__all__ = ['Bond', 'Bonds', 'BondMixin']


class Bond(object):
    """Abstract data structure for `Atom` bonds.

    Parameters
    ----------
    """

    def __init__(self, **kwargs):
        pass

    def __repr__(self):
        """Return string representation of `Bond`."""
        return "Bond()"


class Bonds(object):
    pass


class BondMixin(object):
    pass
