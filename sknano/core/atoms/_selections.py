# -*- coding: utf-8 -*-
"""
===============================================================================
Atom selection classes (:mod:`sknano.core.atoms._selections`)
===============================================================================

.. currentmodule:: sknano.core.atoms._selections

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# from operator import attrgetter

# import numpy as np

from pyparsing import Word, alphas, ParseException, Literal, \
    CaselessLiteral, Combine, Optional, nums, And, Or, Forward, ZeroOrMore, \
    StringEnd, alphanums

from sknano.core import BaseClass
# from ._md_atoms import MDAtom as Atom, MDAtoms as Atoms

__all__ = ['SelectionParser']

point = Literal('.')
e = CaselessLiteral('E')
number = Word(nums)
plus_or_minus_sign = Word('+-', exact=1)
real_number = Combine(Optional(plus_or_minus_sign) +
                      (number + point + Optional(number) | (point + number)) +
                      Optional(e) + Optional(plus_or_minus_sign) + number)

ident = Word(alphas, alphanums + '_')


class SelectionParser(BaseClass):
    """Selection parser class."""

    def __init__(self, selstr=None):
        super().__init__()
        self.selstr = selstr

    def todict(self):
        return dict(selstr=self.selstr)
