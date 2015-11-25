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

from pyparsing import Forward, Word, Literal, Suppress, Combine, Optional, \
    CaselessKeyword, CaselessLiteral, StringEnd, And, Or, ZeroOrMore, \
    ParseException, alphanums, alphas, nums, oneOf

from sknano.core import BaseClass
# from ._md_atoms import MDAtom as Atom, MDAtoms as Atoms

__all__ = ['SelectionException', 'SelectionParser']


class Selection(BaseClass):
    pass


class AllSelection(Selection):
    pass


class SelectionException(ParseException):
    pass


class SelectionParser(BaseClass):
    """Selection parser class."""

    ALL = CaselessKeyword('all')
    NONE = CaselessKeyword('none')
    NAME = CaselessKeyword('name')
    TYPE = CaselessKeyword('type')
    INDEX = CaselessKeyword('index')
    SERIAL = CaselessKeyword('serial')
    ATOMICNUMBER = CaselessKeyword('atomicnumber')
    ELEMENT = CaselessKeyword('element')
    RESIDUE = CaselessKeyword('residue')
    NUMBONDS = CaselessKeyword('numbonds')
    BOOLOP = oneOf('and or', caseless=True)
    BINOP = oneOf('< <= == > >= != LT LE EQ GT GE NE', caseless=True)
    LPAR = Suppress('(')
    RPAR = Suppress(')')

    ATOM_ATTRIBUTE = oneOf(' '.join(('x', 'y', 'z', 'r')))

    point = Literal('.')
    e = CaselessLiteral('E')
    number = Word(nums)
    plus_or_minus_sign = Word('+-', exact=1)
    integer = Combine(Optional(plus_or_minus_sign) + number)
    real_number = Combine(Optional(plus_or_minus_sign) +
                          (number + point + Optional(number) |
                           (point + number)) +
                          Optional(e) + Optional(plus_or_minus_sign) + number)

    NUMBERPATTERN = integer | real_number

    selection_pattern = Forward()
    selection_pattern << Optional(LPAR) + ATOM_ATTRIBUTE + BINOP + \
        NUMBERPATTERN + Optional(RPAR) + \
        ZeroOrMore(BOOLOP + selection_pattern)

    def __init__(self, atoms=None):
        super().__init__()
        self.atoms = atoms
        # self.ATOM_ATTRIBUTE = \
        #     oneOf(' '.join([attr for attr in dir(self.atoms) if not
        #                     attr.startswith('_')]))
        self.fmtstr = "atoms={atoms!r}"
        self.selstr = None

    def parse(self, selstr, **kwargs):
        try:
            return self.selection_pattern.parseString(selstr, parseAll=True)
        except ParseException as e:
            raise SelectionException(e.pstr, e.loc, e.msg, e.parseElement)

    def todict(self):
        return dict(atoms=self.atoms)
