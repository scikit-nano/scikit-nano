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

from abc import abstractmethod

from pyparsing import Forward, Word, Literal, Suppress, Combine, Optional, \
    Group, CaselessKeyword, CaselessLiteral, StringEnd, And, Or, ZeroOrMore, \
    OneOrMore, ParseException, delimitedList, infixNotation, alphanums, \
    alphas, nums, oneOf

from sknano.core import BaseClass
from sknano.core.math import operator_map
# from ._md_atoms import MDAtom as Atom, MDAtoms as Atoms

__all__ = ['SelectionException', 'SelectionParser']


# def selection_wrapper(cls):
#     def selection(s, loc, tokens)

class Selection(BaseClass):
    def __init__(self):
        super().__init__()

    @abstractmethod
    def apply(self):
        return NotImplementedError

    def todict(self, selection):
        return dict(selection=self.selection)


class AndOperator(Selection):

    def __init__(self, left, right):
        self.left = left
        self.right = right

    def apply(self, atoms):
        return self.left.apply(atoms) & self.right.apply(atoms)


class AllSelection(Selection):
    pass


class SelectionException(ParseException):
    pass


class IDSelection(Selection):
    def __init__(self, selection):
        super().__init__()
        self.selection = selection

    def apply(self, atoms):
        return atoms.filtered_ids(self.selection[0].asList())


class AttributeSelection(Selection):
    def __init__(self, selection):
        super().__init__()
        self.selection = selection
        print(selection)

    def apply(self, atoms):
        attr, op, val = self.selection
        return atoms.filtered_ids(op(getattr(atoms, attr), val))


def operator_function(s, l, t):
    return operator_map[t[0]]


def asint(s, l, t):
    return int(t[0])


def asfloat(s, l, t):
    return float(t[0])


class SelectionParser(BaseClass):
    """Selection parser class."""
    ALL = CaselessKeyword('all')
    NONE = CaselessKeyword('none')
    NAME = CaselessKeyword('name')
    TYPE = CaselessKeyword('type')
    INDEX = CaselessKeyword('index')
    ID = CaselessKeyword('id')
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
    integer = Combine(Optional(plus_or_minus_sign) + number) \
        .setParseAction(asint)
    positive_integer = number.setParseAction(asint)
    real_number = \
        Combine(Optional(plus_or_minus_sign) +
                (number + point + Optional(number) | (point + number)) +
                Optional(e) + Optional(plus_or_minus_sign) + number) \
        .setParseAction(lambda s, l, t: float(t[0]))

    NUMBERPATTERN = integer | real_number

    selection_pattern = Forward()

    # select_all_pattern = Forward()

    # boolean_expression =

    id_selection_expression = \
        (Optional(LPAR) + Suppress(ID) +
         Group(delimitedList(OneOrMore(positive_integer), delim=' ')) +
         Optional(RPAR)).setParseAction(IDSelection)

    attr_selection_expression = \
        (Optional(LPAR) +
         ATOM_ATTRIBUTE + BINOP.setParseAction(operator_function) +
         NUMBERPATTERN + Optional(RPAR)).setParseAction(AttributeSelection)

    selection_pattern << Optional(id_selection_expression) + \
        Optional(attr_selection_expression) + \
        ZeroOrMore(BOOLOP + selection_pattern)

    def __init__(self, atoms=None, selstr=None):
        super().__init__()
        self.atoms = atoms
        # self.ATOM_ATTRIBUTE = \
        #     oneOf(' '.join([attr for attr in dir(self.atoms) if not
        #                     attr.startswith('_')]))
        self.selstr = selstr
        self.fmtstr = "atoms={atoms!r}, selstr={selstr!r}"
        # self.selection_pattern = selection_pattern

        # self.selected_atoms = self.atoms.__class__(**self.atoms.kwargs)

        # if selstr is not None:
        #     self.parse(selstr)

    def parse(self, selstr=None, **kwargs):
        if selstr is None and self.selstr is not None:
            selstr = self.selstr

        try:
            selection = \
                self.selection_pattern.parseString(selstr, parseAll=True)[0]
        except ParseException as e:
            raise SelectionException(e.pstr, e.loc, e.msg, e.parseElement)
        else:
            return selection.apply(self.atoms)

    # def id_selection(self, s, l, t):
    #     # [self.selected_atoms.append(atom) for atom in
    #     #  self.atoms.filtered_ids(t[0].asList())]
    #     return self.atoms.filtered_ids(t[0].asList())

    # def attr_selection(self, s, l, t):
    #     print(t)

    def todict(self):
        return dict(atoms=self.atoms, selstr=self.selstr)
