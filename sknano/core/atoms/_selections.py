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

from functools import reduce
from operator import and_, or_

# import numpy as np

# from abc import abstractmethod

from pyparsing import Forward, Word, Suppress, Optional, CaselessKeyword, \
    OneOrMore, ParseException, delimitedList, infixNotation, alphas, oneOf, \
    opAssoc

from sknano.core import BaseClass, binary_operator, integer, \
    number
from sknano.core.geometric_regions import Sphere
# from ._md_atoms import MDAtom as Atom, MDAtoms as Atoms

__all__ = ['SelectionException', 'SelectionParser']


def make_selection(cls):
    def selection(s, loc, tokens):
        return cls(tokens)
    return selection


class Selection(BaseClass):
    def __init__(self, selection):
        super().__init__()
        self.selection = selection
        self.fmtstr = "{selection!r}"

    def apply(self, atoms, filtered):
        klass, kwargs = atoms.__class__, atoms.kwargs
        return klass(atoms=filtered[:], **kwargs)

    def todict(self):
        return dict(selection=self.selection)


class SelectionException(ParseException):
    pass


class AndSelection(Selection):

    def apply(self, atoms):
        reduce_func = \
            lambda lsel, rsel: and_(lsel, rsel)
        return reduce(reduce_func,
                      [sel.apply(atoms) for sel in self.selection[0]])


class OrSelection(Selection):

    def apply(self, atoms):
        reduce_func = \
            lambda lsel, rsel: or_(lsel.apply(atoms), rsel.apply(atoms))
        return reduce(reduce_func, self.selection[0])


class NotSelection(Selection):

    def apply(self, atoms):
        return atoms - self.selection[0][0].apply(atoms)


class AllSelection(Selection):

    def apply(self, atoms):
        return super().apply(atoms, [atom for atom in atoms])


class NoneSelection(Selection):

    def apply(self, atoms):
        return super().apply(atoms, [])


class IDSelection(Selection):

    def apply(self, atoms):
        filtered = [atom for atom in atoms if atom.id
                    in self.selection.asList()]
        return super().apply(atoms, filtered)


class TypeSelection(Selection):

    def apply(self, atoms):
        filtered = [atom for atom in atoms if atom.type
                    in self.selection.asList()]
        return super().apply(atoms, filtered)


class AttributeSelection(Selection):

    def apply(self, atoms):
        try:
            attr, op, val = self.selection
            filtered = [atom for atom in atoms if op(getattr(atom, attr), val)]
        except ValueError:
            attr, val = self.selection
            filtered = [atom for atom in atoms if getattr(atom, attr) == val]
        return super().apply(atoms, filtered)


class WithinSelection(Selection):

    def apply(self, atoms):
        other = self.selection[-1].apply(atoms)
        filtered = atoms.query_ball_tree(other, self.selection[0])
        return super().apply(atoms, filtered)


class ExWithinSelection(Selection):

    def apply(self, atoms):
        other = self.selection[-1].apply(atoms)
        filtered = atoms.query_ball_tree(other, self.selection[0])
        filtered = super().apply(atoms, filtered) - other
        return super().apply(atoms, filtered)


class SelectionParser(BaseClass):
    """Selection parser class."""
    ALL = CaselessKeyword('all')
    NONE = CaselessKeyword('none')
    NAME = CaselessKeyword('name')
    TYPE = CaselessKeyword('type')
    INDEX = CaselessKeyword('index')
    ID = CaselessKeyword('id')
    WITHIN = CaselessKeyword('within')
    EXWITHIN = CaselessKeyword('exwithin')
    SERIAL = CaselessKeyword('serial')
    ATOMICNUMBER = CaselessKeyword('atomicnumber')
    ELEMENT = CaselessKeyword('element')
    RESIDUE = CaselessKeyword('residue')
    NUMBONDS = CaselessKeyword('numbonds')

    OF = CaselessKeyword('of')

    NOT = CaselessKeyword('not')
    AND = CaselessKeyword('and')
    OR = CaselessKeyword('or')

    LPAR, RPAR = map(Suppress, '()')

    ATOM_ATTRIBUTE = \
        oneOf(' '.join(('x', 'y', 'z', 'r', 'vx', 'vy', 'vz', 'v',
                        'fx', 'fy', 'fz', 'f', 'Z', 'element', 'symbol',
                        'mass', 'CN')))

    REGION = ()

    # selection_expression = Forward()

    expr_term = Forward()
    expr = \
        infixNotation(
            expr_term,
            [(Suppress(NOT), 1, opAssoc.RIGHT, make_selection(NotSelection)),
             (Suppress(AND), 2, opAssoc.LEFT, make_selection(AndSelection)),
             (Suppress(OR), 2, opAssoc.LEFT, make_selection(OrSelection))])
    grouped_expr = LPAR + expr_term + RPAR

    id_selection_expression = \
        (Suppress(ID) +
         delimitedList(OneOrMore(integer), delim=' ')
         ).setParseAction(make_selection(IDSelection))

    type_selection_expression = \
        (Suppress(TYPE) +
         delimitedList(OneOrMore(integer), delim=' ')
         ).setParseAction(make_selection(TypeSelection))

    attr_selection_expression = \
        (ATOM_ATTRIBUTE +
         Optional(binary_operator) +
         (number | Word(alphas))
         ).setParseAction(make_selection(AttributeSelection))

    within_selection_expression = \
        (Suppress(WITHIN) + number + Suppress(OF) +
         expr).setParseAction(make_selection(WithinSelection))

    exwithin_selection_expression = \
        (Suppress(EXWITHIN) + number + Suppress(OF) +
         expr).setParseAction(make_selection(ExWithinSelection))

    all_selection_expression = \
        (Suppress(ALL)).setParseAction(make_selection(AllSelection))

    expr_term << \
        (id_selection_expression | type_selection_expression |
         attr_selection_expression | within_selection_expression |
         exwithin_selection_expression | all_selection_expression |
         grouped_expr)

    selection_expression = expr.copy()

    def __init__(self, atoms=None, selstr=None, **kwargs):
        super().__init__(**kwargs)
        self.atoms = atoms
        self.selstr = selstr
        self.fmtstr = "atoms={atoms!r}, selstr={selstr!r}"

        # if selstr is not None:
        #     self.parse(selstr)

    def parse(self, selstr=None, **kwargs):
        if selstr is None and self.selstr is not None:
            selstr = self.selstr

        try:
            selection = \
                self.selection_expression.parseString(selstr, parseAll=True)[0]
            if self.verbose:
                print('selstr: {}'.format(selstr))
                print('selection: {}'.format(selection))
        except ParseException as e:
            raise SelectionException(e.pstr, e.loc, e.msg, e.parseElement)
        else:
            return selection.apply(self.atoms)

    def todict(self):
        return dict(atoms=self.atoms, selstr=self.selstr)
