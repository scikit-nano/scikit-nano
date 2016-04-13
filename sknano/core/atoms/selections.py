# -*- coding: utf-8 -*-
"""
===============================================================================
Atom selection classes (:mod:`sknano.core.atoms.selections`)
===============================================================================

.. currentmodule:: sknano.core.atoms.selections

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from functools import reduce
from importlib import import_module
from operator import and_, or_

import numpy as np

# from abc import abstractmethod

from pyparsing import CaselessKeyword, Forward, Keyword, OneOrMore, Optional, \
    ParseException, Suppress, Word, alphas, delimitedList, infixNotation, \
    oneOf, opAssoc

from sknano.core import BaseClass, binary_operator, integer, kwargs_expr, \
    number
# from sknano.core.geometric_regions import Sphere
# from .md_atoms import MDAtom as Atom, MDAtoms as Atoms

__all__ = ['SelectionException', 'SelectionParser',
           'generate_vmd_selection_string']


def make_selection(cls):
    def selection(s, loc, tokens):
        return cls(tokens)
    return selection


def as_region(s, l, t):
    return \
        getattr(import_module('sknano.core.geometric_regions'), t[0])(**t[1])


class Selection(BaseClass):
    """Base selection class."""
    def __init__(self, selection):
        super().__init__()
        self.selection = selection
        self.fmtstr = "{selection!r}"

    def apply(self, atoms, filtered=None, mask=None):
        if mask is not None and isinstance(mask, np.ndarray) and \
                (mask.dtype in (bool, int)):
            return atoms.filtered(mask)
        elif filtered is not None:
            return atoms.__class__(atoms=filtered[:], **atoms.kwargs)

    def todict(self):
        return dict(selection=self.selection)


class SelectionException(ParseException):
    """Custom :class:`Exception` class for :class:`Selection`\ s."""
    pass


class AndSelection(Selection):
    def apply(self, atoms, as_mask=False):
        mask = None
        try:
            mask = reduce(lambda lsel, rsel: np.bitwise_and(lsel, rsel),
                          [sel.apply(atoms, as_mask=True) for sel in
                           self.selection[0]])
            if as_mask:
                return mask
        except (TypeError, ValueError):
            filtered = reduce(lambda lsel, rsel: and_(lsel, rsel),
                              [sel.apply(atoms) for sel in self.selection[0]])

        if mask is not None:
            return super().apply(atoms, mask=mask)
        else:
            return super().apply(atoms, filtered=filtered)


class OrSelection(Selection):

    def apply(self, atoms, as_mask=False):
        mask = None
        try:
            mask = reduce(lambda lsel, rsel: np.bitwise_or(lsel, rsel),
                          [sel.apply(atoms, as_mask=True) for sel in
                           self.selection[0]])
            if as_mask:
                return mask
        except (TypeError, ValueError):
            filtered = reduce(lambda lsel, rsel: or_(lsel, rsel),
                              [sel.apply(atoms) for sel in self.selection[0]])

        if mask is not None:
            return super().apply(atoms, mask=mask)
        else:
            return super().apply(atoms, filtered=filtered)


class NotSelection(Selection):

    def apply(self, atoms):
        return atoms[:] - self.selection[0][0].apply(atoms[:])


class AllSelection(Selection):
    def apply(self, atoms):
        return atoms[:]


class NoneSelection(Selection):

    def apply(self, atoms):
        # return super().apply(atoms, [])
        return None


class IDSelection(Selection):

    def apply(self, atoms, as_mask=False):
        mask = np.in1d(atoms.ids, self.selection.asList()).nonzero()[0]
        if as_mask:
            return mask
        return super().apply(atoms, mask=mask)


class MolIDSelection(Selection):

    def apply(self, atoms, as_mask=False):
        mask = np.in1d(atoms.mol_ids, self.selection.asList()).nonzero()[0]
        if as_mask:
            return mask
        return super().apply(atoms, mask=mask)


class TypeSelection(Selection):

    def apply(self, atoms, as_mask=False):
        mask = np.in1d(atoms.types, self.selection.asList()).nonzero()[0]
        if as_mask:
            return mask
        return super().apply(atoms, mask=mask)


class AttributeSelection(Selection):

    def apply(self, atoms, as_mask=False):
        try:
            attr, op, val = self.selection
            if hasattr(atoms, attr):
                mask = op(getattr(atoms, attr), val)
            else:
                mask = \
                    np.asarray([op(getattr(atom, attr), val)
                                for atom in atoms])
        except ValueError:
            attr, val = self.selection
            if hasattr(atoms, attr + 's'):
                mask = getattr(atoms, attr + 's') == val
            elif hasattr(atoms, attr):
                mask = getattr(atoms, attr) == val
            else:
                mask = \
                    np.asarray([True if getattr(atom, attr) == val else False
                                for atom in atoms])
        if as_mask:
            return mask
        else:
            return super().apply(atoms, mask=mask)


class WithinSelection(Selection):

    def apply(self, atoms):
        try:
            other = self.selection[-1].apply(atoms)
            filtered = atoms.query_ball_tree(other, self.selection[0])
            return super().apply(atoms, filtered=filtered)
        except AttributeError:
            region = self.selection[-1]
            mask = \
                np.asarray([region.contains(atom.r) for atom in atoms])
            return super().apply(atoms, mask=mask)


class ExWithinSelection(Selection):

    def apply(self, atoms):
        other = self.selection[-1].apply(atoms)
        filtered = atoms.query_ball_tree(other, self.selection[0])
        filtered = super().apply(atoms, filtered=filtered) - other
        return super().apply(atoms, filtered=filtered)


class SelectionParser(BaseClass):
    """Selection parser class."""
    ALL = CaselessKeyword('all')
    NONE = CaselessKeyword('none')
    NAME = CaselessKeyword('name')
    TYPE = CaselessKeyword('type')
    INDEX = CaselessKeyword('index')
    ID = CaselessKeyword('id')
    MOLID = CaselessKeyword('molid') | CaselessKeyword('mol')
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

    PARALLELEPIPED, CUBOID, CUBE, ELLIPSOID, SPHERE, CYLINDER, CONE = \
        map(Keyword, ['Parallelepiped', 'Cuboid', 'Cube', 'Ellipsoid',
                      'Sphere', 'Cylinder', 'Cone'])

    REGIONS = (PARALLELEPIPED | CUBOID | CUBE |
               ELLIPSOID | SPHERE | CYLINDER | CONE)

    region_expression = \
        (REGIONS + LPAR + kwargs_expr + RPAR).setParseAction(as_region)

    # selection_expression = Forward()

    expr_term = Forward()

    expr = \
        infixNotation(expr_term,
                      [(Suppress(NOT), 1, opAssoc.RIGHT, NotSelection),
                       (Suppress(AND), 2, opAssoc.LEFT, AndSelection),
                       (Suppress(OR), 2, opAssoc.LEFT, OrSelection)])

    grouped_expr = LPAR + expr_term + RPAR

    all_selection_expression = (Suppress(ALL)).setParseAction(AllSelection)
    none_selection_expression = (Suppress(NONE)).setParseAction(NoneSelection)

    attr_selection_expression = \
        (ATOM_ATTRIBUTE + Optional(binary_operator) + (number | Word(alphas))
         ).setParseAction(AttributeSelection)

    id_selection_expression = \
        (Suppress(ID) + delimitedList(OneOrMore(integer), delim=' ')
         ).setParseAction(IDSelection)

    molid_selection_expression = \
        (Suppress(MOLID) + delimitedList(OneOrMore(integer), delim=' ')
         ).setParseAction(MolIDSelection)

    type_selection_expression = \
        (Suppress(TYPE) + delimitedList(OneOrMore(integer), delim=' ')
         ).setParseAction(TypeSelection)

    within_selection_expression = \
        (Suppress(WITHIN) +
         ((number + Suppress(OF) + expr) | region_expression)
         ).setParseAction(WithinSelection)

    exwithin_selection_expression = \
        (Suppress(EXWITHIN) +
         number + Suppress(OF) + expr
         ).setParseAction(ExWithinSelection)

    expr_term << \
        (all_selection_expression | none_selection_expression |
         attr_selection_expression | id_selection_expression |
         molid_selection_expression | type_selection_expression |
         within_selection_expression | exwithin_selection_expression |
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
        """Parse `selstr`."""
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


def generate_vmd_selection_string(keyword, iterable):
    """Generate a VMD `keyword` selection string from iterable."""
    return ' '.join((keyword, ' '.join(map(str, iterable))))
