# -*- coding: utf-8 -*-
"""
===============================================================================
Structure selection classes (:mod:`sknano.core.structures.selections`)
===============================================================================

.. currentmodule:: sknano.core.structures.selections

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# from functools import reduce
# from importlib import import_module
# from operator import and_, or_

# import numpy as np

# from abc import abstractmethod

# from pyparsing import CaselessKeyword, Forward, Keyword, OneOrMore, \
#     Optional, ParseException, Suppress, Word, alphas, delimitedList, \
#     infixNotation, oneOf, opAssoc

# from sknano.core import BaseClass, binary_operator, integer, kwargs_expr, \
#     number
# from sknano.core.atoms import AtomsSelectionMixin
# from sknano.core.geometric_regions import Sphere

__all__ = ['StructureSelectionMixin']


# def make_selection(cls):
#     def selection(s, loc, tokens):
#         return cls(tokens)
#     return selection


# def as_region(s, l, t):
#     return \
#         getattr(import_module('sknano.core.geometric_regions'), t[0])(**t[1])


# class StructureSelection(BaseClass):
#     """Base selection class."""
#     def __init__(self, selection):
#         super().__init__()
#         self.selection = selection
#         self.fmtstr = "{selection!r}"

#     def apply(self, structure, filtered=None, mask=None):
#         atoms = structure.atoms
#         if mask is not None and isinstance(mask, np.ndarray) and \
#                 (mask.dtype in (bool, int)):
#             atoms = atoms.filtered(mask)
#         elif filtered is not None:
#             atoms = atoms.__class__(atoms=filtered[:], **atoms.kwargs)
#         structure.atoms = atoms
#         return structure

#     def todict(self):
#         return dict(selection=self.selection)


# class StructureSelectionException(ParseException):
#     """Custom :class:`Exception` class for :class:`StructureSelection`\ s."""
#     pass


# class AndSelection(StructureSelection):
#     def apply(self, structure, as_mask=False):
#         mask = None
#         try:
#             mask = reduce(lambda lsel, rsel: np.bitwise_and(lsel, rsel),
#                           [sel.apply(structure, as_mask=True) for sel in
#                            self.selection[0]])
#             if as_mask:
#                 return mask
#         except (TypeError, ValueError):
#             filtered = reduce(lambda lsel, rsel: and_(lsel, rsel),
#                               [sel.apply(structure) for sel in
#                                self.selection[0]])

#         if mask is not None:
#             return super().apply(structure, mask=mask)
#         else:
#             return super().apply(structure, filtered=filtered)


# class OrSelection(StructureSelection):

#     def apply(self, structure, as_mask=False):
#         mask = None
#         try:
#             mask = reduce(lambda lsel, rsel: np.bitwise_or(lsel, rsel),
#                           [sel.apply(structure, as_mask=True) for sel in
#                            self.selection[0]])
#             if as_mask:
#                 return mask
#         except (TypeError, ValueError):
#             filtered = reduce(lambda lsel, rsel: or_(lsel, rsel),
#                               [sel.apply(structure) for sel in
#                                self.selection[0]])

#         if mask is not None:
#             return super().apply(structure, mask=mask)
#         else:
#             return super().apply(structure, filtered=filtered)


# class NotSelection(StructureSelection):

#     def apply(self, structure):
#         structure.atoms = \
#             structure.atoms[:] - self.selection[0][0].apply(structure).atoms
#         return structure


# class AllSelection(StructureSelection):
#     def apply(self, structure):
#         return structure


# class NoneSelection(StructureSelection):

#     def apply(self, structure):
#         # return super().apply(structure, [])
#         return None


# class IDSelection(StructureSelection):

#     def apply(self, structure, as_mask=False):
#         atoms = structure.atoms
#         mask = np.in1d(atoms.ids, self.selection.asList()).nonzero()[0]
#         if as_mask:
#             return mask
#         return super().apply(structure, mask=mask)


# class AttributeSelection(StructureSelection):

#     def apply(self, structure, as_mask=False):
#         atoms = structure.atoms
#         try:
#             attr, op, val = self.selection
#             if hasattr(atoms, attr):
#                 mask = op(getattr(atoms, attr), val)
#             else:
#                 mask = \
#                     np.asarray([op(getattr(atom, attr), val)
#                                 for atom in atoms])
#         except ValueError:
#             attr, val = self.selection
#             if hasattr(atoms, attr + 's'):
#                 mask = getattr(atoms, attr + 's') == val
#             elif hasattr(atoms, attr):
#                 mask = getattr(atoms, attr) == val
#             else:
#                 mask = \
#                     np.asarray([True if getattr(atom, attr) == val else False
#                                 for atom in atoms])
#         if as_mask:
#             return mask
#         else:
#             return super().apply(structure, mask=mask)


# class WithinSelection(StructureSelection):
#     """:class:`StructureSelection` class for selections within regions or \
#         distance."""
#     def apply(self, structure):
#         atoms = structure.atoms
#         try:
#             other = self.selection[-1].apply(structure).atoms
#             filtered = atoms.query_ball_tree(other, self.selection[0])
#             return super().apply(structure, filtered=filtered)
#         except AttributeError:
#             region = self.selection[-1]
#             mask = \
#                 np.asarray([region.contains(atom.r) for atom in atoms])
#             return super().apply(structure, mask=mask)


# class AtomsSelection(StructureSelection):
#     """Exclusive within :class:`StructureSelection` class."""
#     def apply(self, structure):
#         atoms = structure.atoms
#         other = self.selection[-1].apply(structure).atoms
#         filtered = atoms.query_ball_tree(other, self.selection[0])
#         filtered = super().apply(structure, filtered=filtered).atoms - other
#         return super().apply(structure, filtered=filtered)


# class StructureSelectionParser(BaseClass):
#     """StructureSelection parser class."""
#     ALL = CaselessKeyword('all')
#     NONE = CaselessKeyword('none')
#     NAME = CaselessKeyword('name')
#     ATOMS = CaselessKeyword('atoms')
#     WITHIN = CaselessKeyword('within')
#     EXWITHIN = CaselessKeyword('exwithin')

#     OF = CaselessKeyword('of')

#     NOT = CaselessKeyword('not')
#     AND = CaselessKeyword('and')
#     OR = CaselessKeyword('or')

#     LPAR, RPAR = map(Suppress, '()')

#     ATOM_ATTRIBUTE = \
#         oneOf(' '.join(('x', 'y', 'z', 'r', 'vx', 'vy', 'vz', 'v',
#                         'fx', 'fy', 'fz', 'f', 'Z', 'element', 'symbol',
#                         'mass', 'CN')))

#     PARALLELEPIPED, CUBOID, CUBE, ELLIPSOID, SPHERE, CYLINDER, CONE = \
#         map(Keyword, ['Parallelepiped', 'Cuboid', 'Cube', 'Ellipsoid',
#                       'Sphere', 'Cylinder', 'Cone'])

#     REGIONS = (PARALLELEPIPED | CUBOID | CUBE |
#                ELLIPSOID | SPHERE | CYLINDER | CONE)

#     region_expression = \
#         (REGIONS + LPAR + kwargs_expr + RPAR).setParseAction(as_region)

#     # selection_expression = Forward()

#     expr_term = Forward()

#     expr = \
#         infixNotation(expr_term,
#                       [(Suppress(NOT), 1, opAssoc.RIGHT, NotSelection),
#                        (Suppress(AND), 2, opAssoc.LEFT, AndSelection),
#                        (Suppress(OR), 2, opAssoc.LEFT, OrSelection)])

#     grouped_expr = LPAR + expr_term + RPAR

#     atoms_selection_expression = \
#         (Suppress(ATOMS))

#     all_selection_expression = (Suppress(ALL)).setParseAction(AllSelection)
#     none_selection_expression = (Suppress(NONE)).setParseAction(NoneSelection)

#     attr_selection_expression = \
#         (ATOM_ATTRIBUTE + Optional(binary_operator) + (number | Word(alphas))
#          ).setParseAction(AttributeSelection)

#     id_selection_expression = \
#         (Suppress(ID) + delimitedList(OneOrMore(integer), delim=' ')
#          ).setParseAction(IDSelection)

#     within_selection_expression = \
#         (Suppress(WITHIN) +
#          ((number + Suppress(OF) + expr) | region_expression)
#          ).setParseAction(WithinSelection)

#     exwithin_selection_expression = \
#         (Suppress(EXWITHIN) +
#          number + Suppress(OF) + expr
#          ).setParseAction(ExWithinSelection)

#     expr_term << \
#         (all_selection_expression | none_selection_expression |
#          atoms_selection_expression | grouped_expr)

#     selection_expression = expr.copy()

#     def __init__(self, structure=None, selstr=None, **kwargs):
#         super().__init__(**kwargs)
#         self.structure = structure
#         self.selstr = selstr
#         self.fmtstr = "structure={structure!r}, selstr={selstr!r}"

#         # if selstr is not None:
#         #     self.parse(selstr)

#     def parse(self, selstr=None, **kwargs):
#         """Parse `selstr`."""
#         if selstr is None and self.selstr is not None:
#             selstr = self.selstr

#         try:
#             selection = \
#                 self.selection_expression.parseString(selstr, parseAll=True)[0]
#             if self.verbose:
#                 print('selstr: {}'.format(selstr))
#                 print('selection: {}'.format(selection))
#         except ParseException as e:
#             raise StructureSelectionException(e.pstr, e.loc, e.msg,
#                                               e.parseElement)
#         else:
#             return selection.apply(self.structure)

#     def todict(self):
#         """Return :class:`~python:dict` of constructor parameters."""
#         return dict(structure=self.structure, selstr=self.selstr)


class StructureSelectionMixin:
    """Mixin class for applying selections."""
    def select(self, selstr=None, selstrlist=None, verbose=False):
        """Return `Atom` or `Atoms` from selection command.

        Parameters
        ----------
        selstr : :class:`~python:str`, optional
            optional if `selstrlist` is not `None`
        selstrlist : {`None`, :class:`~python:list`}, optional
            :class:`~python:list` of selection strings.

        Returns
        -------
        :class:`~python:list` of `Atom` or `Atoms` objects
            if `selstrlist` is not `None`
        :class:`Atom` or :class:`Atoms` if `selstr` is not `None`

        """
        return self.atoms.select(selstr, selstrlist, verbose)
