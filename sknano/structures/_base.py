# -*- coding: utf-8 -*-
"""
==============================================================================
Base structure classes (:mod:`sknano.structures._base`)
==============================================================================

.. currentmodule:: sknano.structures._base

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

# import numbers

from sknano.core import BaseClass
from sknano.core.refdata import aCC

__all__ = ['StructureBase']


class StructureBase(BaseClass):
    """Base class for creating abstract representations of nanostructure.

    Parameters
    ----------
    basis : {:class:`~python:str`, :class:`~python:int`, list}, optional
        Element chemical symbols or atomic numbers of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2

    """
    def __init__(self, *args, basis=None, bond=None, **kwargs):

        if basis is None:
            basis = ['C', 'C']
        self.basis = basis[:]

        if 'element1' in kwargs:
            self.basis[0] = kwargs['element1']
            del kwargs['element1']

        if 'element2' in kwargs:
            self.basis[1] = kwargs['element2']
            del kwargs['element2']

        if bond is None:
            bond = aCC
        self.bond = bond

        super().__init__(*args, **kwargs)

    def __getattr__(self, name):
        try:
            return getattr(self.unit_cell, name)
        except AttributeError:
            return super().__getattr__(name)

    @property
    def element1(self):
        "Basis element 1"
        return self.basis[0]

    @element1.setter
    def element1(self, value):
        self.basis[0] = value
        [self.unit_cell.basis.__setitem__(i, value)
         for i in range(0, len(self.unit_cell.basis), 2)]

    @property
    def element2(self):
        return self.basis[1]

    @element2.setter
    def element2(self, value):
        self.basis[1] = value
        [self.unit_cell.basis.__setitem__(i, value)
         for i in range(1, len(self.unit_cell.basis), 2)]
