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

from abc import ABCMeta, abstractmethod

from sknano.core.refdata import aCC

__all__ = ['StructureBase']


class StructureBase(metaclass=ABCMeta):
    """Base class for creating abstract representation of nano structure.

    Parameters
    ----------
    basis : {str, int, list}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2

    """
    def __init__(self, *args, basis=None, bond=None, verbose=False,
                 debug=False, **kwargs):

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

        self.verbose = verbose
        self.debug = debug

        self.fmtstr = ""

        super().__init__(*args, **kwargs)

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    def __getattr__(self, name):
        try:
            return getattr(self.unit_cell, name)
        except AttributeError:
            return super().__getattr__(name)

    @property
    def fmtstr(self):
        return self._fmtstr

    @fmtstr.setter
    def fmtstr(self, value):
        self._fmtstr = value

    @property
    def element1(self):
        return self.basis[0]

    @element1.setter
    def element1(self, value):
        self.basis[0] = value
        self.unit_cell.basis[0] = value

    @property
    def element2(self):
        return self.basis[1]

    @element2.setter
    def element2(self, value):
        self.basis[1] = value
        self.unit_cell.basis[1] = value

    @abstractmethod
    def todict(self):
        return NotImplementedError
