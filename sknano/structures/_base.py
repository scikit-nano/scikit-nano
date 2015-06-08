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

from sknano.core.crystallography import StructureBase as XtalBase
from sknano.core.refdata import aCC

__all__ = ['StructureBase']


class StructureBase:
    """Base class for creating abstract representation of nano structure.

    Parameters
    ----------
    basis : {str, int, list}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2

    """
    def __init__(self, *args, bond=None, verbose=False, debug=False,
                 **kwargs):
        if bond is None:
            bond = aCC
        self.bond = bond
        self.verbose = verbose
        self.debug = debug

        # self.fmtstr = "{lattice!r}, {basis!r}, {coords!r}, cartesian=True"
        super().__init__(*args, **kwargs)


    # def __repr__(self):
    #     return "{}({})".format(self.__class__.__name__,
    #                            self.fmtstr.format(**self.todict()))

    def __getattr__(self, name):
        try:
            return getattr(self.unit_cell, name)
        except AttributeError:
            return super().__getattr__(name)


    @property
    def element1(self):
        return self.basis[0].symbol

    @property
    def element2(self):
        return self.basis[1].symbol

