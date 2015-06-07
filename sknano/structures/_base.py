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


class StructureBase(XtalBase):
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

        super().__init__(*args, **kwargs)
