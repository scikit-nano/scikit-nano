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

import numbers

from sknano.core.refdata import CCbond

__all__ = ['StructureBase']


class StructureBase:
    """Base class for creating abstract representation of nano structure.

    Parameters
    ----------
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

    """
    def __init__(self, *args, element1='C', element2='C', bond=CCbond,
                 verbose=False, debug=False, **kwargs):

        self.element1 = element1
        self.element2 = element2

        if bond is None:
            bond = CCbond

        self.bond = bond
        self.verbose = verbose
        self.debug = debug

        super(StructureBase, self).__init__(*args, **kwargs)

    @property
    def bond(self):
        """Bond length in **\u212b**."""
        return self._bond

    @bond.setter
    def bond(self, value):
        """Set bond length in **\u212b**."""
        if not (isinstance(value, numbers.Real) or value > 0):
            raise TypeError('Expected a real, positive number.')
        self._bond = float(value)

    @bond.deleter
    def bond(self):
        del self._bond
