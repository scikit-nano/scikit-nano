# -*- coding: utf-8 -*-
"""
===============================================================================
Atom sub-class for POAV analysis (:mod:`sknano.core.atoms._poav_atom`)
===============================================================================

An `Atom` sub-class for POAV analysis.

.. currentmodule:: sknano.core.atoms._poav_atom

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from sknano.utils.analysis import POAV1, POAV2, POAVR
from ._kdtree_atom import KDTAtom

__all__ = ['POAVAtomMixin', 'POAVAtom']


class POAVAtomMixin(object):
    """Mixin class for `POAVAtom`."""
    @property
    def POAV1(self):
        """:class:`~sknano.utils.analysis.POAV1` instance."""
        return self._POAV1

    @POAV1.setter
    def POAV1(self, value):
        """Set :class:`~sknano.utils.analysis.POAV1` instance."""
        if not isinstance(value, POAV1):
            raise TypeError('Expected a `POAV1` instance.')
        self._POAV1 = value

    @property
    def POAV2(self):
        """:class:`~sknano.utils.analysis.POAV2` instance."""
        return self._POAV2

    @POAV2.setter
    def POAV2(self, value):
        """Set :class:`~sknano.utils.analysis.POAV2` instance."""
        if not isinstance(value, POAV2):
            raise TypeError('Expected a `POAV2` instance.')
        self._POAV2 = value

    @property
    def POAVR(self):
        """:class:`~sknano.utils.analysis.POAVR` instance."""
        return self._POAVR

    @POAVR.setter
    def POAVR(self, value):
        """Set :class:`~sknano.utils.analysis.POAVR` instance."""
        if not isinstance(value, POAVR):
            raise TypeError('Expected a `POAVR` instance.')
        self._POAVR = value


class POAVAtom(POAVAtomMixin, KDTAtom):
    """An `Atom` sub-class for POAV analysis."""
    _atomattrs = KDTAtom._atomattrs + ['POAV1', 'POAV2', 'POAVR']

    def __init__(self, **kwargs):
        super(POAVAtom, self).__init__(**kwargs)
        self._POAV1 = None
        self._POAV2 = None
        self._POAVR = None
