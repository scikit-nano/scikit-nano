# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class for POAV analysis (:mod:`sknano.core.atoms._poav_atom`)
===============================================================================

An `Atom` class for structure analysis.

.. currentmodule:: sknano.core.atoms._poav_atom

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numbers

from ._bonds import Bonds
from ._kdtree_atom import KDTAtom

__all__ = ['POAVAtom']


class POAVAtom(KDTAtom):
    """An `Atom` class for POAV analysis.

    """
    _atomattrs = KDTAtom._atomattrs + \
        ['pyramidalization_angle', 'sigma_bond_angle', 'poav', 'poma']

    def __init__(self, **kwargs):
        super(POAVAtom, self).__init__(**kwargs)

        self._pyramidalization_angle = None
        self._sigma_bond_angle = None
        self._misalignment_angles = None
        self.poav = None
        self.poma = []

    def __str__(self):
        return repr(self)

    def __repr__(self):
        """Return canonical string representation of `POAVAtom`."""
        return super(POAVAtom, self).__repr__()

    @property
    def sigma_bond_angle(self):
        return self._sigma_bond_angle

    @sigma_bond_angle.setter
    def sigma_bond_angle(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._sigma_bond_angle = value

    @property
    def pyramidalization_angle(self):
        """Return the pyramidalization angle :math:`\\theta_P`"""
        return self._pyramidalization_angle

    @pyramidalization_angle.setter
    def pyramidalization_angle(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._pyramidalization_angle = value

    #@property
    #def poav_misalignment_angle(self):
    #    return self._poav_misalignment_angle

    #@poav_misalignment_angle.setter
    #def poav_misalignment_angle(self, value):
    #    if not isinstance(value, ):
    #        raise TypeError('Expected a number')
    #    self._poav_misalignment_angle = value
