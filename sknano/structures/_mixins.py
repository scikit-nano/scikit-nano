# -*- coding: utf-8 -*-
"""
==============================================================================
Mixin structure classes (:mod:`sknano.structures._mixins`)
==============================================================================

.. currentmodule:: sknano.structures._mixins

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from sknano.core.atoms import vdw_radius_from_basis

__all__ = ['BasisMixin']


class BasisMixin:
    """Mixin class for structure atoms basis."""
    @property
    def vdw_radius(self):
        if self._vdw_radius is not None:
            return self._vdw_radius
        else:
            return vdw_radius_from_basis(self.basis[0], self.basis[1])

    @vdw_radius.setter
    def vdw_radius(self, value):
        self._vdw_radius = value

    @property
    def vdw_distance(self):
        return 2 * self.vdw_radius
