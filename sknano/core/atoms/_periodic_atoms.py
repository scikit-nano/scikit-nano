# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin classes for PBC (:mod:`sknano.core.atoms._periodic_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._periodic_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

__all__ = ['PBCAtomMixin', 'PBCAtomsMixin']


class PBCAtomMixin:
    """Mixin class for PBC."""
    @property
    def xperiodic(self):
        return self._xperiodic

    @xperiodic.setter
    def xperiodic(self, value):
        self._xperiodic = bool(value)

    @property
    def yperiodic(self):
        return self._yperiodic

    @yperiodic.setter
    def yperiodic(self, value):
        self._yperiodic = bool(value)

    @property
    def zperiodic(self):
        return self._zperiodic

    @zperiodic.setter
    def zperiodic(self, value):
        self._zperiodic = bool(value)


class PBCAtomsMixin:
    """Mixin :class:`~sknano.core.atoms.Atoms` class for PBC."""
    pass
