# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin class for PBCs (:mod:`sknano.core.atoms._periodic_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._periodic_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numpy as np

__all__ = ['PBCAtomsMixin']


class PBCAtomsMixin:
    """:class:`~sknano.core.atoms.Atoms` class for PBC."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._xperiodic = False
        self._yperiodic = False
        self._zperiodic = False

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

    @property
    def pbc(self):
        return np.asarray([self.xperiodic, self.yperiodic, self.zperiodic],
                          dtype=bool)

    @pbc.setter
    def pbc(self, value):
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self.xperiodic, self.yperiodic, self.zperiodic = value

    def set_pbc(self, xperiodic=True, yperiodic=True, zperiodic=True):
        self.xperiodic = xperiodic
        self.yperiodic = yperiodic
        self.zperiodic = zperiodic

    def unset_pbc(self):
        self.xperiodic = False
        self.yperiodic = False
        self.zperiodic = False

    def wrap_coords(self, pbc=None):
        try:
            [setattr(atom, 'r', self.lattice.wrap_cartesian_coordinate(
                     atom.r, pbc=pbc if pbc is not None else self.pbc))
             for atom in self]
        except AttributeError:
            pass
