# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes for PBCs (:mod:`sknano.core.atoms._periodic_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._periodic_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from functools import total_ordering
from operator import attrgetter

import numpy as np

from ._atoms import Atom, Atoms

__all__ = ['PBCAtom', 'PBCAtoms']


@total_ordering
class PBCAtom(Atom):
    """Atom class with PBC attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    xperiodic, yperiodic, zperiodic : :class:`~python:bool`, optional
        PBC along each :math:`x, y, z` axis
    """
    def __init__(self, *args, xperiodic=False, yperiodic=False,
                 zperiodic=False, **kwargs):
        super().__init__(*args, **kwargs)

        self.xperiodic = xperiodic
        self.yperiodic = yperiodic
        self.zperiodic = zperiodic

        self.fmtstr = super().fmtstr + \
            ", xperiodic={xperiodic!r}, yperiodic={yperiodic!r}" + \
            ", zperiodic={zperiodic!r}"

    def __eq__(self, other):
        return np.allclose(self.pbc, other.pbc) and super().__eq__(other)

    def __lt__(self, other):
        return (self.pbc < other.pbc and super().__le__(other)) or \
            (self.pbc <= other.pbc and super().__lt__(other))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['xperiodic', 'yperiodic', 'zperiodic'])
        return attrs

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

    def todict(self):
        pdict = super().todict()
        pdict.update(dict(xperiodic=self.xperiodic,
                          yperiodic=self.yperiodic,
                          zperiodic=self.zperiodic))
        return pdict


class PBCAtoms(Atoms):
    """:class:`~sknano.core.atoms.Atoms` class for PBC."""
    @property
    def __atom_class__(self):
        return PBCAtom

    def sort(self, key=attrgetter('pbc'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def pbc(self):
        return np.asarray([atom.pbc for atom in self])

    def set_pbc(self, xperiodic=True, yperiodic=True, zperiodic=True):
        [setattr(atom, 'xperiodic', xperiodic) for atom in self]
        [setattr(atom, 'yperiodic', yperiodic) for atom in self]
        [setattr(atom, 'zperiodic', zperiodic) for atom in self]

    def unset_pbc(self):
        [setattr(atom, 'xperiodic', False) for atom in self]
        [setattr(atom, 'yperiodic', False) for atom in self]
        [setattr(atom, 'zperiodic', False) for atom in self]

    def wrap_coords(self, pbc=None):
        try:
            [setattr(atom, 'r', self.lattice.wrap_cartesian_coordinate(
                     atom.r, pbc=pbc if pbc is not None else atom.pbc))
             for atom in self]
        except AttributeError:
            pass
