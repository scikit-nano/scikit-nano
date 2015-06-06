# -*- coding: utf-8 -*-
"""
===============================================================================
Atoms class with velocity attributes (:mod:`sknano.core.atoms._velocity_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._velocity_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

# from sknano.core.math import Vector, transformation_matrix
from ._atoms import Atoms
from ._velocity_atom import VelocityAtom

__all__ = ['VelocityAtoms']


class VelocityAtoms(Atoms):
    """An eXtended `Atoms` class.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.VelocityAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `VelocityAtoms`}, optional
        if not `None`, then a list of `VelocityAtom` instance objects or an
        existing `VelocityAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return VelocityAtom

    def sort(self, key=attrgetter('v'), reverse=False):
        self.data.sort(key=key, reverse=reverse)
        super().sort(reverse=reverse)

    @property
    def velocities(self):
        """:class:`~numpy:numpy.ndarray` of `VelocityAtom` velocities."""
        return np.asarray([atom.v for atom in self])
