# -*- coding: utf-8 -*-
"""
===============================================================================
Topology Mixin (:mod:`sknano.core.atoms.mixins._topology`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins._topology

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# from operator import attrgetter

import numpy as np

# from sknano.core import BaseClass, UserList, cyclic_pairs, dedupe
# from sknano.core.math import Vector, vector as vec

# import sknano.core.atoms

# from sknano.core import flatten
from sknano.core.math import angle, dot

__all__ = ['AngleMixin', 'BondMixin', 'DihedralMixin', 'ImproperMixin',
           'AtomTopologyMixin', 'AtomsTopologyMixin',
           'get_angle', 'get_dihedral']


def get_angle(atoms):
    # if all([p is None for p in (atoms, bonds)]):
    #     raise ValueError('Expected a sequence of `atoms` or `bonds`')
    # if bonds is not None:
    #     atoms = Bonds(bonds).atoms
    atom1, atom2, atom3 = atoms
    r21 = atom1.r - atom2.r
    r23 = atom3.r - atom2.r
    return angle(r21, r23)


def get_dihedral(atoms):
    atom1, atom2, atom3, atom4 = atoms
    r21 = atom1.r - atom2.r
    r23 = atom3.r - atom2.r
    r34 = atom4.r - atom3.r
    m = r21.cross(r23)
    n = r23.cross(r34)

    cos1234 = dot(m, n) / (m.norm * n.norm)
    sin1234 = dot(n, r21) * r34.norm / (m.norm * n.norm)
    return -np.arctan2(sin1234, cos1234)

    # phi = angle(m, n)
    # return (phi if scalar_triple_product(r21, r23, r34) <= 0.0 else -phi)


class AngleMixin:
    def get_angle(self, *atoms):
        pass


class BondMixin:
    @property
    def bonds(self):
        """Return atom `Bonds` instance."""
        from .._bonds import Bond, Bonds
        try:
            return Bonds([Bond(self, nn) for nn in self.NN])
        except (AttributeError, TypeError):
            return Bonds()


class DihedralMixin:
    """Mixin class for dihedral angle calculations."""

    def get_dihedral(self, *atoms):
        pass


class ImproperMixin:

    def get_improper(self, *atoms):
        pass


class AtomTopologyMixin(ImproperMixin, DihedralMixin, BondMixin, AngleMixin):
    pass


class AtomsTopologyMixin:

    @property
    def bonds(self):
        from .._bonds import Bonds
        return Bonds([bond for atom in self for bond in atom.bonds])

    @property
    def max_bond(self):
        return np.max(self.bonds.lengths)
