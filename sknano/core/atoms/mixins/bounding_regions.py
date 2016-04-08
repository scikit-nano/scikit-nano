# -*- coding: utf-8 -*-
"""
========================================================================================
Mixin classes for bounding regions (:mod:`sknano.core.atoms.mixins.bounding_regions`)
========================================================================================

.. currentmodule:: sknano.core.atoms.mixins.bounding_regions

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numpy as np
# import pandas as pd

from sknano.core.geometric_regions import Cuboid, Sphere
# from sknano.core.math import Vector

__all__ = ['BoundingRegionsMixin']


class BoundingRegionsMixin:
    """Mixin `Atoms` class for computing bounding regions."""

    @property
    def bounding_box(self):
        """Axis-aligned bounding box of `Atoms`.

        Returns
        -------
        :class:`~sknano.core.geometric_regions.Cuboid`

        """
        try:
            return self.bounding_region.bounding_box
        except AttributeError:
            None

    @property
    def coordinates_bounding_box(self):
        """Bounding box of atom coordinates.

        Returns
        -------

        """
        try:
            pmin, pmax = self.r.minmax
            return Cuboid(pmin=pmin, pmax=pmax)
        except AttributeError:
            return None

    @property
    def bounding_region(self):
        """Bounding :class:`~sknano.core.geometric_regions.Geometric3DRegion`.

        Returns
        -------
        :class:`~sknano.core.geometric_regions.Geometric3DRegion`
            The :class:`~sknano.core.geometric_regions.Parallelepiped` attached
            to the :class:`~sknano.core.crystallography.Crystal3DLattice` of
            the :class:`~sknano.core.atoms.LatticeAtoms` class.

        """
        try:
            return self.lattice.region
        except AttributeError:
            return self.coordinates_bounding_box

    @property
    def bounds(self):
        """Alias for :attr:`~BoundingRegionsMixin.bounding_region`."""
        return self.bounding_region

    @property
    def bounding_sphere(self):
        """Bounding :class:`Sphere` of `Atoms`.

        Returns
        -------
        :class:`~sknano.core.geometric_regions.Sphere`

        """
        return Sphere(center=self.centroid,
                      r=np.max((self.r - self.centroid).norms))

    @property
    def lattice_region(self):
        """:attr:`~sknano.core.atoms.LatticeAtoms.lattice` region.

        Returns
        -------
        :class:`~sknano.core.geometric_regions.Parallelepiped`
            The :class:`~sknano.core.geometric_regions.Parallelepiped` attached
            to the :class:`~sknano.core.crystallography.Crystal3DLattice` of
            the :class:`~sknano.core.atoms.LatticeAtoms` class.

        """
        try:
            return self.lattice.region
        except AttributeError:
            return None

    @property
    def volume(self):
        """Volume of region containing atoms."""
        try:
            return self._volume
        except AttributeError:
            return self.bounding_region.volume

    @volume.setter
    def volume(self, value):
        self._volume = float(value)

    def within_region(self, region):
        """Returns new `Atoms` object containing atoms within `region`.

        Similar to :meth:`XYZAtom.clip_bounds`, but returns a new
        `Atoms` object.

        Parameters
        ----------
        region : :class:`~sknano.core.geometric_regions.Geometric3DRegion`

        """
        return \
            self.__class__([atom for atom in self if region.contains(atom.r)],
                           update_item_class=False, **self.kwargs)
