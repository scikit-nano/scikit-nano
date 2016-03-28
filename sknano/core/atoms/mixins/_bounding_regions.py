# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin classes for bounding regions (:mod:`sknano.core.atoms.mixins._geometric_regions`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins._geometric_regions

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
            return self.lattice.region.bounding_box
        except AttributeError:
            return Cuboid(pmin=[self.x.min(), self.y.min(), self.z.min()],
                          pmax=[self.x.max(), self.y.max(), self.z.max()])

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
    def bounding_region(self):
        """Bounding :class:`~sknano.core.geometric_regions.Geometric3DRegion`.

        Returns
        -------
        :class:`~sknano.core.geometric_regions.Geometric3DRegion`

        """
        try:
            return self.lattice.region
        except AttributeError:
            return self.bounding_box

    @property
    def bounds(self):
        """Alias for :attr:`~BoundingRegionsMixin.bounding_box`."""
        return self.bounding_box

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
