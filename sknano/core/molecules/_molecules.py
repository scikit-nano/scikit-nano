# -*- coding: utf-8 -*-
"""
==============================================================================
Base class for structure molecules (:mod:`sknano.core.molecules._molecules`)
==============================================================================

.. currentmodule:: sknano.core.molecules._molecules

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import OrderedDict
from operator import attrgetter

import numpy as np

from sknano.core import UserList, xyz
from sknano.core.math import Vector, transformation_matrix
from sknano.core.geometric_regions import Cuboid  # , Rectangle

__all__ = ['Molecules']


class Molecules(UserList):
    """Base class for collection of `Molecule` objects.

    Parameters
    ----------
    molecules : {None, sequence, `Molecules`}, optional
        if not `None`, then a list of `Molecule` instance objects or an
        existing `Molecules` instance object.

    """
    _moleculeattrs = []

    def __init__(self, molecules=None):
        super().__init__(initlist=molecules)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        """Return canonical string representation of `Molecules`."""
        return "Molecules(molecules={!r})".format(self.data)

    def sort(self, key=attrgetter('id'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def Nmolecules(self):
        """Number of molecules in `Molecules`."""
        return len(self)

    @property
    def CM(self):
        """Center-of-Mass coordinates of `Molecules`.

        Returns
        -------
        ndarray
            3-element ndarray specifying center-of-mass coordinates of
            `Molecules`.

        """
        masses = np.asarray([self.masses])
        coords = self.coords
        MxR = masses.T * coords
        return Vector(np.sum(MxR, axis=0) / np.sum(masses))

    @property
    def M(self):
        """Total mass of `Molecules`."""
        #return math.fsum(self.masses)
        return self.masses.sum()

    @property
    def coords(self):
        """Return list of `Molecule` coordinates."""
        return np.asarray([molecule.r for molecule in self])

    @property
    def masses(self):
        """Return list of `Molecule` masses."""
        return np.asarray([molecule.m for molecule in self])

    @property
    def symbols(self):
        """Return list of `Molecule` symbols."""
        return np.asarray([molecule.symbol for molecule in self])

    @property
    def x(self):
        """Return :math:`x` coordinates of `Molecule` objects as array."""
        return self.coords[:,0]

    @property
    def y(self):
        """Return :math:`y` coordinates of `Molecule` objects as array."""
        return self.coords[:,1]

    @property
    def z(self):
        """Return :math:`z` coordinates of `Molecule` objects as array."""
        return self.coords[:,2]

    @property
    def bounds(self):
        """Return bounds of `Molecules`."""
        return Cuboid(pmin=[self.x.min(), self.y.min(), self.z.min()],
                      pmax=[self.x.max(), self.y.max(), self.z.max()])

    def center_CM(self, axes=None):
        """Center molecules on CM coordinates."""
        dr = -self.CM
        self.translate(dr)

    def clip_bounds(self, region, center_before_clipping=False):
        """Remove molecules outside the given limits along given dimension.

        Parameters
        ----------
        region : :class:`~sknano.core.geometric_regions.`GeometricRegion`

        """
        CM0 = None
        if center_before_clipping:
            CM0 = self.CM
            self.translate(-CM0)

        self.data = \
            np.asarray(self)[np.logical_and(
                np.logical_and(
                    self.x <= region.limits['x']['max'],
                    np.logical_and(
                        self.y <= region.limits['y']['max'],
                        self.z <= region.limits['z']['max'])),
                np.logical_and(
                    self.x >= region.limits['x']['min'],
                    np.logical_and(
                        self.y >= region.limits['y']['min'],
                        self.z >= region.limits['z']['min'])))].tolist()

        if CM0 is not None:
            self.translate(CM0)

    def filter(self, condition, invert=False):
        """Filter `Molecules` by `condition`.

        Parameters
        ----------
        condition : array_like, bool
        invert : bool, optional

        Returns
        -------
        filtered_molecules : `Molecules`

        """
        return self.__class__(molecules=np.asarray(self)[condition].tolist())

    def get_molecules(self, asarray=False):
        """Return list of `Molecules`.

        Parameters
        ----------
        asarray : bool, optional

        Returns
        -------
        sequence or ndarray

        """
        if asarray:
            return np.asarray(self)
        else:
            return self

    def get_coords(self, asdict=False):
        """Return molecule coords.

        Parameters
        ----------
        asdict : bool, optional

        Returns
        -------
        coords : :py:class:`python:~collections.OrderedDict` or ndarray

        """
        coords = self.coords
        if asdict:
            return OrderedDict(list(zip(xyz, coords.T)))
        else:
            return coords

    def rezero_coords(self, epsilon=1.0e-10):
        """Alias for :meth:`Molecules.rezero`."""
        self.rezero(epsilon=epsilon)

    def rezero_xyz(self, epsilon=1.0e-10):
        self.rezero(epsilon=epsilon)

    def rezero(self, epsilon=1.0e-10):
        """Set really really small coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        [molecule.rezero(epsilon=epsilon) for molecule in self]

    def rotate(self, angle=None, axis=None, anchor_point=None,
               rot_point=None, from_vector=None, to_vector=None,
               degrees=False, transform_matrix=None, verbose=False, **kwargs):
        """Rotate `Molecule` position vectors.

        Parameters
        ----------
        angle : float
        axis : :class:`~sknano.core.math.Vector`, optional
        anchor_point : :class:`~sknano.core.math.Point`, optional
        rot_point : :class:`~sknano.core.math.Point`, optional
        from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
        degrees : bool, optional
        transform_matrix : :class:`~numpy:numpy.ndarray`

        """
        if transform_matrix is None:
            transform_matrix = \
                transformation_matrix(angle=angle, axis=axis,
                                      anchor_point=anchor_point,
                                      rot_point=rot_point,
                                      from_vector=from_vector,
                                      to_vector=to_vector, degrees=degrees,
                                      verbose=verbose, **kwargs)
        [molecule.rotate(transform_matrix=transform_matrix)
         for molecule in self]

    def translate(self, t, fix_anchor_points=True):
        """Translate `Molecule` position vectors by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_points : bool, optional

        """
        [molecule.translate(t, fix_anchor_point=fix_anchor_points)
         for molecule in self]
