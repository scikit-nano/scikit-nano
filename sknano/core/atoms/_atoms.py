# -*- coding: utf-8 -*-
"""
==============================================================================
Base class for structure data atoms (:mod:`sknano.core.atoms._atoms`)
==============================================================================

.. currentmodule:: sknano.core.atoms._atoms

"""
from __future__ import absolute_import, division, print_function
from six.moves import zip
__docformat__ = 'restructuredtext en'

from collections import OrderedDict
from operator import attrgetter

import numpy as np

from sknano.core import UserList, xyz
from sknano.core.math import Vector, transformation_matrix
from sknano.utils.geometric_shapes import Cuboid  # , Rectangle

__all__ = ['Atoms']


class Atoms(UserList):
    """Base class for collection of `Atom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instance objects or an
        existing `Atoms` instance object.
    copylist : bool, optional
        perform shallow copy of atoms list
    deepcopy : bool, optional
        perform deepcopy of atoms list


    """
    _atomattrs = ['symbol', 'Z', 'm', 'r', 'x', 'y', 'z', 'dr']

    def __init__(self, atoms=None, copylist=True, deepcopy=False):
        super(Atoms, self).__init__(initlist=atoms,
                                    copylist=copylist,
                                    deepcopy=deepcopy)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        """Return canonical string representation of `Atoms`."""
        return "Atoms(atoms={!r})".format(self.data)

    def sort(self, key=None, reverse=False):

        if key is None:
            self.data.sort(key=attrgetter('element', 'Z', 'z'),
                           reverse=reverse)
        else:
            self.data.sort(key=key, reverse=reverse)

    @property
    def Natoms(self):
        """Number of atoms in `Atoms`."""
        return len(self)

    @property
    def CM(self):
        """Center-of-Mass coordinates of `Atoms`.

        Returns
        -------
        ndarray
            3-element ndarray specifying center-of-mass coordinates of `Atoms`.

        """
        masses = np.asarray([self.masses])
        coords = self.coords
        MxR = masses.T * coords
        return Vector(np.sum(MxR, axis=0) / np.sum(masses))

    @property
    def M(self):
        """Total mass of `Atoms`."""
        #return math.fsum(self.masses)
        return self.masses.sum()

    @property
    def coords(self):
        """Return list of `Atom` coordinates."""
        return np.asarray([atom.r for atom in self])

    @property
    def masses(self):
        """Return list of `Atom` masses."""
        return np.asarray([atom.m for atom in self])

    @property
    def symbols(self):
        """Return list of `Atom` symbols."""
        return np.asarray([atom.symbol for atom in self])

    @property
    def x(self):
        """Return :math:`x` coordinates of `Atom` objects as array."""
        return self.coords[:,0]

    @property
    def y(self):
        """Return :math:`y` coordinates of `Atom` objects as array."""
        return self.coords[:,1]

    @property
    def z(self):
        """Return :math:`z` coordinates of `Atom` objects as array."""
        return self.coords[:,2]

    @property
    def bounds(self):
        """Return bounds of `Atoms`."""
        return Cuboid(pmin=[self.x.min(), self.y.min(), self.z.min()],
                      pmax=[self.x.max(), self.y.max(), self.z.max()])

    def center_CM(self, axes=None):
        """Center atoms on CM coordinates."""
        dr = -self.CM
        self.translate(dr)

    def clip_bounds(self, region, center_before_clipping=False):
        """Remove atoms outside the given limits along given dimension.

        Parameters
        ----------
        region : :class:`~sknano.utils.geometric_shapes.`GeometricRegion`

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

        #for dim, limits in region.limits.iteritems():
        #    atoms = atoms[np.where(getattr(self, dim) <= limits['max'])]
        #    atoms = atoms[np.where(getattr(self, dim) >= limits['min'])]
        #    self = atoms.tolist()

        if CM0 is not None:
            self.translate(CM0)

    def filter(self, condition, invert=False):
        """Filter `Atoms` by `condition`.

        Parameters
        ----------
        condition : array_like, bool
            Boolean index array having same shape as the initial dimensions
            of the list of `Atoms` being indexed.
        invert : bool, optional
            If `True`, the boolean array `condition` is inverted element-wise.

        Returns
        -------
        filtered_atoms : `Atoms`
            If `invert` is `False`, return the elements where `condition`
            is `True`.

            If `invert` is `True`, return the elements where `~condition`
            (i.e., numpy.invert(condition)) is `True`.

        """
        if invert:
            condition = ~condition
        return self.__class__(atoms=np.asarray(self)[condition].tolist())

    def get_atoms(self, asarray=False):
        """Return list of `Atoms`.

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
        """Return atom coords.

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
        """Alias for :meth:`Atoms.rezero`."""
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
        [atom.rezero(epsilon=epsilon) for atom in self]

    def rotate(self, angle=None, rot_axis=None, anchor_point=None,
               deg2rad=False, transform_matrix=None):
        """Rotate atom coordinates about arbitrary axis.

        Parameters
        ----------
        angle : float

        """
        if transform_matrix is None:
            transform_matrix = \
                transformation_matrix(angle, rot_axis=rot_axis,
                                      anchor_point=anchor_point,
                                      deg2rad=deg2rad)
        [atom.rotate(transform_matrix=transform_matrix) for atom in self]

    def translate(self, t):
        """Translate atom coordinates.

        Parameters
        ----------
        t : array_like
            3-elment array of :math:`x,y,z` components of translation vector
        """
        [atom.translate(t) for atom in self]
