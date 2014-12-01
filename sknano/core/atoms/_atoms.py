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

        Computes the position vector of the center-of-mass coordinates:

        .. math::

           \\mathbf{R}_{CM} = \\frac{1}{M}\\sum_{i=1}^{N_{\\mathrm{atoms}}}
           m_i\\mathbf{r}_i

        Returns
        -------
        CM : :class:`~sknano.core.math.Vector`
            The position vector of the center of mass coordinates.

        """
        masses = np.asarray([self.masses])
        coords = self.coords
        MxR = masses.T * coords
        CM = Vector(np.sum(MxR, axis=0) / np.sum(masses))
        CM.rezero()
        return CM

    @property
    def M(self):
        """Total mass of `Atoms`."""
        #return math.fsum(self.masses)
        return self.masses.sum()

    @property
    def centroid(self):
        """Centroid of `Atoms`.

        Computes the position vector of the centroid of the `Atoms`
        coordinates.

        .. math::
           \\mathbf{C} =
           \\frac{\\sum_{i=1}^{N_{\\mathrm{atoms}}}
           m_i\\mathbf{r}_i}{\\sum_{i=1}^{N_{\\mathrm{atoms}}}m_i}

        Returns
        -------
        C : `~sknano.core.math.Vector`
            The position vector of the centroid coordinates.
        """
        C = Vector(np.mean(self.coords, axis=0))
        C.rezero()
        return C

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
    def inertia_tensor(self):
        """Return the inertia tensor."""
        Ixx = (self.masses * (self.y**2 + self.z**2)).sum()
        Iyy = (self.masses * (self.x**2 + self.z**2)).sum()
        Izz = (self.masses * (self.x**2 + self.y**2)).sum()
        Ixy = Iyx = (-self.masses * self.x * self.y).sum()
        Ixz = Izx = (-self.masses * self.x * self.z).sum()
        Iyz = Izy = (-self.masses * self.y * self.z).sum()
        return np.array([[Ixx, Ixy, Ixz], [Iyx, Iyy, Iyz], [Izx, Izy, Izz]])

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
        """Alias for :meth:`Atoms.rezero`."""
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
               rot_point=None, from_vector=None, to_vector=None,
               deg2rad=False, transform_matrix=None, verbose=False):
        """Rotate `Atom` position vectors.

        Parameters
        ----------
        angle : float
        rot_axis : :class:`~sknano.core.math.Vector`, optional
        anchor_point : :class:`~sknano.core.math.Point`, optional
        rot_point : :class:`~sknano.core.math.Point`, optional
        from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
        deg2rad : bool, optional
        transform_matrix : :class:`~numpy:numpy.ndarray`

        """
        if transform_matrix is None:
            transform_matrix = \
                transformation_matrix(angle=angle, rot_axis=rot_axis,
                                      anchor_point=anchor_point,
                                      rot_point=rot_point,
                                      from_vector=from_vector,
                                      to_vector=to_vector, deg2rad=deg2rad,
                                      verbose=verbose)
        [atom.rotate(transform_matrix=transform_matrix) for atom in self]

    def translate(self, t, fix_anchor_points=True):
        """Translate `Atom` position vectors by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_points : bool, optional

        """
        [atom.translate(t, fix_anchor_point=fix_anchor_points)
         for atom in self]
