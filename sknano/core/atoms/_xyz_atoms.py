# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes with x, y, z attributes (:mod:`sknano.core.atoms._xyz_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._xyz_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from collections import OrderedDict
from functools import total_ordering
from operator import attrgetter
import numbers
import numpy as np

from sknano.core import xyz
from sknano.core.math import Vector
from sknano.core.geometric_regions import Cuboid  # , Rectangle

from ._atoms import Atom, Atoms

__all__ = ['XYZAtom', 'XYZAtoms']


@total_ordering
class XYZAtom(Atom):
    """An `Atom` class with x, y, z attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    x, y, z : float, optional
        :math:`x, y, z` components of `XYZAtom` position vector relative to
        origin.

    """
    def __init__(self, *args, x=None, y=None, z=None, **kwargs):
        super().__init__(*args, **kwargs)

        self._r0 = Vector([x, y, z])
        self._r = Vector([x, y, z])
        self.fmtstr = super().fmtstr + ", x={x:.6f}, y={y:.6f}, z={z:.6f}"

    def __eq__(self, other):
        return self.r == other.r and super().__eq__(other)

    def __lt__(self, other):
        return (self.r < other.r and super().__le__(other)) or \
            (self.r <= other.r and super().__lt__(other))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['x', 'y', 'z'])
        return attrs

    @property
    def x(self):
        """:math:`x`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`x`-coordinate in units of **Angstroms**.

        """
        return self._r.x

    @x.setter
    def x(self, value):
        """Set `Atom` :math:`x`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`x`-coordinate in units of **Angstroms**.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._r.x = value

    @property
    def y(self):
        """:math:`y`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`y`-coordinate in units of **Angstroms**.

        """
        return self._r.y

    @y.setter
    def y(self, value):
        """Set `Atom` :math:`y`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`y`-coordinate in units of **Angstroms**.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._r.y = value

    @property
    def z(self):
        """:math:`z`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`z`-coordinate in units of **Angstroms**.

        """
        return self._r.z

    @z.setter
    def z(self, value):
        """Set `Atom` :math:`z`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`z`-coordinate in units of **Angstroms**.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._r.z = value

    @property
    def r(self):
        """:math:`x, y, z` components of `Atom` position vector.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x, y, z`] coordinates of `Atom`.

        """
        return self._r

    @r.setter
    def r(self, value):
        """Set :math:`x, y, z` components of `Atom` position vector.

        Parameters
        ----------
        value : array_like
            :math:`x, y, z` coordinates of `Atom` position vector relative to
            the origin.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._r[:] = Vector(value, nd=3)

    @property
    def r0(self):
        """:math:`x, y, z` components of `Atom` position vector at t=0.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x, y, z`] coordinates of `Atom`.

        """
        return self._r0

    @r0.setter
    def r0(self, value):
        """Set :math:`x, y, z` components of `Atom` position vector at t=0.

        Parameters
        ----------
        value : array_like
            :math:`x, y, z` coordinates of `Atom` position vector relative to
            the origin.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._r0[:] = Vector(value, nd=3)

    @property
    def dr(self):
        """:math:`x, y, z` components of `Atom` displacement vector.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x, y, z`] coordinates of `Atom`.

        """
        return self.r - self.r0

    def get_coords(self, asdict=False):
        """Return atom coords.

        Parameters
        ----------
        asdict : bool, optional

        Returns
        -------
        coords : :class:`~python:collections.OrderedDict` or ndarray

        """
        if asdict:
            return OrderedDict(list(zip(xyz, self.r)))
        else:
            return self.r

    def rezero(self, epsilon=1.0e-10):
        """Alias for :meth:`Atom.rezero_xyz`, but calls `super` class \
            `rezero` method as well."""
        self.r.rezero(epsilon)
        super().rezero(epsilon)

    def rezero_coords(self, epsilon=1.0e-10):
        """Alias for :meth:`Atom.rezero_xyz`."""
        self.rezero_xyz(epsilon=epsilon)

    def rezero_xyz(self, epsilon=1.0e-10):
        """Re-zero position vector components.

        Set position vector components with absolute value less than
        `epsilon` to zero.

        Unlike the :meth:`Atom.rezero` method, this method does **not**
        call the super class `rezero` method.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        self.r.rezero(epsilon=epsilon)

    def rotate(self, **kwargs):
        """Rotate `Atom` position vector.

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
        self.r.rotate(**kwargs)
        self.r0.rotate(**kwargs)
        super().rotate(**kwargs)

    def translate(self, t, fix_anchor_point=True):
        """Translate `Atom` position vector by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_point : bool, optional

        """
        # TODO compare timing benchmarks for translation of position vector.
        self.r.translate(t, fix_anchor_point=fix_anchor_point)
        self.r0.translate(t, fix_anchor_point=fix_anchor_point)
        super().translate(t, fix_anchor_point=fix_anchor_point)
        # self.r += t

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(x=self.x, y=self.y, z=self.z))
        return super_dict


class XYZAtoms(Atoms):
    """An `Atoms` sub-class for `XYZAtom`\ s.

    A container class for :class:`~sknano.core.atoms.XYZAtom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `XYZAtoms`}, optional
        if not `None`, then a list of `XYZAtom` instance objects or an
        existing `XYZAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return XYZAtom

    def sort(self, key=attrgetter('r'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def center_of_mass(self):
        """Center-of-Mass coordinates of `Atoms`.

        Computes the position vector of the center-of-mass coordinates:

        .. math::

           \\mathbf{R}_{COM} = \\frac{1}{M}\\sum_{i=1}^{N_{\\mathrm{atoms}}}
           m_i\\mathbf{r}_i

        Returns
        -------
        com : :class:`~sknano.core.math.Vector`
            The position vector of the center of mass coordinates.

        """
        masses = np.asarray([self.masses])
        coords = self.coords
        MxR = masses.T * coords
        com = Vector(np.sum(MxR, axis=0) / np.sum(masses))
        com.rezero()
        return com

    @property
    def com(self):
        """Alias for :attr:`~XYZAtoms.center_of_mass`."""
        return self.center_of_mass

    @property
    def CM(self):
        """Alias for :attr:`~XYZAtoms.center_of_mass`."""
        return self.center_of_mass

    @property
    def centroid(self):
        """Centroid of `Atoms`.

        Computes the position vector of the centroid of the `Atoms`
        coordinates.

        .. math::
           \\mathbf{C} =
           \\frac{\\sum_{i=1}^{N_{\\mathrm{atoms}}}
           \\mathbf{r}_i}{N_{\\mathrm{atoms}}}

        Returns
        -------
        C : :class:`~sknano.core.math.Vector`
            The position vector of the centroid coordinates.
        """
        C = Vector(np.mean(self.coords, axis=0))
        C.rezero()
        return C

    @property
    def bounds(self):
        """Bounds of `Atoms`.

        Returns
        -------
        :class:`~sknano.core.geometric_regions.Cuboid`"""
        return Cuboid(pmin=[self.x.min(), self.y.min(), self.z.min()],
                      pmax=[self.x.max(), self.y.max(), self.z.max()])

    @property
    def coords(self):
        """Alias for :attr:`Atoms.r`."""
        return self.r

    @property
    def r(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`Atom.r` position \
            `Vector`\ s"""
        return np.asarray([atom.r for atom in self])

    @property
    def dr(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`Atom.dr` displacement \
            `Vector`\ s"""
        return np.asarray([atom.dr for atom in self])

    @property
    def x(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`x` coordinates."""
        return self.r[:, 0]

    @property
    def y(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`y` coordinates."""
        return self.r[:, 1]

    @property
    def z(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`z` coordinates."""
        return self.r[:, 2]

    @property
    def inertia_tensor(self):
        """Inertia tensor."""
        Ixx = (self.masses * (self.y**2 + self.z**2)).sum()
        Iyy = (self.masses * (self.x**2 + self.z**2)).sum()
        Izz = (self.masses * (self.x**2 + self.y**2)).sum()
        Ixy = Iyx = (-self.masses * self.x * self.y).sum()
        Ixz = Izx = (-self.masses * self.x * self.z).sum()
        Iyz = Izy = (-self.masses * self.y * self.z).sum()
        return np.array([[Ixx, Ixy, Ixz], [Iyx, Iyy, Iyz], [Izx, Izy, Izz]])

    def center_center_of_mass(self, axis=None):
        """Center atoms on center-of-mass coordinates."""
        self.translate(-self.center_of_mass)

    def center_com(self, axis=None):
        """Alias for :attr:`~XYZAtoms.center_center_of_mass`."""
        self.center_center_of_mass(axis=axis)

    def center_CM(self, axis=None):
        """Alias for :attr:`~XYZAtoms.center_center_of_mass`."""
        self.center_center_of_mass(axis=axis)

    def center_centroid(self):
        """Center :attr:`~XYZAtoms.centroid` on origin."""
        self.translate(-self.centroid)

    def clip_bounds(self, region, center_before_clipping=False):
        """Remove atoms outside the given region.

        Parameters
        ----------
        region : :class:`~sknano.core.geometric_regions.`GeometricRegion`

        """
        centroid0 = None
        if center_before_clipping:
            centroid0 = self.centroid
            self.translate(-centroid0)

        self.data = [atom for atom in self if region.contains(atom.r)]

        if centroid0 is not None:
            self.translate(centroid0)

    def get_coords(self, asdict=False):
        """Return atom coords.

        Parameters
        ----------
        asdict : :class:`~python:bool`, optional

        Returns
        -------
        coords : :class:`~python:collections.OrderedDict` or \
            :class:`~numpy:numpy.ndarray`

        """
        coords = self.coords
        if asdict:
            return OrderedDict(list(zip(xyz, coords.T)))
        else:
            return coords

    def rezero_coords(self, epsilon=1.0e-10):
        """Alias for :meth:`Atoms.rezero_xyz`."""
        self.rezero(epsilon=epsilon)

    def rezero_xyz(self, epsilon=1.0e-10):
        """Rezero position vector components with absolute value less than \
            `epsilon`.

        Calls the :meth:`XYZAtom.rezero_xyz` method on each `atom` in `self`.

        Parameters
        ----------
        epsilon : :class:`~python:float`
            values with absolute value less than `epsilon` are set to zero.

        """
        [atom.rezero_xyz(epsilon=epsilon) for atom in self]
