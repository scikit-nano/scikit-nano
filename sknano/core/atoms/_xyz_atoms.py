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

import numbers

from collections import OrderedDict
from math import fsum
from operator import attrgetter

import numpy as np

from sknano.core import rezero_array, xyz
from sknano.core.math import Vector, Vectors

from ._atoms import Atom, Atoms

__all__ = ['XYZAtom', 'XYZAtoms']


class XYZAtom(Atom):
    """An `Atom` sub-class with x, y, z attributes.

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
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.r == other.r and super().__eq__(other)

    def __le__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.r > other.r or not super().__le__(other):
            return False
        return True

    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.r >= other.r or not super().__lt__(other):
            return False
        return True

    def __ge__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.r < other.r or not super().__ge__(other):
            return False
        return True

    def __gt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.r <= other.r or not super().__gt__(other):
            return False
        return True

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
        :class:`Vector`
            3D :class:`Vector` of [:math:`x, y, z`] coordinates of `Atom`.

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

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
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
    def coords(self):
        """Alias for :attr:`Atoms.r`."""
        return self.r

    @property
    def positions(self):
        """Alias for :attr:`Atoms.r`."""
        return self.r

    @property
    def r(self):
        """:class:`Vectors` of :attr:`Atom.r` position `Vector`\ s"""
        return Vectors([atom.r for atom in self])

    @property
    def dr(self):
        """:class:`Vectors` of :attr:`Atom.dr` displacement `Vector`\ s"""
        return Vectors([atom.dr for atom in self])

    @property
    def x(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`XYZAtom.x` coordinates."""
        # return np.asarray(self.r)[:, 0]
        return self.r.x

    @property
    def y(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`XYZAtom.y` coordinates."""
        # return np.asarray(self.r)[:, 1]
        return self.r.y

    @property
    def z(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`XYZAtom.z` coordinates."""
        # return np.asarray(self.r)[:, 2]
        return self.r.z

    @property
    def inertia_tensor(self):
        """Return inertia tensor about the origin."""
        masses = self.masses
        r = self.r - self.center_of_mass
        x, y, z = r.x, r.y, r.z

        Ixx = fsum(masses * (y**2 + z**2))
        Iyy = fsum(masses * (x**2 + z**2))
        Izz = fsum(masses * (x**2 + y**2))
        Ixy = Iyx = fsum(-masses * x * y)
        Ixz = Izx = fsum(-masses * x * z)
        Iyz = Izy = fsum(-masses * y * z)
        I = np.array([[Ixx, Ixy, Ixz],
                      [Iyx, Iyy, Iyz],
                      [Izx, Izy, Izz]])
        return rezero_array(I, epsilon=1e-10)

    @property
    def moment_of_inertia(self):
        """Alias for :attr:`~XYZAtoms.inertia_tensor`."""
        return self.inertia_tensor

    @property
    def principal_axes(self):
        """Return principal axes of rotation computed from the inertia tensor.

        Since the :attr:`~XYZAtoms.inertia_tensor` forms a :math:`3\\times3`
        real symmetric matrix :math:`[I]`, then there is a matrix
        :math:`[P]=[\\mathbf{p}_1\\,\\mathbf{p}_2\\,\\mathbf{p}_3]`
        that orthogonally diagonalizes :math:`[I]` such that
        :math:`[\\Lambda]=[P]^T [I] [P]` is a diagonal matrix:

        .. math::

           [\\Lambda] =  \\begin{bmatrix}
           \\lambda_1 & 0 & 0\\\\
           0 & \\lambda_2 & 0\\\\
           0 & 0 & \\lambda_3
           \\end{bmatrix}

        where :math:`\\lambda_1,\\lambda_2,\\lambda_3` are the eigenvalues
        of :math:`[I]` corresponding to the eigenvectors
        :math:`\\mathbf{p}_1, \\mathbf{p}_2, \\mathbf{p}_3`, which form the
        column vectors of matrix :math:`[P]` and are the principal axes
        of the moment of inertia. In other words, the matrix :math:`[P]`
        is a rotation matrix which transforms the coordinates
        :math:`\\mathbf{r}`
        of the atoms such that the moment of inertia i

        Returns
        -------
        p : :class:`~sknano.core.math.Vectors`
            :class:`~sknano.core.math.Vectors` object where `p[0]`, `p[1]`,
            and `p[2]` are the 3 principal axes corresponding to the
            eigenvectors that form the successive columns of :math:`P`.

        """
        w, v = np.linalg.eig(self.inertia_tensor)
        v = rezero_array(v[:, np.argsort(w)[::-1]], epsilon=1e-10)
        # v = rezero_array(v, epsilon=1e-10)
        p = Vectors(list(map(Vector, [v[:, i] for i in range(3)])))
        return p

    @property
    def principal_moments_of_inertia(self):
        """Return principal moments of inertia."""
        w = np.linalg.eigvals(self.inertia_tensor)
        return w[np.argsort(w)[::-1]]

    @property
    def radius_of_gyration(self):
        """Return radius of gyration.

        .. math::

           R_g = \\sqrt{\\frac{1}{M}\\sum_{i=1}^{N}
           m_i(\\mathbf{r}_i-\\mathbf{R})^2}

        """
        r = self.r
        masses = self.masses
        R = self.center_of_mass
        M = self.M
        return np.sqrt(np.sum(masses * (r - R).dot(r - R)) / M)

    def align_principal_axis(self, index, vector):
        """Align :attr:`~XYZAtoms.principal_axes`[`index`] along `vector`.

        Parameters
        ----------
        index : :class:`~python:int`
        vector : {:class:`~python:list`, :class:`~numpy:numpy.ndarray`}

        Examples
        --------
        >>> from sknano.generators import SWNTGenerator
        >>> swnt = SWNTGenerator(10, 5)

        """
        self.rotate(from_vector=self.principal_axes[index],
                    to_vector=Vector(vector))

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
        region : :class:`~sknano.core.geometric_regions.GeometricRegion`

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
