# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class with extended feature set (:mod:`sknano.core.atoms._extended_atom`)
===============================================================================

An "eXtended" `Atom` class for structure analysis.

.. currentmodule:: sknano.core.atoms._extended_atom

"""
from __future__ import absolute_import, division, print_function
from six.moves import zip

__docformat__ = 'restructuredtext en'

from collections import OrderedDict
from functools import total_ordering
import numbers
import numpy as np

from sknano.core import xyz
from sknano.core.math import Point, Vector
from ._atom import Atom

__all__ = ['XAtom']


@total_ordering
class XAtom(Atom):
    """An `Atom` class with an eXtended set of attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    x, y, z : float, optional
        :math:`x, y, z` components of `XAtom` position vector relative to
        origin.
    atomID : int, optional
        atom ID
    moleculeID : int, optional
        molecule ID
    atomtype : int, optional
        atom type
    q : {int, float}, optional
        Net charge of `XAtom`.
    vx, vy, vz : float, optional
        :math:`v_x, v_y, v_z` components of `XAtom` velocity.
    nx, ny, nz : int, optional
        :math:`n_x, n_y, n_z` image flags

    """
    _atomattrs = Atom._atomattrs + \
        ['atomID', 'moleculeID', 'atomtype', 'q', 'dr',
         'r', 'x', 'y', 'z', 'v', 'vx', 'vy', 'vz',
         'f', 'fx', 'fy', 'fz', 'n', 'nx', 'ny', 'nz',
         'pe', 'ke', 'etotal']

    __hash__ = Atom.__hash__

    def __init__(self, element=None, atomID=0, moleculeID=0, atomtype=1,
                 q=0., mass=None, x=None, y=None, z=None,
                 vx=None, vy=None, vz=None, fx=None, fy=None, fz=None,
                 nx=None, ny=None, nz=None, pe=None, ke=None, etotal=None,
                 **kwargs):
        if mass is None and 'm' in kwargs:
            mass = kwargs['m']

        super(XAtom, self).__init__(element=element, mass=mass)

        self.r = Vector([x, y, z])
        #self._p0 = Point([x, y, z])
        self.dr = Vector(np.zeros(3), p0=[x, y, z])

        self._v = Vector([vx, vy, vz])
        self._f = Vector([fx, fy, fz])
        self._n = Point([nx, ny, nz], dtype=int)

        self._atomID = int(atomID)
        self._moleculeID = int(moleculeID)
        self._atomtype = int(atomtype)
        self._q = q

        self.pe = pe
        self.ke = ke
        self.etotal = etotal

    def __str__(self):
        """Return nice string representation of `XAtom`."""
        strrep = "Atom(element={element!r}, atomID={atomID!r}, " + \
            "moleculeID={moleculeID!r}, atomtype={atomtype!r})"
        parameters = dict(element=self.element, atomID=self.atomID,
                          moleculeID=self.moleculeID, atomtype=self.atomtype)

        return strrep.format(**parameters)

    def __repr__(self):
        """Return canonical string representation of `XAtom`."""
        reprstr = "Atom(element={element!r}, atomID={atomID!r}, " + \
            "moleculeID={moleculeID!r}, atomtype={atomtype!r}, " + \
            "q={q!r}, mass={mass!r}, x={x:.6f}, y={y:.6f}, z={z:.6f}, " \
            "pe={pe!r}, ke={ke!r}, etotal={etotal!r})"
        parameters = dict(element=self.element, atomID=self.atomID,
                          moleculeID=self.moleculeID, atomtype=self.atomtype,
                          q=self.q, mass=self.mass,
                          x=self.x, y=self.y, z=self.z,
                          pe=self.pe, ke=self.ke, etotal=self.etotal)

        return reprstr.format(**parameters)

    def __eq__(self, other):
        return super(XAtom, self).__eq__(other)

    def __lt__(self, other):
        if self.atomID < other.atomID:
            return True
        else:
            return super(XAtom, self).__lt__(other)

    @property
    def q(self):
        """Charge :math:`q` as multiple of elementary charge :math:`e`.

        """
        return self._q

    @q.setter
    def q(self, value):
        """Set `XAtom` charge :math:`q`.

        Parameters
        ----------
        value : {int, float}
            net charge on atom as a multiple of the elementary charge
            :math:`e`.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._q = value

    @property
    def atomID(self):
        """:attr:`~XAtom.atomID`."""
        return self._atomID

    @atomID.setter
    def atomID(self, value):
        """Set atom ID of atom.

        Parameters
        ----------
        value : int
            atom ID

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._atomID = int(value)

    @property
    def moleculeID(self):
        """:attr:`~XAtom.moleculeID`."""
        return self._moleculeID

    @moleculeID.setter
    def moleculeID(self, value):
        """Set :attr:`~XAtom.moleculeID` attribute.

        Parameters
        ----------
        value : int
            molecule ID

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._moleculeID = int(value)

    @property
    def atomtype(self):
        """:attr:`~XAtom.atomtype`."""
        return self._atomtype

    @atomtype.setter
    def atomtype(self, value):
        """Set :attr:`~XAtom.atomtype` attribute.

        Parameters
        ----------
        value : int
            atom type

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._atomtype = int(value)

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
        self._r = Vector(value)

    @property
    def dr(self):
        """:math:`x, y, z` components of `Atom` displacement vector.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x, y, z`] coordinates of `Atom`.

        """
        return self._dr

    @dr.setter
    def dr(self, value):
        """Set :math:`x, y, z` components of `Atom` displacement vector.

        Parameters
        ----------
        value : array_like
            :math:`x, y, z` components of `Atom` displacement vector.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._dr = Vector(value)

    def get_coords(self, asdict=False):
        """Return atom coords.

        Parameters
        ----------
        asdict : bool, optional

        Returns
        -------
        coords : :class:`python:~collections.OrderedDict` or ndarray

        """
        if asdict:
            return OrderedDict(list(zip(xyz, self.r)))
        else:
            return self.r

    @property
    def nx(self):
        """:math:`n_x` image flag.

        Returns
        -------
        int
            :math:`n_x` image flag.

        """
        return self._n.x

    @nx.setter
    def nx(self, value):
        """Set :math:`n_x` image flag.

        Parameters
        ----------
        value : float
            :math:`n_x` image flag.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._n.x = int(value)

    @property
    def ny(self):
        """:math:`n_y` image flag.

        Returns
        -------
        int
            :math:`n_y` image flag.

        """
        return self._n.y

    @ny.setter
    def ny(self, value):
        """Set :math:`n_y` image flag.

        Parameters
        ----------
        value : float
            :math:`n_y` image flag.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._n.y = int(value)

    @property
    def nz(self):
        """:math:`n_z` image flag.

        Returns
        -------
        int
            :math:`n_z` image flag.

        """
        return self._n.z

    @nz.setter
    def nz(self, value):
        """Set `XAtom` :math:`n_z` image flag.

        Parameters
        ----------
        value : float
            :math:`n_z` image flag.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._n.z = int(value)

    @property
    def n(self):
        """:math:`n_x, n_y, n_z` image flags of `XAtom`.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`n_x`, :math:`n_y`, :math:`n_z`]
            image flags of `XAtom`.

        """
        return self._n

    @n.setter
    def n(self, value):
        """Set :math:`n_x, n_y, n_z` image flags of `XAtom`.

        Parameters
        ----------
        value : array_like
            3-element ndarray of [:math:`n_x`, :math:`n_y`, :math:`n_z`]
            image flags of `XAtom`.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._n = Point(value, dtype=int)

    @property
    def vx(self):
        """:math:`x` component of `XAtom` velocity vector"""
        return self._v.x

    @vx.setter
    def vx(self, value):
        """Set :math:`v_x`.

        Set :math:`v_x`, the :math:`x` component of `XAtom` velocity
        vector.

        Parameters
        ----------
        value : float
            :math:`v_x` component of velocity

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._v.x = value

    @property
    def vy(self):
        """:math:`x` component of `XAtom` velocity vector"""
        return self._v.y

    @vy.setter
    def vy(self, value):
        """Set :math:`v_y`.

        Set :math:`v_y`, the :math:`y` component of `XAtom` velocity
        vector.

        Parameters
        ----------
        value : float
            :math:`v_y` component of velocity

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._v.y = value

    @property
    def vz(self):
        """:math:`z` component of `XAtom` velocity vector"""
        return self._v.z

    @vz.setter
    def vz(self, value):
        """Set :math:`v_z`.

        Set :math:`v_z`, the :math:`z` component of `XAtom` velocity
        vector.

        Parameters
        ----------
        value : float
            :math:`v_z` component of velocity

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._v.z = value

    @property
    def v(self):
        """:math:`v_x, v_y, v_z` array of velocity components.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`v_x`, :math:`v_y`, :math:`v_z`]
            velocity components of `XAtom`.

        """
        return self._v

    @v.setter
    def v(self, value):
        """Set :math:`x, y, z` components of `XAtom` velocity.

        Parameters
        ----------
        value : array_like
            3-element ndarray of [:math:`v_x`, :math:`v_y`, :math:`v_z`]
            velocity components of `XAtom`.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._v = Vector(value)

    @property
    def fx(self):
        """:math:`x` component of `XAtom` force vector"""
        return self._f.x

    @fx.setter
    def fx(self, value):
        """Set :math:`f_x`.

        Set :math:`f_x`, the :math:`x` component of `XAtom` force vector.

        Parameters
        ----------
        value : float
            :math:`f_x` component of force vector.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._f.x = value

    @property
    def fy(self):
        """:math:`x` component of `XAtom` force vector"""
        return self._f.y

    @fy.setter
    def fy(self, value):
        """Set :math:`f_y`.

        Set :math:`f_y`, the :math:`y` component of `XAtom` force vector.

        Parameters
        ----------
        value : float
            :math:`f_y` component of force vector.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._f.y = value

    @property
    def fz(self):
        """:math:`z` component of `XAtom` force vector"""
        return self._f.z

    @fz.setter
    def fz(self, value):
        """Set :math:`f_z`.

        Set :math:`f_z`, the :math:`z` component of `XAtom` force vector.

        Parameters
        ----------
        value : float
            :math:`f_z` component of force vector.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._f.z = value

    @property
    def f(self):
        """:math:`f_x, f_y, f_z` array of force vector components.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`f_x`, :math:`f_y`, :math:`f_z`]
            force vector components of `XAtom`.

        """
        return self._f

    @f.setter
    def f(self, value):
        """Set :math:`x, y, z` components of `XAtom` force vector.

        Parameters
        ----------
        value : array_like
            3-element ndarray of [:math:`f_x`, :math:`f_y`, :math:`f_z`]
            force components of `XAtom`.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._f = Vector(value)

    @property
    def pe(self):
        """`Atom` potential energy."""
        return self._pe

    @pe.setter
    def pe(self, value):
        if not isinstance(value, (float, type(None))):
            raise TypeError('Expected a number')
        if value is None:
            value = 0.0
        self._pe = value

    @property
    def ke(self):
        """`Atom` kinetic energy."""
        return self._ke

    @ke.setter
    def ke(self, value):
        if not isinstance(value, (float, type(None))):
            raise TypeError('Expected a number')
        if value is None:
            value = 0.0
        self._ke = value

    @property
    def etotal(self):
        """`Atom` total energy (:attr:`~XAtoms.pe` + :attr:`~XAtoms.ke`)."""
        return self._etotal

    @etotal.setter
    def etotal(self, value):
        if not isinstance(value, (float, type(None))):
            raise TypeError('Expected a number')
        if value is None:
            try:
                value = self.pe + self.ke
            except TypeError:
                value = 0.0
        self._etotal = value

    def rezero(self, epsilon=1.0e-10):
        """Re-zero position vector components.

        Set position vector components with absolute value less than
        `epsilon` to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        self._r.rezero(epsilon=epsilon)

    def rezero_coords(self, epsilon=1.0e-10):
        """Alias for :meth:`Atom.rezero`."""
        self.rezero(epsilon=epsilon)

    def rotate(self, angle=None, rot_axis=None, anchor_point=None,
               rot_point=None, from_vector=None, to_vector=None, deg2rad=False,
               transform_matrix=None, verbose=False):
        """Rotate `Atom` position vector.

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
        self.r.rotate(angle=angle, rot_axis=rot_axis,
                      anchor_point=anchor_point,
                      rot_point=rot_point, from_vector=from_vector,
                      to_vector=to_vector, transform_matrix=transform_matrix,
                      deg2rad=deg2rad, verbose=verbose)

    def translate(self, t, fix_anchor_point=True):
        """Translate `Atom` position vector by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_point : bool, optional

        """
        #TODO compare timing benchmarks for translation of position vector.
        self.r.translate(t, fix_anchor_point=fix_anchor_point)
        #self.r += t

    def todict(self):
        return dict(element=self.element, atomID=self.atomID,
                    moleculeID=self.moleculeID, atomtype=self.atomtype,
                    q=self.q, mass=self.mass,
                    x=self.x, y=self.y, z=self.z,
                    vx=self.vx, vy=self.vy, vz=self.vz,
                    fx=self.fx, fy=self.fy, fz=self.fz,
                    nx=self.nx, ny=self.ny, nz=self.nz,
                    pe=self.pe, ke=self.ke, etotal=self.etotal)
