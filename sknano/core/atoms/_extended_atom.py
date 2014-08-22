# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class with extended feature set (:mod:`sknano.core.atoms._extended_atom`)
===============================================================================

An "eXtended" `Atom` class for structure analysis.

.. currentmodule:: sknano.core.atoms._extended_atom

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core.math import Vector
from ._atom import Atom

__all__ = ['XAtom']


class XAtom(Atom):
    """An eXtended `Atom` class for structure analysis.

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
    CN : int, optional
        `XAtom` coordination number.
    NN : sequence, optional
        List of nearest-neighbor `XAtom` objects instances

    """

    def __init__(self, element=None, atomID=0, moleculeID=0, atomtype=1, q=0.,
                 mass=None, x=None, y=None, z=None, vx=None, vy=None, vz=None,
                 nx=None, ny=None, nz=None, CN=None, NN=None):

        super(XAtom, self).__init__(element=element, m=mass, x=x, y=y, z=z)

        self._v = Vector([vx, vy, vz])

        self._n = np.zeros(3, dtype=int)
        for i, ni in enumerate((nx, ny, nz)):
            if ni is not None:
                self._n[i] = ni

        self._atomID = int(atomID)
        self._moleculeID = int(moleculeID)
        self._atomtype = int(atomtype)
        self._q = q
        self._CN = CN
        self._NN = NN

        self._attributes.extend(
            ['q','v', 'atomID', 'moleculeID', 'atomtype', 'CN', 'NN'])

    @property
    def CN(self):
        """Return `XAtom` coordination number."""
        return self._CN

    @CN.setter
    def CN(self, value=int):
        """Set `XAtom` coordination number."""
        self._CN = int(value)

    @property
    def NN(self):
        """Return list of nearest-neighbor `XAtom` objects."""
        return self._NN

    @NN.setter
    def NN(self, value=list):
        """Set list of nearest-neighbor `XAtom` objects."""
        self._NN = value

    @property
    def q(self):
        """Charge :math:`q` as multiple of elementary charge :math:`e`.

        """
        return self._q

    @q.setter
    def q(self, value=float):
        """Set `XAtom` charge :math:`q`.

        Parameters
        ----------
        value : {int, float}
            net charge on atom as a multiple of the elementary charge
            :math:`e`.

        """
        self._q = value

    @property
    def atomID(self):
        """:attr:`~XAtom.atomID`."""
        return self._atomID

    @atomID.setter
    def atomID(self, value=int):
        """Set atom ID of atom.

        Parameters
        ----------
        value : int
            atom ID

        """
        self._atomID = int(value)

    @property
    def moleculeID(self):
        """:attr:`~XAtom.moleculeID`."""
        return self._moleculeID

    @moleculeID.setter
    def moleculeID(self, value=int):
        """Set :attr:`~XAtom.moleculeID` attribute.

        Parameters
        ----------
        value : int
            molecule ID

        """
        self._moleculeID = int(value)

    @property
    def atomtype(self):
        """:attr:`~XAtom.atomtype`."""
        return self._atomtype

    @atomtype.setter
    def atomtype(self, value=int):
        """Set :attr:`~XAtom.atomtype` attribute.

        Parameters
        ----------
        value : int
            atom type

        """
        self._atomtype = int(value)

    @property
    def nx(self):
        """:math:`n_x` image flag.

        Returns
        -------
        int
            :math:`n_x` image flag.

        """
        return self._n[0]

    @nx.setter
    def nx(self, value=int):
        """Set :math:`n_x` image flag.

        Parameters
        ----------
        value : float
            :math:`n_x` image flag.

        """
        self._n[0] = int(value)

    @property
    def ny(self):
        """:math:`n_y` image flag.

        Returns
        -------
        int
            :math:`n_y` image flag.

        """
        return self._n[1]

    @ny.setter
    def ny(self, value=int):
        """Set :math:`n_y` image flag.

        Parameters
        ----------
        value : float
            :math:`n_y` image flag.

        """
        self._n[1] = int(value)

    @property
    def nz(self):
        """:math:`n_z` image flag.

        Returns
        -------
        int
            :math:`n_z` image flag.

        """
        return self._n[2]

    @nz.setter
    def nz(self, value=int):
        """Set `XAtom` :math:`n_z` image flag.

        Parameters
        ----------
        value : float
            :math:`n_z` image flag.

        """
        self._n[2] = int(value)

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
    def n(self, value=np.ndarray):
        """Set :math:`n_x, n_y, n_z` image flags of `XAtom`.

        Parameters
        ----------
        value : array_like
            3-element ndarray of [:math:`n_x`, :math:`n_y`, :math:`n_z`]
            image flags of `XAtom`.

        """
        self._n[:] = value

    @property
    def vx(self):
        """:math:`x` component of `XAtom` velocity vector"""
        return self._v.x

    @vx.setter
    def vx(self, value=float):
        """Set :math:`v_x`.

        Set :math:`v_x`, the :math:`x` component of `XAtom` velocity
        vector.

        Parameters
        ----------
        value : float
            :math:`v_x` component of velocity

        """
        self._v.x = value

    @property
    def vy(self):
        """:math:`x` component of `XAtom` velocity vector"""
        return self._v.y

    @vy.setter
    def vy(self, value=float):
        """Set :math:`v_y`.

        Set :math:`v_y`, the :math:`y` component of `XAtom` velocity
        vector.

        Parameters
        ----------
        value : float
            :math:`v_y` component of velocity

        """
        self._v.y = value

    @property
    def vz(self):
        """:math:`z` component of `XAtom` velocity vector"""
        return self._v.z

    @vz.setter
    def vz(self, value=float):
        """Set :math:`v_z`.

        Set :math:`v_z`, the :math:`z` component of `XAtom` velocity
        vector.

        Parameters
        ----------
        value : float
            :math:`v_z` component of velocity

        """
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
    def v(self, value=np.ndarray):
        """Set :math:`x, y, z` components of `XAtom` velocity.

        Parameters
        ----------
        value : array_like
            3-element ndarray of [:math:`v_x`, :math:`v_y`, :math:`v_z`]
            velocity components of `XAtom`.

        """
        self._v[:] = value

    def compute_pyramidalization_angle(self):
        pass
