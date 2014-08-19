# -*- coding: utf-8 -*-
"""
============================================================================
Structure Atom (:mod:`sknano.io.atoms._structure_atom`)
============================================================================

`Atom` class for :class:`~sknano.generator.StructureGenerator` classes.

.. currentmodule:: sknano.io.atoms._structure_atom

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numpy as np

from ._atom import Atom

__all__ = ['StructureAtom']


class StructureAtom(Atom):
    """`Atom` class for :class:`~sknano.nanogen.StructureGenerator` classes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    x, y, z : float, optional
        :math:`x, y, z` components of `Atom` position vector relative to
        origin.
    atomID : int, optional
        atom ID
    moleculeID : int, optional
        molecule ID
    atomtype : int, optional
        atom type
    q : {int, float}, optional
        Net charge of `Atom`.
    vx, vy, vz : float, optional
        :math:`v_x, v_y, v_z` components of `Atom` velocity.
    r_units : str, optional
        Units of position components.
    v_units : str, optional
        Units of velocity components.
    nx, ny, nz : int, optional
        :math:`n_x, n_y, n_z` image flags
    CN : int, optional
        `Atom` coordination number.
    NN : sequence, optional
        List of nearest-neighbor `Atom` objects instances

    """

    def __init__(self, element=None, atomID=0, moleculeID=0, atomtype=1, q=0.,
                 mass=None, x=None, y=None, z=None, vx=None, vy=None, vz=None,
                 with_units=False, r_units=None, v_units=None,
                 nx=None, ny=None, nz=None, CN=None, NN=None):

        super(StructureAtom, self).__init__(element=element, m=mass,
                                            x=x, y=y, z=z)

        self._v = Vector(x=vx, y=vy, z=vz)

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
        """Return `Atom` coordination number."""
        return self._CN

    @CN.setter
    def CN(self, value=int):
        """Set `Atom` coordination number."""
        self._CN = int(value)

    @property
    def NN(self):
        """Return list of nearest-neighbor `Atom` objects."""
        return self._NN

    @NN.setter
    def NN(self, value=list):
        """Set list of nearest-neighbor `Atom` objects."""
        self._NN = value

    @property
    def q(self):
        """Charge :math:`q` as multiple of elementary charge :math:`e`.

        """
        return self._q

    @q.setter
    def q(self, value=float):
        """Set `Atom` charge :math:`q`.

        Parameters
        ----------
        value : {int, float}
            net charge on atom as a multiple of the elementary charge
            :math:`e`.

        """
        self._q = value

    @property
    def atomID(self):
        """:attr:`~StructureAtom.atomID`."""
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
        """:attr:`~StructureAtom.moleculeID`."""
        return self._moleculeID

    @moleculeID.setter
    def moleculeID(self, value=int):
        """Set molecule ID of atom.

        Parameters
        ----------
        value : int
            molecule ID

        """
        self._moleculeID = int(value)

    @property
    def atomtype(self):
        """:attr:`~StructureAtom.atomtype`."""
        return self._atomtype

    @atomtype.setter
    def atomtype(self, value=int):
        """Set atom type of atom.

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
        """Set `Atom` :math:`n_z` image flag.

        Parameters
        ----------
        value : float
            :math:`n_z` image flag.

        """
        self._n[2] = int(value)

    @property
    def n(self):
        """:math:`n_x, n_y, n_z` image flags of `Atom`.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`n_x`, :math:`n_y`, :math:`n_z`]
            image flags of `Atom`.

        """
        return self._n

    @n.setter
    def n(self, value=np.ndarray):
        """Set :math:`n_x, n_y, n_z` image flags of `Atom`.

        Parameters
        ----------
        value : array_like
            3-element ndarray of [:math:`n_x`, :math:`n_y`, :math:`n_z`]
            image flags of `Atom`.

        """
        for i, ni in enumerate(value):
            self._n[i] = ni

    @property
    def vx(self):
        """:math:`v_x` component in units of `v_units`.

        Returns
        -------
        float
            :math:`v_x` component in units of `v_units`.

        """
        return self._v.x

    @vx.setter
    def vx(self, value=float):
        """Set `Atom` :math:`v_x` component.

        Parameters
        ----------
        value : float
            :math:`v_x` component in units of `v_units`.

        """
        self._v.x = value

    @property
    def vy(self):
        """:math:`v_y` component in units of `v_units`.

        Returns
        -------
        float
            :math:`v_y` component in units of `v_units`.

        """
        return self._v.y

    @vy.setter
    def vy(self, value=float):
        """Set `Atom` :math:`v_y` component.

        Parameters
        ----------
        value : float
            :math:`v_y` component in units of `v_units`.

        """
        self._v.y = value

    @property
    def vz(self):
        """:math:`v_z` component in units of `v_units`.

        Returns
        -------
        float
            :math:`v_z` component in units of `v_units`.

        """
        return self._v.z

    @vz.setter
    def vz(self, value=float):
        """Set `Atom` :math:`v_z` component.

        Parameters
        ----------
        value : float
            :math:`v_z` component in units of `v_units`.

        """
        self._v.z = value

    @property
    def v(self):
        """:math:`v_x, v_y, v_z` velocity components in default units.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`v_x`, :math:`v_y`, :math:`v_z`]
            velocity components of `Atom`.

        """
        return self._v.components

    @v.setter
    def v(self, value=np.ndarray):
        """Set :math:`x, y, z` components of `Atom` velocity.

        Parameters
        ----------
        value : array_like
            3-element ndarray of [:math:`v_x`, :math:`v_y`, :math:`v_z`]
            velocity components of `Atom`.

        """
        self._v.components = value