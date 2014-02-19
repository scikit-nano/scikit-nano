# -*- coding: utf-8 -*-
"""
==================================================================
Abstract representation of an Atom (:mod:`sknano.chemistry._atom`)
==================================================================

.. currentmodule:: sknano.chemistry._atom

"""
from __future__ import division, absolute_import, print_function
__docformat__ = 'restructuredtext'

import numpy as np

from pkshared.tools.refdata import atomic_masses, atomic_mass_symbol_map, \
    atomic_numbers, atomic_number_symbol_map, element_symbols


class AtomAttributes(object):
    pass


class Atom(object):

    """Class for creating abstract object representing an Atom.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    atomID : int, optional
        atom ID, a LAMMPS atom attribute.
    moleculeID : int, optional
        molecule ID, a LAMMPS atom attribute.
    atomtype : int, optional
        atom type, a LAMMPS atom attribute.
    q : {int, float}, optional
        Net charge of Atom.
    x, y, z : float, optional
        :math:`x, y, z` components of position.
    vx, vy, vz : float, optional
        :math:`v_x, v_y, v_z` components of velocity.
    r_units : str, optional
        Units of position components.
    v_units : str, optional
        Units of velocity components.
    nx, ny, nz : int, optional
        :math:`n_x, n_y, n_z` image flags for LAMMPS atoms.
    CN : int, optional
        Atom coordination number.

    """

    def __init__(self, element=None, atomID=0, moleculeID=0, atomtype=1, q=0.,
                 mass=None, x=None, y=None, z=None, vx=None, vy=None, vz=None,
                 r_units=None, v_units=None, nx=None, ny=None, nz=None,
                 CN=None):

        self._symbol = None
        self._Z = None
        self._m = None

        self._r_units = r_units
        self._r = np.zeros(3, dtype=float)
        for i, ri in enumerate((x, y, z)):
            if ri is not None:
                self._r[i] = ri

        self._v = np.zeros(3, dtype=float)
        self._v_units = v_units
        for i, vi in enumerate((vx, vy, vz)):
            if vi is not None:
                self._v[i] = vi

        self._n = np.zeros(3, dtype=int)
        for i, ni in enumerate((nx, ny, nz)):
            if ni is not None:
                self._n[i] = ni

        if isinstance(element, (int, float)):
            self._Z = int(element)
            idx = self._Z - 1
            try:
                self._symbol = element_symbols[idx]
                self._m = atomic_masses[self._symbol]
            except KeyError as e:
                print(e)
                print('unrecognized element number: {}'.format(element))
        elif isinstance(element, str):
            self._symbol = element
            try:
                self._Z = atomic_numbers[self._symbol]
                self._m = atomic_masses[self._symbol]
            except KeyError as e:
                print(e)
                print('Unrecognized atomic symbol: {}'.format(element))
        else:
            self._symbol = None
            self._Z = None
            if mass is not None and isinstance(mass, (int, float)):
                try:
                    if isinstance(mass, float):
                        self._symbol = atomic_mass_symbol_map[mass]
                    elif isinstance(mass, int):
                        self._symbol = atomic_number_symbol_map[int(mass / 2)]
                    self._Z = atomic_numbers[self._symbol]
                    self._m = atomic_masses[self._symbol]
                except KeyError as e:
                    self._symbol = None
                    self._Z = None
                    self._m = mass
            else:
                self._m = 0

        self._check_type(q, (int, float))
        self._q = q

        self._check_type(atomID, (int, float))
        self._atomID = int(atomID)

        self._check_type(moleculeID, (int, float))
        self._moleculeID = int(moleculeID)

        self._check_type(atomtype, (int, float))
        self._atomtype = int(atomtype)

        self._CN = CN

        self._attributes = ['symbol', 'Z', 'm', 'q', 'r', 'v', 'atomID',
                            'moleculeID', 'atomtype', 'CN']

    def __str__(self):
        """Return string representation of atom."""
        atom_str = ''
        for attr in self._attributes:
            atom_str += \
                'Atom {}: {}\n'.format(attr, getattr(self, '_' + attr))
        return atom_str

    def _check_type(self, value, valid_types):
        if not isinstance(value, valid_types):
            raise TypeError('{} not valid type.\n'.format(value) +
                            '(Valid Types: {})'.format(valid_types))

    @property
    def CN(self):
        """Coordination number."""
        return self._CN

    @CN.setter
    def CN(self, value=int):
        """Set coordination number."""
        self._CN = int(value)

    @property
    def Z(self):
        """Atomic number :math:`Z`.

        Returns
        -------
        int
            Atomic number :math:`Z`.
        """
        return self._Z

    @property
    def symbol(self):
        """Element symbol.

        Returns
        -------
        str
            Element symbol.
        """
        return self._symbol

    @property
    def m(self):
        """Atomic mass :math:`m_a` in atomic mass units.

        Returns
        -------
        float
            Atomic mass :math:`m_a` in atomic mass units.
        """
        return self._m

    @property
    def x(self):
        """:math:`x`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`x`-coordinate in units of **Angstroms**.

        """
        return self._r[0]

    @x.setter
    def x(self, value=float):
        """Set Atom :math:`x`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`x`-coordinate in units of **Angstroms**.

        """
        self._check_type(value, (int, float))
        self._r[0] = float(value)

    @property
    def y(self):
        """:math:`y`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`y`-coordinate in units of **Angstroms**.

        """
        return self._r[1]

    @y.setter
    def y(self, value=float):
        """Set Atom :math:`y`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`y`-coordinate in units of **Angstroms**.

        """
        self._check_type(value, (int, float))
        self._r[1] = float(value)

    @property
    def z(self):
        """:math:`z`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`z`-coordinate in units of **Angstroms**.

        """
        return self._r[2]

    @z.setter
    def z(self, value=float):
        """Set Atom :math:`z`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`z`-coordinate in units of **Angstroms**.

        """
        self._check_type(value, (int, float))
        self._r[2] = float(value)

    @property
    def r(self):
        """:math:`x, y, z` position in units of **Angstroms**.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x, y, z`] coordinates
            of Atom.

        """
        return self._r

    @r.setter
    def r(self, value=np.ndarray):
        """Set :math:`x, y, z` coordinates of atom.

        Parameters
        ----------
        value : array_like
            3-element array of :math:`x, y, z`-coordinates in units of
            **Angstroms**.

        """
        self._check_type(value, np.ndarray)
        for i, ri in enumerate(value):
            self._r[i] = ri

    @property
    def q(self):
        """Charge :math:`q` as multiple of elementary charge :math:`e`.

        """
        return self._q

    @q.setter
    def q(self, value=float):
        """Set Atom charge :math:`q`.

        Parameters
        ----------
        value : {int, float}
            net charge on atom as a multiple of the elementary charge
            :math:`e`.

        """
        self._check_type(value, (int, float))
        self._q = value

    @property
    def atomID(self):
        """Atom ID (*LAMMPS atom attribute*)."""
        return self._atomID

    @atomID.setter
    def atomID(self, value=int):
        """Set atom ID of atom.

        Parameters
        ----------
        value : int
            atom ID

        """
        self._check_type(value, (int, float))
        self._atomID = int(value)

    @property
    def moleculeID(self):
        """Molecule ID (*LAMMPS atom attribute*)."""
        return self._moleculeID

    @moleculeID.setter
    def moleculeID(self, value=int):
        """Set molecule ID of atom.

        Parameters
        ----------
        value : int
            molecule ID

        """
        self._check_type(value, (int, float))
        self._moleculeID = int(value)

    @property
    def atomtype(self):
        """Atom type (*LAMMPS atom attribute*)."""
        return self._atomtype

    @atomtype.setter
    def atomtype(self, value=int):
        """Set atom type of atom.

        Parameters
        ----------
        value : int
            atom type

        """
        self._check_type(value, (int, float))
        self._atomtype = int(value)

    @property
    def vx(self):
        """:math:`v_x` component in units of ``v_units``.

        Returns
        -------
        float
            :math:`v_x` component in units of ``v_units``.

        """
        return self._v[0]

    @vx.setter
    def vx(self, value=float):
        """Set Atom :math:`v_x` component.

        Parameters
        ----------
        value : float
            :math:`v_x` component in units of ``v_units``.

        """
        self._check_type(value, (int, float))
        self._v[0] = float(value)

    @property
    def vy(self):
        """:math:`v_y` component in units of ``v_units``.

        Returns
        -------
        float
            :math:`v_y` component in units of ``v_units``.

        """
        return self._v[1]

    @vy.setter
    def vy(self, value=float):
        """Set Atom :math:`v_y` component.

        Parameters
        ----------
        value : float
            :math:`v_y` component in units of ``v_units``.

        """
        self._check_type(value, (int, float))
        self._v[1] = float(value)

    @property
    def vz(self):
        """:math:`v_z` component in units of ``v_units``.

        Returns
        -------
        float
            :math:`v_z` component in units of ``v_units``.

        """
        return self._v[2]

    @vz.setter
    def vz(self, value=float):
        """Set Atom :math:`v_z` component.

        Parameters
        ----------
        value : float
            :math:`v_z` component in units of ``v_units``.

        """
        self._check_type(value, (int, float))
        self._v[2] = float(value)

    @property
    def v(self):
        """:math:`v_x, v_y, v_z` velocity components in default units.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`v_x`, :math:`v_y`, :math:`v_z`]
            velocity components of Atom.

        """
        return self._v

    @v.setter
    def v(self, value=np.ndarray):
        """Set :math:`x, y, z` components of Atom velocity.

        Parameters
        ----------
        value : array_like
            3-element ndarray of [:math:`v_x`, :math:`v_y`, :math:`v_z`]
            velocity components of Atom.

        """
        self._check_type(value, np.ndarray)
        for i, vi in enumerate(value):
            self._v[i] = vi

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
        self._check_type(value, (int, float))
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
        self._check_type(value, (int, float))
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
        """Set Atom :math:`n_z` image flag.

        Parameters
        ----------
        value : float
            :math:`n_z` image flag.

        """
        self._check_type(value, (int, float))
        self._n[2] = int(value)

    @property
    def n(self):
        """:math:`n_x, n_y, n_z` image flags of Atom.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`n_x`, :math:`n_y`, :math:`n_z`]
            image flags of Atom.

        """
        return self._n

    @n.setter
    def n(self, value=np.ndarray):
        """Set :math:`n_x, n_y, n_z` image flags of Atom.

        Parameters
        ----------
        value : array_like
            3-element ndarray of [:math:`n_x`, :math:`n_y`, :math:`n_z`]
            image flags of Atom.

        """
        self._check_type(value, np.ndarray)
        for i, ni in enumerate(value):
            self._n[i] = ni

    def fix_minus_zero_coords(self, epsilon=1.0e-10):
        """Set really really small negative coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon zero so we don't end up with -0.00000
        coordinates in structure data output.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any component of atom's
            position

        """
        r = self._r.tolist()
        for i, ri in enumerate(r[:]):
            if abs(ri) < epsilon:
                r[i] = 0.0
        self._r[0], self._r[1], self._r[2] = r

    def rezero_coords(self, epsilon=1.0e-10):
        """Set really really small coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        r = self._r.tolist()
        for i, ri in enumerate(r[:]):
            if abs(ri) < epsilon:
                r[i] = 0.0
        self._r[0], self._r[1], self._r[2] = r
