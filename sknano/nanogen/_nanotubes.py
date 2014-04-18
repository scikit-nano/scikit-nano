# -*- coding: utf-8 -*-
"""
=============================================================
Nanotube structure tools (:mod:`sknano.nanogen._nanotubes`)
=============================================================

.. currentmodule:: sknano.nanogen._nanotubes

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from fractions import gcd
from collections import OrderedDict

import numpy as np
np.seterr(all='raise')

try:
    from pint import UnitRegistry
    ureg = UnitRegistry()
    Qty = ureg.Quantity
except ImportError:
    Qty = None

from ..chemistry import Atom
from ..tools import comparison_symbol_operator_mappings
from ..tools.refdata import CCbond, dVDW, grams_per_Da

from ._parameter_luts import param_units, param_symbols, param_strfmt

__all__ = ['Nanotube', 'NanotubeBundle']


class Nanotube(object):
    u"""Class for creating interactive Nanotube objects.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.chemistry.Atoms` 1 and 2
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b
    Lx, Ly, Lz : float, optional
        Spatial extent of nanotube(s) in :math:`x, y, z` dimensions
        in units of **nanometers**. By construction, the nanotube is oriented
        with its principle axis along the :math:`z`-axis. Therefore,
        :math:`L_z` represents the length of the nanotube(s).
        Overrides the :math:`n_x, n_y, n_z` cell values.

        .. versionadded:: 0.2.5

    tube_length : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides `nz` value.

        .. deprecated:: 0.2.5
           Use `Lz` instead

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.

        .. versionadded:: 0.2.6

    units : None, optional
        System of physical units to attach to quantities.
        **This parameter is not yet fully implemented or supported.
        Use at your own risk!**
    with_units : bool, optional
        Attach `units` to physical quantities.
        **This parameter is not yet fully implemented or supported.
        Use at your own risk!**
    verbose : bool, optional
        verbose output

    Examples
    --------

    >>> from sknano.nanogen import Nanotube

    Create a Nanotube with :math:`\\mathbf{C}_{h} = (10, 10)` chirality.

    >>> nt = Nanotube(n=10, m=10, verbose=True)
    n: 10
    m: 10
    t₁: 1
    t₂: -1
    d: 10
    dR: 30
    N: 20
    R: (1, 0)
    θc: 30.00°
    Ch: 42.63 Å
    T: 2.46 Å
    dt: 13.57 Å
    rt: 6.78 Å

    Change the chirality to :math:`\\mathbf{C}_{h} = (20, 10)`.

    >>> nt.n = 20
    n: 20
    m: 10
    t₁: 4
    t₂: -5
    d: 10
    dR: 10
    N: 140
    R: (1, -1)
    θc: 19.11°
    Ch: 65.12 Å
    T: 11.28 Å
    dt: 20.73 Å
    rt: 10.36 Å

    Change the chirality to :math:`\\mathbf{C}_{h} = (20, 0)`.

    >>> nt.m = 0
    n: 20
    m: 0
    t₁: 1
    t₂: -2
    d: 20
    dR: 20
    N: 40
    R: (1, -1)
    θc: 0.00°
    Ch: 49.22 Å
    T: 4.26 Å
    dt: 15.67 Å
    rt: 7.83 Å

    """
    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1, element1='C',
                 element2='C', bond=CCbond, tube_length=None,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False,
                 with_units=False, units=None, verbose=False):

        if tube_length is not None and Lz is None:
            Lz = tube_length

        self._params = OrderedDict()

        # add each parameter in the order I want them to appear in
        # verbose output mode
        self._params['n'] = {}
        self._params['m'] = {}
        self._params['t1'] = {}
        self._params['t2'] = {}
        self._params['d'] = {}
        self._params['dR'] = {}
        self._params['N'] = {}
        #self._params['M'] = {}
        self._params['R'] = {}
        #self._params['bond'] = {}
        self._params['chiral_angle'] = {}
        self._params['Ch'] = {}
        self._params['T'] = {}
        self._params['dt'] = {}
        self._params['rt'] = {}
        self._params['electronic_type'] = {}

        self._n = int(n)
        self._m = int(m)
        self._element1 = element1
        self._element2 = element2

        if bond is None:
            bond = CCbond

        if with_units and Qty is None:
            with_units = False
        self._with_units = with_units

        if with_units and units is None:
            units = 'angstroms'
        self._units = units

        if with_units and isinstance(bond, float):
            self._bond = Qty(bond, units)
        else:
            self._bond = bond

        #self._bond = bond_lengths[element1][element2]

        self._verbose = verbose

        self._L0 = Lz  # store initial value of Lz

        self._fix_Lz = fix_Lz
        self._assume_integer_unit_cells = True
        if fix_Lz:
            self._assume_integer_unit_cells = False

        self._d = None
        self._dR = None
        self._t1 = None
        self._t2 = None
        self._N = None
        self._chiral_angle = None
        self._Ch = None
        self._T = None
        self._dt = None
        self._rt = None
        self._M = None
        self._R = None
        self._p = None
        self._q = None
        self._Lz = None
        self._electronic_type = None
        self._Natoms = None
        self._Natoms_per_tube = None

        self._unit_cell_mass = None
        self._linear_mass_density = None
        self._tube_mass = None

        self._Ntubes = 1
        self._nx = int(nx)
        self._ny = int(ny)

        self._Lx = Lx
        self._Ly = Ly

        if Lz is not None:
            self._T = self.compute_T(n=self._n,
                                     m=self._m,
                                     bond=self._bond,
                                     with_units=with_units,
                                     magnitude=False)
            self._Lz = float(Lz)
            if with_units:
                self._Lz = Qty(self._Lz, 'nanometers')
                self._nz = \
                    self._Lz.to('angstroms').magnitude / self._T.magnitude
            else:
                self._nz = 10 * self._Lz / self._T
        else:
            self._nz = nz

        if self._assume_integer_unit_cells:
            self._nz = int(np.ceil(self._nz))

        for k, v in self.__dict__.iteritems():
            p = k.strip('_')
            if p in self._params.keys():
                self._params[p]['units'] = param_units.get(p)
                self._params[p]['strfmt'] = param_strfmt.get(p)
                if param_symbols.get(p) is not None:
                    self._params[p]['var'] = param_symbols[p]
                else:
                    self._params[p]['var'] = p

        self.compute_tube_params()

    def compute_tube_params(self):
        """Compute/update nanotube parameters."""
        self._d = self.compute_d(n=self._n, m=self._m)
        self._dR = self.compute_dR(n=self._n, m=self._m)
        self._t1 = self.compute_t1(n=self._n, m=self._m)
        self._t2 = self.compute_t2(n=self._n, m=self._m)

        # Compute geometric/physical properties
        self._chiral_angle = \
            self.compute_chiral_angle(n=self._n, m=self._m)

        self._Ch = \
            self.compute_Ch(n=self._n, m=self._m, bond=self._bond,
                            with_units=self._with_units)

        self._T = \
            self.compute_T(n=self._n, m=self._m, bond=self._bond,
                           with_units=self._with_units)

        self._dt = \
            self.compute_dt(n=self._n, m=self._m, bond=self._bond,
                            with_units=self._with_units)

        self._rt = \
            self.compute_rt(n=self._n, m=self._m, bond=self._bond,
                            with_units=self._with_units)

        self._M = \
            self.compute_M(n=self._n, m=self._m)

        self._Lx = \
            self.compute_Lx(n=self._n, m=self._m, nx=self._nx, bond=self._bond,
                            with_units=self._with_units)

        self._Ly = \
            self.compute_Ly(n=self._n, m=self._m, ny=self._ny, bond=self._bond,
                            with_units=self._with_units)

        self._Lz = \
            self.compute_Lz(n=self._n, m=self._m, nz=self._nz, bond=self._bond,
                            with_units=self._with_units)

        self._unit_cell_mass = \
            self.compute_unit_cell_mass(n=self._n, m=self._m,
                                        element1=self._element1,
                                        element2=self._element2,
                                        with_units=self._with_units)

        self._tube_mass = \
            self.compute_tube_mass(n=self._n, m=self._m, nz=self._nz,
                                   element1=self._element1,
                                   element2=self._element2,
                                   with_units=self._with_units)
        self._electronic_type = \
            self.compute_electronic_type(n=self._n, m=self._m)

        # Compute symmetry properties
        self._R = \
            self.compute_R(n=self._n, m=self._m)

        # Compute atomistic properties
        self._N = \
            self.compute_N(n=self._n, m=self._m)

        self._Natoms = \
            self.compute_Natoms(n=self._n, m=self._m)

        self._Natoms_per_tube = \
            self.compute_Natoms_per_tube(n=self._n, m=self._m, nz=self._nz)

        self._linear_mass_density = \
            self.compute_linear_mass_density(n=self._n, m=self._m,
                                             element1=self._element1,
                                             element2=self._element2,
                                             with_units=self._with_units)

        for k, v in self.__dict__.iteritems():
            p = k.strip('_')
            if p in self._params.keys():
                self._params[p]['val'] = v

        if self._verbose:
            for p, pdict in self._params.iteritems():
                pvar = pdict['var']
                try:
                    pval = pdict['val'].magnitude
                except AttributeError:
                    pval = pdict['val']
                punits = pdict['units']
                pstrfmt = pdict['strfmt']
                if punits is not None:
                    if pstrfmt is not None:
                        print(u"{}: {}{}".format(pvar,
                                                 pstrfmt.format(pval),
                                                 punits))
                    else:
                        print(u"{}: {}{}".format(pvar, pval, punits))
                else:
                    print(u"{}: {}".format(pvar, pval))
            print()

    @property
    def bond(self):
        u"""Bond length in **\u212b**."""
        return self._bond

    @bond.setter
    def bond(self, value):
        u"""Set bond length in **\u212b**."""
        self._bond = float(value)
        self.compute_tube_params()

    @property
    def element1(self):
        """Element symbol of :class:`~sknano.chemistry.Atom` 1."""
        return self._element1

    @element1.setter
    def element1(self, value):
        """Set element symbol of :class:`~sknano.chemistry.Atom` 1."""
        self._element1 = value
        self.compute_tube_params()

    @property
    def element2(self):
        """Element symbol of :class:`~sknano.chemistry.Atom` 2."""
        return self._element2

    @element2.setter
    def element2(self, value):
        """Set element symbol of :class:`~sknano.chemistry.Atom` 2."""
        self._element2 = value
        self.compute_tube_params()

    @property
    def n(self):
        """Chiral index :math:`n`.

        The component of the chiral vector :math:`\\mathbf{C}_h`
        along :math:`\\mathbf{a}_1`:

        .. math::

           \\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)

        """
        return self._n

    @n.setter
    def n(self, value):
        """Set chiral index :math:`n`"""
        self._n = int(value)
        self.compute_tube_params()

    @property
    def m(self):
        """Chiral index :math:`m`.

        The component of the chiral vector :math:`\\mathbf{C}_h`
        along :math:`\\mathbf{a}_2`:

        .. math::

           \\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)

        """
        return self._m

    @m.setter
    def m(self, value):
        """Set chiral index :math:`m`"""
        self._m = int(value)
        self.compute_tube_params()

    @property
    def d(self):
        """:math:`d=\\gcd{(n, m)}`

        :math:`d` is the **G**\ reatest **C**\ ommon **D**\ ivisor of
        :math:`n` and :math:`m`.

        """
        return self._d

    @classmethod
    def compute_d(cls, n=int, m=int):
        """Compute :math:`d=\\gcd{(n, m)}`

        :math:`d` is the **G**\ reatest **C**\ ommon **D**\ ivisor of
        :math:`n` and :math:`m`.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

        Returns
        -------
        int
            Greatest Common Divisor of :math:`n` and :math:`m`

        """
        return gcd(n, m)

    @property
    def dR(self):
        """:math:`d_R=\\gcd{(2n + m, 2m + n)}`

        :math:`d_R` is the **G**\ reatest **C**\ ommon **D**\ ivisor of
        :math:`2n + m` and :math:`2m + n`.

        """
        return self._dR

    @classmethod
    def compute_dR(cls, n=int, m=int):
        """Compute :math:`d_R=\\gcd{(2n + m, 2m + n)}`

        :math:`d_R` is the **G**\ reatest **C**\ ommon **D**\ ivisor of
        :math:`2n + m` and :math:`2m + n`.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

        Returns
        -------
        int
            greatest common divisor of :math:`2n+m` and :math:`2m+n`

        """
        return gcd(2 * m + n, 2 * n + m)

    @property
    def N(self):
        """Number of graphene hexagons in nanotube *unit cell*.

        .. math::

           N = \\frac{4(n^2 + m^2 + nm)}{d_R}

        """
        return self._N

    @classmethod
    def compute_N(cls, n=int, m=int):
        """Compute :math:`N = \\frac{2(n^2+m^2+nm)}{d_R}`.

        :math:`N` is the number of graphene hexagons mapped to a nanotube
        *unit cell*.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

        Returns
        -------
        int
            Number of hexagons per nanotube unit cell:
            :math:`N = \\frac{2(n^2+m^2+nm)}{d_R}`.

        """
        dR = Nanotube.compute_dR(n=n, m=m)
        try:
            return int(2 * (n**2 + m**2 + n * m) / dR)
        except ZeroDivisionError:
            return 0

    @property
    def t1(self):
        """:math:`t_{1} = \\frac{2m + n}{d_{R}}`

        where :math:`d_R = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_1`:

        .. math::

           \\mathbf{T} = t_1\\mathbf{a}_{1} + t_2\\mathbf{a}_2

        """
        return self._t1

    @classmethod
    def compute_t1(cls, n=int, m=int):
        """Compute :math:`t_1 = \\frac{2m + n}{d_R}`

        where :math:`d_R = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_1`:

        .. math::

           \\mathbf{T} = t_1\\mathbf{a}_{1} + t_2\\mathbf{a}_2

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

        Returns
        -------
        int
            :math:`t_1`

        """
        dR = Nanotube.compute_dR(n=n, m=m)
        try:
            return int((2 * m + n) / dR)
        except ZeroDivisionError:
            return 0

    @property
    def t2(self):
        """:math:`t_2 = -\\frac{2n + m}{d_R}`

        where :math:`d_R = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_2`:

        .. math::

           \\mathbf{T} = t_1\\mathbf{a}_1 + t_2\\mathbf{a}_2

        """
        return self._t2

    @classmethod
    def compute_t2(cls, n=int, m=int):
        """Compute :math:`t_2 = -\\frac{2n + m}{d_R}`

        where :math:`d_R = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_2`:

        .. math::

           \\mathbf{T} = t_1\\mathbf{a}_1 + t_2\\mathbf{a}_2

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

        Returns
        -------
        int
            :math:`t_2`

        """
        dR = Nanotube.compute_dR(n=n, m=m)
        try:
            return -int((2 * n + m) / dR)
        except ZeroDivisionError:
            return 0

    @property
    def Ch(self):
        u"""Nanotube circumference :math:`|\\mathbf{C}_h|` in **\u212b**"""
        return self._Ch

    @classmethod
    def compute_Ch(cls, n=int, m=int, bond=None, with_units=False,
                   units='angstrom', magnitude=True):
        u"""Compute nanotube circumference :math:`|\\mathbf{C}_{h}|` in
        **\u212b**.

        .. math::

           |\\mathbf{C}_h| = a\\sqrt{n^2 + m^2 + nm} =
           \\sqrt{3}a_{\\mathrm{CC}}\\sqrt{n^2 + m^2 + nm}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        bond : float, optional
            Distance between nearest neighbor atoms (i.e., bond length).
            Must be in units of **\u212b**. Default value is
            the carbon-carbon bond length in graphite, approximately
            :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

        Returns
        -------
        float
            Nanotube circumference :math:`|\\mathbf{C}_h|` in \u212b.

        """
        if bond is None:
            bond = CCbond

        if with_units and isinstance(bond, float) and Qty is not None:
            bond = Qty(bond, units)

        if magnitude and with_units:
            try:
                return bond.magnitude * np.sqrt(3 * (n**2 + m**2 + n * m))
            except AttributeError:
                return bond * np.sqrt(3 * (n**2 + m**2 + n * m))
        else:
            return bond * np.sqrt(3 * (n**2 + m**2 + n * m))

    @property
    def chiral_angle(self):
        """Chiral angle :math:`\\theta_c` in **degrees**.

        .. math::

           \\theta_c = \\tan^{-1}\\left(\\frac{\\sqrt{3} m}{2n + m}\\right)

        """
        return self._chiral_angle

    @classmethod
    def compute_chiral_angle(cls, n=int, m=int, rad2deg=True):
        """Compute chiral angle :math:`\\theta_c`.

        .. math::

           \\theta_c = \\tan^{-1}\\left(\\frac{\\sqrt{3} m}{2n + m}\\right)

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        rad2deg : bool, optional
            If `True`, return angle in degrees.

        Returns
        -------
        float
            Chiral angle :math:`\\theta_{c}` in
            degrees (default) or radians (if `rad2deg=False`).

        """
        theta = np.arctan(np.sqrt(3) * m / (2 * n + m))
        #return np.arccos((2*n + m) / (2 * np.sqrt(n**2 + m**2 + n*m)))
        if rad2deg:
            return np.degrees(theta)
        else:
            return theta

    @property
    def T(self):
        u"""Length of nanotube unit cell :math:`|\\mathbf{T}|` in \u212b.

        .. math::

           |\\mathbf{T}| = \\frac{\\sqrt{3} |\\mathbf{C}_{h}|}{d_{R}}

        """
        return self._T

    @classmethod
    def compute_T(cls, n=None, m=None, bond=None, with_units=False,
                  units='angstrom', length=True, magnitude=True):
        u"""Compute length of nanotube unit cell :math:`|\\mathbf{T}|` in
        \u212b.

        .. math::

           |\\mathbf{T}| = \\frac{\\sqrt{3} |\\mathbf{C}_{h}|}{d_{R}}
           = \\frac{\\sqrt{3}a\\sqrt{n^2 + m^2 + nm}}{d_{R}}
           = \\frac{3a_{\\mathrm{CC}}\\sqrt{n^2 + m^2 + nm}}{d_{R}}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        bond : float, optional
            Distance between nearest neighbor atoms (i.e., bond length).
            Must be in units of **\u212b**. Default value is
            the carbon-carbon bond length in graphite, approximately
            :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b
        length : bool, optional
            Compute the magnitude (i.e., length) of the translation vector.
        magnitude : bool, optional
            Return the length of the translation vector **without** units.
        with_units : bool, optional
        units : str, optional

        Returns
        -------
        float or Qty or 2-tuple of ints
            If `length` is `True` and `with_units` is `False`, then
            return the length of unit cell in \u212b.

            If `length` is `True` and `with_units` is `True` and
            magnitude is `True`, then return the length of unit cell
            as a float, but with units in `units`.

            If `length` is `True` and `with_units` is `True` and
            magnitude is `False`, then return the length of unit cell
            as a `Quantity` instance and with units in `units`.

            If `length` is `False`, return the componets of the
            translation vector as a 2-tuple of ints
            (:math:`t_1`, :math:`t_2`).

        """

        if length:
            if bond is None:
                bond = CCbond

            if with_units and isinstance(bond, float) and Qty is not None:
                bond = Qty(bond, units)

            Ch = Nanotube.compute_Ch(n=n, m=m, bond=bond,
                                     with_units=with_units, units=units,
                                     magnitude=magnitude)
            dR = Nanotube.compute_dR(n=n, m=m)

            try:
                return np.sqrt(3) * Ch / dR
            except (FloatingPointError, ZeroDivisionError):
                return 0
        else:
            t1 = Nanotube.compute_t1(n=n, m=m)
            t2 = Nanotube.compute_t2(n=n, m=m)

            return (t1, t2)

    @property
    def dt(self):
        u"""Nanotube diameter :math:`d_t = \\frac{|\\mathbf{C}_h|}{\\pi}`
        in \u212b."""
        return self._dt

    @classmethod
    def compute_dt(cls, n=int, m=int, bond=None, with_units=False,
                   units='angstrom', magnitude=True):
        u"""Compute nanotube diameter :math:`d_t` in \u212b.

        .. math::

           d_t = \\frac{|\\mathbf{C}_h|}{\\pi}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        bond : float, optional
            Distance between nearest neighbor atoms (i.e., bond length).
            Must be in units of **\u212b**. Default value is
            the carbon-carbon bond length in graphite, approximately
            :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

        Returns
        -------
        float
            Nanotube diameter :math:`d_t` in \u212b.

        """
        Ch = Nanotube.compute_Ch(n, m, bond=bond, with_units=with_units,
                                 units=units, magnitude=magnitude)
        return Ch / np.pi

    @property
    def rt(self):
        u"""Nanotube radius :math:`r_t = \\frac{|\\mathbf{C}_h|}{2\\pi}`
        in \u212b."""
        return self._rt

    @classmethod
    def compute_rt(cls, n=int, m=int, bond=None, with_units=False,
                   units='angstrom', magnitude=True):
        u"""Compute nanotube radius :math:`r_t` in \u212b.

        .. math::

           r_t = \\frac{|\\mathbf{C}_h|}{2\\pi}


        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        bond : float, optional
            Distance between nearest neighbor atoms (i.e., bond length).
            Must be in units of **\u212b**. Default value is
            the carbon-carbon bond length in graphite, approximately
            :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

        Returns
        -------
        float
            Nanotube radius :math:`r_t` in \u212b.

        """
        Ch = Nanotube.compute_Ch(n=n, m=m, bond=bond, with_units=with_units,
                                 units=units, magnitude=magnitude)
        return Ch / (2 * np.pi)

    @property
    def M(self):
        """:math:`M = np - nq`

        :math:`M` is the number of multiples of the translation vector
        :math:`\\mathbf{T}` in the vector :math:`N\\mathbf{R}`.

        """
        return self._M

    @classmethod
    def compute_M(cls, n=int, m=int):
        """Compute :math:`M = mp - nq`

        :math:`M` is the number of multiples of the translation vector
        :math:`\\mathbf{T}` in the vector :math:`N\\mathbf{R}`.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

        Returns
        -------
        int
            :math:`M = mp - nq`

        """
        p, q = Nanotube.compute_R(n=n, m=m)
        return m * p - n * q

    @property
    def R(self):
        """Symmetry vector :math:`\\mathbf{R} = (p, q)`.

        .. math::

           \\mathbf{R} = p\\mathbf{a}_1 + q\\mathbf{a}_2

        """
        return self._R

    @classmethod
    def compute_R(cls, n=int, m=int, bond=None, with_units=False,
                  units='angstrom', length=False, magnitude=True):
        u"""Compute symmetry vector :math:`\\mathbf{R} = (p, q)`.

        The *symmetry vector* is any lattice vector of the unfolded graphene
        layer that represents a *symmetry operation* of the nanotube. The
        symmetry vector :math:`\\mathbf{R}` can be written as:

        .. math::

           \\mathbf{R} = p\\mathbf{a}_1 + q\\mathbf{a}_2

        where :math:`p` and :math:`q` are integers.
        The *symmetry vector* represents a *symmetry operation* of the nanotube
        which arises as a *screw translation*, which is a combination of
        a rotation :math:`\\psi` and translation :math:`\\tau`. The symmetry
        operation of the nanotube can be written as:

        .. math::

           R = (\\psi|\\tau)

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        bond : float, optional
            Distance between nearest neighbor atoms (i.e., bond length).
            Must be in units of **\u212b**. Default value is
            the carbon-carbon bond length in graphite, approximately
            :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b
        length : bool, optional
            If `True`, return :math:`|\\mathbf{R}|`.
        magnitude : bool, optional
            If `True`, return :math:`|\\mathbf{R}|` in units of **\u212b**.

        Returns
        -------
        (p, q) : tuple
            2-tuple of ints -- components of :math:`\\mathbf{R}`.
        float
            Length of :math:`\\mathbf{R}` (:math:`|\\mathbf{R}|`) if `length`
            is `True` in units of **\u212b**.

        """
        N = Nanotube.compute_N(n=n, m=m)
        t1 = Nanotube.compute_t1(n=n, m=m)
        t2 = Nanotube.compute_t2(n=n, m=m)

        p = q = 0
        for i in xrange(0, t1 + n + 1):
            for j in xrange(t2, m + 1):
                R = t1 * j - t2 * i
                if R == 1:
                    M = m * i - n * j
                    if M > 0 and M <= N:
                        p = i
                        q = j

        if length:
            if bond is None:
                bond = CCbond

            if with_units and isinstance(bond, float) and Qty is not None:
                bond = Qty(bond, units)

            if magnitude and with_units:
                try:
                    return bond.magnitude * np.sqrt(3 * (p**2 + q**2 + p * q))
                except AttributeError:
                    return bond * np.sqrt(3 * (p**2 + q**2 + p * q))
            else:
                return bond * np.sqrt(3 * (p**2 + q**2 + p * q))
        else:
            return (p, q)

    @classmethod
    def compute_R_chiral_angle(cls, n=int, m=int, rad2deg=True):
        """Compute "chiral angle" of symmetry vector :math:`\\theta_R`.

        .. math::

           \\theta_R = \\tan^{-1}\\left(\\frac{\\sqrt{3}q}{2p + q}\\right)

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        rad2deg : bool, optional
            If `True`, return angle in degrees

        Returns
        -------
        float
            Chiral angle of *symmetry vector* :math:`\\theta_R` in
            degrees (default) or radians (if `rad2deg=False`).


        """
        p, q = Nanotube.compute_R(n=n, m=m)
        theta = np.arctan((np.sqrt(3) * q) / (2 * p + q))
        if rad2deg:
            return np.degrees(theta)
        else:
            return theta

    @classmethod
    def compute_so(cls, n=int, m=int, bond=None):
        u"""Compute symmetry operation :math:`(\\psi|\\tau)`.

        The *symmetry vector* `R` represents a *symmetry
        operation* of the nanotube which arises as a *screw translation*, which
        is a combination of a rotation :math:`\\psi` and translation
        :math:`\\tau`.
        The symmetry operation of the nanotube can be written as:

        .. math::

           R = (\\psi|\\tau)

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        bond : float, optional
            Distance between nearest neighbor atoms (i.e., bond length).
            Must be in units of **\u212b**. Default value is
            the carbon-carbon bond length in graphite, approximately
            :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

        Returns
        -------
        (psi, tau) : tuple
            2-tuple of floats -- :math:`\\psi` in radians and
            :math:`\\tau` in \u212b.

        """
        psi = Nanotube.compute_psi(n=n, m=m)
        tau = Nanotube.compute_tau(n=n, m=m, bond=bond)
        return (psi, tau)

    @classmethod
    def compute_psi(cls, n=int, m=int):
        """Compute rotation component of symmetry operation
        :math:`\\psi` in **radians**.

        .. math::

           \\psi = \\frac{2\\pi}{N}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

        Returns
        -------
        float
            Rotation component of symmetry operation :math:`\\psi`
            in **radians**.

        """
        N = Nanotube.compute_N(n=n, m=m)
        try:
            return 2 * np.pi / N
        except (FloatingPointError, ZeroDivisionError):
            return 0

    @classmethod
    def compute_tau(cls, n=int, m=int, bond=None, with_units=False,
                    units='angstrom', magnitude=True):
        u"""Compute translation component of symmetry operation
        :math:`\\tau` in **\u212b**.

        .. math::

           \\tau = \\frac{M|\\mathbf{T}|}{N}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        bond : float, optional
            Distance between nearest neighbor atoms (i.e., bond length).
            Must be in units of **\u212b**. Default value is
            the carbon-carbon bond length in graphite, approximately
            :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

        Returns
        -------
        float
            Translation component of symmetry operation :math:`\\tau`
            in **\u212b**.

        """
        M = Nanotube.compute_M(n=n, m=m)
        N = Nanotube.compute_N(n=n, m=m)
        T = Nanotube.compute_T(n=n, m=m, bond=bond, with_units=with_units,
                               units=units, length=True, magnitude=magnitude)
        try:
            return M * T / N
        except ZeroDivisionError:
            return 0

    @property
    def nx(self):
        """Number of nanotubes along the :math:`x`-axis."""
        return int(self._nx)

    @nx.setter
    def nx(self, value=int):
        """Set :math:`n_x`"""
        self._nx = value

    @property
    def ny(self):
        """Number of nanotubes along the :math:`y`-axis."""
        return int(self._ny)

    @ny.setter
    def ny(self, value=int):
        """Set :math:`n_y`"""
        self._ny = value

    @property
    def nz(self):
        """Number of nanotube unit cells along the :math:`z`-axis."""
        return self._nz

    @nz.setter
    def nz(self, value):
        """Set number of nanotube unit cells along the :math:`z`-axis."""
        if self._assume_integer_unit_cells:
            self._nz = int(value)
        else:
            self._nz = value
        self.compute_tube_params()

    @property
    def Lx(self):
        """Spatial extent of nanotubes along :math:`x`-axis in **nm**."""
        return self._Lx

    @classmethod
    def compute_Lx(cls, n=int, m=int, nx=int, bond=None, with_units=False,
                   units='nanometer', magnitude=True):
        u"""Compute :math:`L_x` in **nanometers**.

        .. math::

           L_x = n_x (d_t + d_{\\mathrm{vdW}})

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        nx : int
            Number of nanotubes along :math:`x`-axis.
        bond : float, optional
            Distance between nearest neighbor atoms (i.e., bond length).
            Must be in units of **\u212b**. Default value is
            the carbon-carbon bond length in graphite, approximately
            :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

        Returns
        -------
        float
            :math:`L_x` in **nanometers**

        """
        dt = Nanotube.compute_dt(n=n, m=m, bond=bond, with_units=with_units,
                                 units=units, magnitude=magnitude)
        Lx = nx * (dt + dVDW)
        if with_units:
            try:
                Lx.ito(units)
            except Exception:
                Lx.ito('nanometer')
        else:
            Lx = Lx / 10

        if magnitude and with_units:
            try:
                return Lx.magnitude
            except AttributeError:
                return Lx
        else:
            return Lx

    @property
    def Ly(self):
        """Spatial extent of nanotubes along :math:`y`-axis in **nm**."""
        return self._Ly

    @classmethod
    def compute_Ly(cls, n=int, m=int, ny=int, bond=None, with_units=False,
                   units='nanometer', magnitude=True):
        """Compute :math:`L_y` in **nanometers**.

        .. math::

           L_y = n_y (d_t + d_{\\mathrm{vdW}})

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        ny : int
            Number of nanotubes along :math:`y`-axis.
        bond : float, optional
            Distance between nearest neighbor atoms (i.e., bond length).
            Must be in units of **\u212b**. Default value is
            the carbon-carbon bond length in graphite, approximately
            :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

        Returns
        -------
        float
            :math:`L_y` in **nanometers**

        """
        dt = Nanotube.compute_dt(n=n, m=m, bond=bond, with_units=with_units,
                                 units=units, magnitude=magnitude)
        Ly = ny * (dt + dVDW)
        if with_units:
            try:
                Ly.ito(units)
            except Exception:
                Ly.ito('nanometer')
        else:
            Ly = Ly / 10

        if magnitude and with_units:
            try:
                return Ly.magnitude
            except AttributeError:
                return Ly
        else:
            return Ly

    @property
    def Lz(self):
        """Nanotube length :math:`L_z = L_{\\mathrm{tube}}` in
        **nanometers**."""
        return self._Lz

    @Lz.setter
    def Lz(self, value=float):
        """Set tube length"""
        self._Lz = value

    @classmethod
    def compute_Lz(cls, n=int, m=int, nz=float, bond=None, with_units=False,
                   units='nanometer', magnitude=True):
        """Compute :math:`L_z = L_{\\mathrm{tube}}` in **nanometers**.

        .. math::

           L_z = n_z |\\mathbf{T}|


        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        nz : {int, float}
            Number of nanotube unit cells
        bond : float, optional
            Distance between nearest neighbor atoms (i.e., bond length).
            Must be in units of **\u212b**. Default value is
            the carbon-carbon bond length in graphite, approximately
            :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

        Returns
        -------
        float
            :math:`L_z = L_{\\mathrm{tube}}` in **nanometers**

        """
        T = Nanotube.compute_T(n=n, m=m, bond=bond, with_units=with_units,
                               units=units, length=True, magnitude=False)
        Lz = nz * T
        if with_units:
            try:
                Lz.ito(units)
            except Exception:
                Lz.ito('nanometer')
        else:
            Lz = Lz / 10

        if magnitude and with_units:
            try:
                return Lz.magnitude
            except AttributeError:
                return Lz
        else:
            return Lz

    @property
    def electronic_type(self):
        """Nanotube electronic type.

        .. versionadded:: 0.2.7

        The electronic type is determined as follows:

        if :math:`(2n + m)\\,\\mathrm{mod}\\,3=0`, the nanotube is
        **metallic**.

        if :math:`(2n + m)\\,\\mathrm{mod}\\,3=1`, the nanotube is
        **semiconducting, type 1**.

        if :math:`(2n + m)\\,\\mathrm{mod}\\,3=2`, the nanotube is
        **semiconducting, type 2**.

        The :math:`x\\,\\mathrm{mod}\\,y` notation is mathematical
        shorthand for the *modulo* operation, which computes the
        **remainder** of the division of :math:`x` by :math:`y`.
        So, for example, all *armchair* nanotubes must be metallic
        since the chiral indices satisfy: :math:`2n + m = 2n + n = 3n` and
        therefore :math:`3n\\,\\mathrm{mod}\\,3` i.e. the remainder of the
        division of :math:`3n/3=n` is always zero.

        .. note::
           Mathematically, :math:`(2n + m)\\,\\mathrm{mod}\\,3` is equivalent
           to :math:`(n - m)\\,\\mathrm{mod}\\,3` when distinguishing
           between metallic and semiconducting. However, when
           distinguishing between semiconducting types,
           one must be careful to observe the following convention:

           * Semiconducting, **type 1** means:

             * :math:`(2n + m)\\,\\mathrm{mod}\\,3=1`
             * :math:`(n - m)\\,\\mathrm{mod}\\,3=2`

           * Semiconducting, **type 2** means:

             * :math:`(2n + m)\\,\\mathrm{mod}\\,3=2`
             * :math:`(n - m)\\,\\mathrm{mod}\\,3=1`

        """
        return self._electronic_type

    @classmethod
    def compute_electronic_type(cls, n=int, m=int):
        """Compute nanotube electronic type.

        .. versionadded:: 0.2.7

        The electronic type is determined as follows:

        if :math:`(2n + m)\\,\\mathrm{mod}\\,3=0`, the nanotube is
        **metallic**.

        if :math:`(2n + m)\\,\\mathrm{mod}\\,3=1`, the nanotube is
        **semiconducting, type 1**.

        if :math:`(2n + m)\\,\\mathrm{mod}\\,3=2`, the nanotube is
        **semiconducting, type 2**.

        The :math:`x\\,\\mathrm{mod}\\,y` notation is mathematical
        shorthand for the *modulo* operation, which computes the
        **remainder** of the division of :math:`x` by :math:`y`.
        So, for example, all *armchair* nanotubes must be metallic
        since the chiral indices satisfy: :math:`2n + m = 2n + n = 3n` and
        therefore :math:`3n\\,\\mathrm{mod}\\,3` i.e. the remainder of the
        division of :math:`3n/3=n` is always zero.

        .. note::
           Mathematically, :math:`(2n + m)\\,\\mathrm{mod}\\,3` is equivalent
           to :math:`(n - m)\\,\\mathrm{mod}\\,3` when distinguishing
           between metallic and semiconducting. However, when
           distinguishing between semiconducting types,
           one must be careful to observe the following convention:

           * Semiconducting, **type 1** means:

             * :math:`(2n + m)\\,\\mathrm{mod}\\,3=1`
             * :math:`(n - m)\\,\\mathrm{mod}\\,3=2`

           * Semiconducting, **type 2** means:

             * :math:`(2n + m)\\,\\mathrm{mod}\\,3=2`
             * :math:`(n - m)\\,\\mathrm{mod}\\,3=1`

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

        Returns
        -------
        str

        """
        if (2 * n + m) % 3 == 1:
            return 'semiconducting, type 1'
        elif (2 * n + m) % 3 == 2:
            return 'semiconducting, type 2'
        else:
            return 'metallic'

    @property
    def Natoms(self):
        """Number of atoms in nanotube *unit cell*.

        .. math::

           N_{\\mathrm{atoms}} = 2N = \\frac{4(n^2 + m^2 + nm)}{d_R}

        where :math:`N` is the number of graphene hexagons mapped to the
        nanotube unit cell.

        """
        return self._Natoms

    @classmethod
    def compute_Natoms(cls, n=int, m=int):
        """Compute :math:`N_{\mathrm{atoms/cell}} = 2N`.

        .. math::

           N_{\\mathrm{atoms}} = 2N = \\frac{4(n^2 + m^2 + nm)}{d_R}

        where :math:`N` is the number of graphene hexagons mapped to the
        nanotube unit cell.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

        Returns
        -------
        `2N` : int
            Number of atoms in nanotube *unit cell*:
            N_{\\mathrm{atoms}} = 2N = \\frac{4(n^2 + m^2 + nm)}{d_R}


        """
        N = Nanotube.compute_N(n=n, m=m)
        return 2 * N

    @property
    def Natoms_per_tube(self):
        """Number of atoms in nanotube :math:`N_{\\mathrm{atoms/tube}}`."""
        return self._Natoms_per_tube

    @classmethod
    def compute_Natoms_per_tube(cls, n=int, m=int, nz=float):
        """Compute :math:`N_{\\mathrm{atoms/tube}}`

        .. math::

           N_{\\mathrm{atoms/tube}} = N_{\\mathrm{atoms/cell}} \\times
           n_{z-\\mathrm{cells}}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        nz : {int, float}
            Number of nanotube unit cells

        Returns
        -------
        int
            :math:`N_{\\mathrm{atoms/tube}}`
        """
        Natoms = Nanotube.compute_Natoms(n=n, m=m)
        return int(Natoms * nz)

    @property
    def unit_cell_mass(self):
        """Unit cell mass in atomic mass units."""
        return self._unit_cell_mass

    @classmethod
    def compute_unit_cell_mass(cls, n=int, m=int, element1=None, element2=None,
                               with_units=False, units='Da', magnitude=True):
        """Compute nanotube unit cell mass in atomic mass units.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        element1, element2 : {str, int}, optional
            Element symbol or atomic number of basis
            :class:`~sknano.chemistry.Atoms` 1 and 2
        with_units : bool, optional
            Attach `units` to physical quantities.
            **This parameter is not yet fully implemented or supported.
            Use at your own risk!**
        units : None, optional
            System of physical units to attach to quantities.
            **This parameter is not yet fully implemented or supported.
            Use at your own risk!**

        Returns
        -------
        float
            Unit cell mass in **Daltons**.

        Notes
        -----

        .. todo::

           Handle different elements and perform accurate calculation by
           determining number of atoms of each element.

        """
        N = Nanotube.compute_N(n=n, m=m)

        if element1 is None:
            element1 = 'C'
        if element2 is None:
            element2 = 'C'

        mass = N * (Atom(element1).m + Atom(element2).m)
        if with_units and Qty is not None:
            mass = Qty(mass, 'Da')
            if units is not None and units not in ('Da', 'amu'):
                mass.ito(units)
            if magnitude:
                try:
                    return mass.magnitude
                except AttributeError:
                    return mass
            else:
                return mass
        else:
            return mass

    @property
    def linear_mass_density(self):
        """Linear mass density of nanotube in g/nm."""
        return self._linear_mass_density

    @classmethod
    def compute_linear_mass_density(cls, n=int, m=int, bond=None,
                                    element1=None, element2=None,
                                    with_units=False, units='grams/nm',
                                    magnitude=True):
        """Compute nanotube linear mass density (mass per unit length) in
        **grams/nm**.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        element1, element2 : {str, int}, optional
            Element symbol or atomic number of basis
            :class:`~sknano.chemistry.Atoms` 1 and 2
        bond : float, optional
            Distance between nearest neighbor atoms (i.e., bond length).
            Must be in units of **\u212b**. Default value is
            the carbon-carbon bond length in graphite, approximately
            :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b
        with_units : bool, optional
            Attach `units` to physical quantities.
            **This parameter is not yet fully implemented or supported.
            Use at your own risk!**
        units : None, optional
            System of physical units to attach to quantities.
            **This parameter is not yet fully implemented or supported.
            Use at your own risk!**

        Returns
        -------
        float
            Linear mass density in units of **g/nm**.

        """
        mass = Nanotube.compute_unit_cell_mass(n=n, m=m,
                                               element1=element1,
                                               element2=element2,
                                               with_units=False)
        T = Nanotube.compute_T(n=n, m=m, bond=bond,
                               with_units=False, length=True)

        try:
            linear_mass_density = mass / T

            if with_units and Qty is not None:
                linear_mass_density = Qty(linear_mass_density, 'Da/angstrom')
                if units is not None:
                    linear_mass_density.ito(units)

                if magnitude:
                    try:
                        return linear_mass_density.magnitude
                    except AttributeError:
                        return linear_mass_density
                else:
                    return linear_mass_density
            else:
                # there are 1.6605e-24 grams / Da and 10 angstroms / nm
                linear_mass_density *= 10 * grams_per_Da
                return linear_mass_density
        except ZeroDivisionError:
            return 0

    @property
    def Ntubes(self):
        """Number of nanotubes."""
        return int(self._Ntubes)

    @classmethod
    def compute_symmetry_chiral_angle(cls, n=int, m=int, rad2deg=True):
        """Alias for :meth:`Nanotube.compute_R_chiral_angle`."""
        return Nanotube.compute_R_chiral_angle(n=n, m=m, rad2deg=rad2deg)

    @classmethod
    def compute_tube_diameter(cls, n=int, m=int, bond=None, with_units=False,
                              units='angstrom', magnitude=True):
        """Alias for :meth:`Nanotube.compute_dt`"""
        return Nanotube.compute_dt(n=n, m=m, bond=bond, with_units=with_units,
                                   units=units, magnitude=magnitude)

    @classmethod
    def compute_tube_radius(cls, n=int, m=int, bond=None, with_units=False,
                            units='angstrom', magnitude=True):
        """Alias for :meth:`Nanotube.compute_rt`"""
        return Nanotube.compute_rt(n, m, bond=bond, with_units=with_units,
                                   units=units, magnitude=magnitude)

    @property
    def tube_length(self):
        """Alias for :attr:`Nanotube.Lz`"""
        return self.Lz

    @classmethod
    def compute_tube_length(cls, n=int, m=int, nz=float, bond=None,
                            with_units=False, units='nanometer',
                            magnitude=True):
        """Alias for :meth:`Nanotube.compute_Lz`"""
        return Nanotube.compute_Lz(n=n, m=m, nz=nz, bond=bond,
                                   with_units=with_units, units=units,
                                   magnitude=magnitude)

    @property
    def tube_mass(self):
        """Nanotube mass in **grams**."""
        return self._tube_mass

    @classmethod
    def compute_tube_mass(cls, n=int, m=int, nz=float, element1=None,
                          element2=None, with_units=False, units='grams',
                          magnitude=True):
        """Compute nanotube mass in **grams**.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        nz : {int, float}
            Number of nanotube unit cells
        element1, element2 : {str, int}, optional
            Element symbol or atomic number of basis
            :class:`~sknano.chemistry.Atoms` 1 and 2
        with_units : bool, optional
            Attach `units` to physical quantities.
            **This parameter is not yet fully implemented or supported.
            Use at your own risk!**
        units : None, optional
            System of physical units to attach to quantities.
            **This parameter is not yet fully implemented or supported.
            Use at your own risk!**

        Returns
        -------
        float
            Nanotube mass in **grams**.

        Notes
        -----

        .. todo::

           Handle different elements and perform accurate calculation by
           determining number of atoms of each element.

        """
        Natoms_per_tube = \
            Nanotube.compute_Natoms_per_tube(n=n, m=m, nz=nz)

        if element1 is None:
            element1 = 'C'
        if element2 is None:
            element2 = 'C'

        atom1 = Atom(element1)
        atom2 = Atom(element2)

        mass = Natoms_per_tube * (atom1.m + atom2.m) / 2
        if with_units and Qty is not None:
            mass = Qty(mass, 'Da')
            if units is not None and units not in ('Da', 'amu'):
                mass.ito(units)
            if magnitude:
                try:
                    return mass.magnitude
                except AttributeError:
                    return mass
            else:
                return mass
        else:
            # there are 1.6605e-24 grams / Da
            mass *= grams_per_Da
            return mass

    @classmethod
    def filter_Ch_list(cls, Ch_list=None, property_filters=None):
        """Filter list of chiralities by list of properties.

        Parameters
        ----------
        Ch_list : sequence
            List of chiralities
        property_filters : sequence

        """
        if Ch_list is not None and property_filters is not None:
            filtered_list = Ch_list[:]
            try:
                for filter_index, (prop, cmp_symbol, value) in \
                        enumerate(property_filters, start=1):
                    cmp_op = comparison_symbol_operator_mappings[cmp_symbol]
                    tmp_list = []
                    for Ch in filtered_list:
                        n, m = Ch
                        nanotube = Nanotube(n=n, m=m)
                        try:
                            if cmp_op(getattr(nanotube, prop), value):
                                tmp_list.append(Ch)
                        except AttributeError:
                            break
                    filtered_list = tmp_list[:]
            except ValueError as e:
                print(e)
            finally:
                return filtered_list
        else:
            return None


class NanotubeBundle(Nanotube):
    u"""Class for creating interactive Nanotube bundle objects.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.chemistry.Atoms` 1 and 2
    bond : float, optional
        Bond length between nearest neighbor atoms.
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b
    Lx, Ly, Lz : float, optional
        Length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.

        .. versionadded:: 0.2.5

    tube_length : float, optional
        Length of nanotube in units of **nanometers**.

        .. deprecated:: 0.2.5
           Use `Lz` instead

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.

        .. versionadded:: 0.2.6

    verbose : bool, optional
        Verbose output.

    """
    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond, tube_length=None,
                 vdw_spacing=dVDW, bundle_packing=None, bundle_geometry=None,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False, with_units=False,
                 units=None, verbose=False):

        if tube_length is not None and Lz is None:
            Lz = tube_length

        super(NanotubeBundle, self).__init__(
            n=n, m=m, nx=nx, ny=ny, nz=nz,
            element1=element1, element2=element2,
            bond=bond, Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz,
            with_units=with_units, units=units, verbose=verbose)

        self._vdw_spacing = vdw_spacing
        self._bundle_packing = bundle_packing
        self._bundle_geometry = bundle_geometry

        self._Natoms_per_bundle = None

        self._bundle_mass = None
        self._bundle_density = None

        self.compute_bundle_params()

    def compute_bundle_params(self, d_vdw=None):
        """Compute/update nanotube bundle parameters."""
        super(NanotubeBundle, self).compute_tube_params()
        self._Ntubes = self.compute_Ntubes(nx=self._nx, ny=self._ny)
        self._bundle_mass = \
            self.compute_bundle_mass(n=self._n, m=self._m,
                                     nx=self._nx, ny=self._ny, nz=self._nz,
                                     with_units=False)
        self._Natoms_per_bundle = self._Ntubes * self._Natoms_per_tube
        if d_vdw is None:
            d_vdw = self._vdw_spacing

        if self._with_units and isinstance(d_vdw, float):
            d_vdw = Qty(d_vdw, 'angstrom')

        self._bundle_density = \
            self.compute_bundle_density(n=self.n, m=self.m, d_vdw=d_vdw,
                                        bond=self.bond, element1=self.element1,
                                        element2=self.element2,
                                        with_units=self._with_units)

    @property
    def Natoms_per_bundle(self):
        """Number of atoms in nanotube bundle.

        .. versionadded:: 0.2.5

        """
        return self._Natoms_per_bundle

    @Natoms_per_bundle.setter
    def Natoms_per_bundle(self, value=int):
        """Set Natoms per bundle.

        .. versionadded:: 0.2.5

        """
        self._Natoms_per_bundle = value

    @classmethod
    def compute_Natoms_per_bundle(cls, n=int, m=int, nz=float, Ntubes=None,
                                  nx=None, ny=None):
        """Compute number of atoms in nanotube bundle.

        .. versionadded:: 0.2.5

        .. math::

           N_{\\mathrm{atoms/bundle}} = N_{\\mathrm{tubes}}\\times
           N_{\\mathrm{atoms/tube}}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        nz : {int, float}
        Ntubes : {None, int}, optional
        nx, ny : {None, int}, optional

        Returns
        -------
        Natoms_per_bundle : int
            Number of atoms in a nanotube bundle.

        Raises
        ------
        TypeError
            If `Ntubes` is `None` and `nx` or `ny` are `None`

        """
        Natoms_per_tube = \
            Nanotube.compute_Natoms_per_tube(n=n, m=m, nz=nz)
        if Ntubes is None and (nx is None or ny is None):
            raise TypeError('`Ntubes` or `nx` and `ny` must be specified as '
                            'integers')
        elif Ntubes is None and nx is not None and ny is not None:
            Ntubes = nx * ny
        Natoms_per_bundle = Ntubes * Natoms_per_tube

        return Natoms_per_bundle

    @property
    def Ntubes(self):
        """Number of nanotubes in nanotube bundle: :math:`n_x\\times n_y`."""
        return int(self._Ntubes)

    @Ntubes.setter
    def Ntubes(self, value=int):
        """Set Ntubes."""
        self._Ntubes = value
        self._nx = value
        self._ny = 1
        self.compute_bundle_params()

    @classmethod
    def compute_Ntubes(cls, nx=int, ny=int):
        """Compute number of nanotubes in nanotube bundle.

        Parameters
        ----------
        nx, ny : int, optional
            Number of repeat unit cells along :math:`x, y`-axes.

        Returns
        -------
        int
            Number of nanotubes in nanotube bundle: :math:`n_x\\times n_y`

        """
        return int(nx * ny)

    @property
    def bundle_density(self):
        """Nanotube bundle mass density :math:`\\rho_{\\mathrm{bundle}}` in
        :math:`\\mathrm{g/cm^3}`."""
        return self._bundle_density

    @classmethod
    def compute_bundle_density(cls, n=int, m=int, d_vdw=None, bond=None,
                               element1=None, element2=None,
                               with_units=False, units='g/cm**3',
                               magnitude=True):
        u"""Compute nanotube bundle mass density
        :math:`\\rho_{\\mathrm{bundle}}(n, m)` in :math:`\\mathrm{g/cm^3}`.

        .. math::

           \\rho_{\\mathrm{bundle}}(n, m) = \\frac{8\\pi^2 m_{\\mathrm{C}}
           \\sqrt{n^2 + m^2 + nm}}{9\\sqrt{3}a_{\\mathrm{CC}}^3 \\times
           \\left(\\sqrt{n^2 + m^2 + nm} +
           \\frac{\\pi d_{\\mathrm{vdW}}}{\\sqrt{3}a_{\\mathrm{CC}}}\\right)^2}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        d_vdw : int
            van der Waals distance between nearest-neighbor tubes
        bond : float, optional
            Bond length.

        Returns
        -------
        float
            :math:`\\rho_{\\mathrm{bundle}}` in units of
            :math:`\\mathrm{\\frac{g}{cm^3}}`

        """
        if bond is None:
            bond = CCbond
        else:
            try:
                bond = bond.magnitude
            except AttributeError:
                pass

        if d_vdw is None:
            if n == m:
                d_vdw = 3.38
            elif (m == 0) or (n == 0):
                d_vdw = 3.41
            else:
                d_vdw = 3.39
        else:
            try:
                d_vdw = d_vdw.magnitude
            except AttributeError:
                pass

        if element1 is None:
            element1 = 'C'
        if element2 is None:
            element2 = 'C'

        if element1 == element2:
            bundle_density = 8 * np.pi**2 * Atom(element1).m * \
                np.sqrt(n**2 + m**2 + n*m) / \
                (9 * np.sqrt(3) * bond**3 *
                    (np.sqrt(n**2 + m**2 + n*m) +
                        np.pi * d_vdw / (np.sqrt(3) * bond))**2)
        else:
            bundle_density = 0

        if with_units and Qty is not None:
            bundle_density = Qty(bundle_density, 'Da/angstroms**3')
            if units is not None:
                bundle_density.ito(units)

            if magnitude:
                try:
                    return bundle_density.magnitude
                except AttributeError:
                    return bundle_density
            else:
                return bundle_density
        else:
            # there are 1.6605e-24 grams / Da and 1e-8 cm / angstrom
            bundle_density *= grams_per_Da / (1e-8)**3
            return bundle_density

    @property
    def bundle_mass(self):
        """Nanotube bundle mass :math:`M_{\\mathrm{bundle}}` in **grams**."""
        return self._bundle_mass

    @classmethod
    def compute_bundle_mass(cls, n=int, m=int, nx=int, ny=int, nz=float,
                            element1=None, element2=None, with_units=False,
                            units='grams', magnitude=True):
        """Compute nanotube bundle mass :math:`M_{\\mathrm{bundle}}`
        in **grams**.

        .. math::

           M_{\\mathrm{bundle}} = N_{\\mathrm{tubes}} M_{\\mathrm{tube}}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
        nx, ny, nz : int, optional
            Number of repeat unit cells in the :math:`x, y, z` dimensions
        element1, element2 : {str, int}, optional
            Element symbol or atomic number of basis
            :class:`~sknano.chemistry.Atoms` 1 and 2

        Returns
        -------
        float
            Nanotube bundle mass :math:`M_{\\mathrm{bundle}}` in **grams**.

        """
        Ntubes = NanotubeBundle.compute_Ntubes(nx=nx, ny=ny)
        tube_mass = Nanotube.compute_tube_mass(n=n, m=m, nz=nz,
                                               element1=element1,
                                               element2=element2,
                                               with_units=with_units,
                                               units=units,
                                               magnitude=magnitude)
        return Ntubes * tube_mass
