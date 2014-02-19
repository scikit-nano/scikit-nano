# -*- coding: utf-8 -*-
"""
=============================================================
Nanotube structure tools (:mod:`sknano.nanogen._nanotube`)
=============================================================

.. currentmodule:: sknano.nanogen._nanotube

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

#import itertools
import warnings
warnings.filterwarnings('ignore')  # to suppress the Pint UnicodeWarning

from fractions import gcd
from collections import OrderedDict
from pint import UnitRegistry
ureg = UnitRegistry()
Qty = ureg.Quantity

import numpy as np

from pkshared.tools.refdata import CCbond

from ..chemistry import Atom, Atoms
param_units = {}
param_units['dt'] = \
    param_units['rt'] = \
    param_units['Ch'] = \
    param_units['T'] = \
    param_units['bond'] = u' \u212B'
param_units['chiral_angle'] = u'\u00b0'

param_symbols = {}
param_symbols['dt'] = u'd\u209C'
param_symbols['rt'] = u'r\u209C'
param_symbols['Ch'] = u'C\u2095'
param_symbols['t1'] = u't\u2081'
param_symbols['t2'] = u't\u2082'
param_symbols['chiral_angle'] = u'\u03b8\u1d04'

param_strfmt = {}
param_strfmt['Ch'] = \
    param_strfmt['T'] = \
    param_strfmt['dt'] = \
    param_strfmt['rt'] = \
    param_strfmt['chiral_angle'] = '{:.2f}'
param_strfmt['bond'] = '{:.3f}'

__all__ = ['param_units', 'param_symbols', 'param_strfmt',
           'NanotubeError', 'ChiralityError',
           'Nanotube', 'NanotubeBundle']


class NanotubeError(Exception):
    """Base class for nanotube module exceptions."""
    pass


class ChiralityError(NanotubeError):
    """Exception raised for errors in the chirality indices (n, m)."""
    pass


class TubeLengthError(NanotubeError):
    """Exception raised for errors in the length."""
    pass


class CellError(NanotubeError):
    """Exception raised for errors in the number of unit cells."""
    pass


class Nanotube(object):
    u"""Class for creating interactive Nanotube objects.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        bond length between nearest neighbor atoms.
        Must be in units of **Angstroms**. Default value is
        the carbon-carbon bond length in graphite:
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.421` \u212b ([SoFaCNTs]_)
    Lx, Ly, Lz : float, optional
        Length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.

        .. versionadded:: 0.2.5

    tube_length : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the ``nz`` value.

        .. deprecated:: 0.2.5
           Use ``Lz`` instead

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If ``True``, then
        non integer :math:`n_z` cells are permitted.

        .. versionadded:: 0.2.6

    verbose : bool, optional
        verbose output

    References
    ----------
    .. [SoFaCNTs] Science of Fullerenes and Carbon Nanotubes,
       M. Dresselhaus, G. Dresselhaus, P. Eklund, 1996, p. 760.

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
    M: 10
    R: (1, 0)
    bond: 1.421 Å
    Cₕ: 42.63 Å
    T: 2.46 Å
    dₜ: 13.57 Å
    rₜ: 6.78 Å
    θᴄ: 30.00°

    Change the chirality to :math:`\\mathbf{C}_{h} = (20, 10)`.

    >>> nt.n = 20
    n: 20
    m: 10
    t₁: 4
    t₂: -5
    d: 10
    dR: 10
    N: 140
    M: 30
    R: (1, -1)
    bond: 1.421 Å
    Cₕ: 65.12 Å
    T: 11.28 Å
    dₜ: 20.73 Å
    rₜ: 10.36 Å
    θᴄ: 19.11°

    Change the chirality to :math:`\\mathbf{C}_{h} = (20, 0)`.

    >>> nt.m = 0
    n: 20
    m: 0
    t₁: 1
    t₂: -2
    d: 20
    dR: 20
    N: 40
    M: 20
    R: (1, -1)
    bond: 1.421 Å
    Cₕ: 49.22 Å
    T: 4.26 Å
    dₜ: 15.67 Å
    rₜ: 7.83 Å
    θᴄ: 0.00°

    """
    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1, element1='C',
                 element2='C', bond=CCbond, tube_length=None,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False,
                 with_units=False, verbose=False):

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
        self._params['M'] = {}
        self._params['R'] = {}
        self._params['bond'] = {}
        self._params['Ch'] = {}
        self._params['T'] = {}
        self._params['dt'] = {}
        self._params['rt'] = {}
        self._params['chiral_angle'] = {}
        self._params['electronic_type'] = {}

        self._n = int(n)
        self._m = int(m)
        self._element1 = element1
        self._element2 = element2

        if bond is None:
            bond = CCbond

        if with_units and isinstance(bond, float):
            self._bond = Qty(bond, 'angstroms')
        else:
            self._bond = bond

        #self._bond = bond_lengths[element1][element2]

        self._with_units = with_units
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
        self._Ch = None
        self._T = None
        self._dt = None
        self._rt = None
        self._chiral_angle = None
        self._N = None
        self._M = None
        self._R = None
        self._p = None
        self._q = None
        self._Natoms = None
        self._Lz = None
        self._Natoms_per_tube = None
        self._electronic_type = None

        try:
            self._a = np.sqrt(3) * self._bond.magnitude
        except AttributeError:
            self._a = np.sqrt(3) * self._bond

        self._a1 = np.zeros(2, dtype=float)
        self._a2 = np.zeros(2, dtype=float)

        self._a1[0] = self._a2[0] = np.sqrt(3) / 2 * self._a
        self._a1[1] = 1 / 2 * self._a
        self._a2[1] = -self._a1[1]

        self._b1 = np.zeros(2, dtype=float)
        self._b2 = np.zeros(2, dtype=float)

        self._b1[0] = self._b2[0] = \
            1 / np.sqrt(3) * 2 * np.pi / self._a
        self._b1[1] = 2 * np.pi / self._a
        self._b2[1] = -self._b1[1]

        if with_units:
            self._a1 = Qty(self._a1, 'angstrom')
            self._a2 = Qty(self._a2, 'angstrom')

            self._b1 = Qty(self._b1, '1/angstrom')
            self._b2 = Qty(self._b2, '1/angstrom')

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
        """Compute nanotube parameters."""
        self._d = self.compute_d(n=self._n, m=self._m)
        self._dR = self.compute_dR(n=self._n, m=self._m)
        self._t1 = self.compute_t1(n=self._n, m=self._m)
        self._t2 = self.compute_t2(n=self._n, m=self._m)

        # Compute geometric properties
        self._Ch = self.compute_Ch(n=self._n, m=self._m,
                                   with_units=self._with_units,
                                   bond=self._bond)
        self._T = self.compute_T(n=self._n, m=self._m,
                                 with_units=self._with_units,
                                 bond=self._bond)
        self._dt = self.compute_dt(n=self._n, m=self._m,
                                   with_units=self._with_units,
                                   bond=self._bond)
        self._rt = self.compute_rt(n=self._n, m=self._m,
                                   with_units=self._with_units,
                                   bond=self._bond)
        self._chiral_angle = self.compute_chiral_angle(n=self._n, m=self._m)
        self._M = self.compute_M(n=self._n, m=self._m)

        # Compute physical properties
        self._Lz = \
            self.compute_Lz(n=self._n, m=self._m, nz=self._nz,
                            with_units=self._with_units,
                            bond=self._bond)
        self._tube_mass = \
            self.compute_tube_mass(n=self._n, m=self._m, nz=self._nz,
                                   with_units=self._with_units)
        self._electronic_type = \
            self.compute_electronic_type(n=self._n, m=self._m)

        # Compute symmetry properties
        self._R = self.compute_R(n=self._n, m=self._m)

        # Compute atomistic properties
        self._N = self.compute_N(n=self._n, m=self._m)
        self._Natoms = self.compute_Natoms(n=self._n, m=self._m)
        self._Natoms_per_tube = \
            self.compute_Natoms_per_tube(n=self._n,
                                         m=self._m,
                                         nz=self._nz)

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
    def n(self):
        """Chiral index :math:`n`"""
        return self._n

    @n.setter
    def n(self, value):
        """Set chiral index :math:`n`"""
        self._n = int(value)
        self.compute_tube_params()

    @property
    def m(self):
        """Chiral index :math:`m`"""
        return self._m

    @m.setter
    def m(self, value):
        """Set chiral index :math:`m`"""
        self._m = int(value)
        self.compute_tube_params()

    @property
    def bond(self):
        """Bond length in **Angstroms**."""
        return self._bond

    @property
    def a(self):
        """Length of graphene unit cell vector."""
        return self._a

    @property
    def a1(self):
        """:math:`a_1` unit vector."""
        return self._a1

    @property
    def a2(self):
        """:math:`a_2` unit vector."""
        return self._a2

    @property
    def b1(self):
        """:math:`b_1` reciprocal lattice vector."""
        return self._b1

    @property
    def b2(self):
        """:math:`b_2` reciprocal lattice vector."""
        return self._b2

    @property
    def t1(self):
        """:math:`t_{1} = \\frac{2m + n}{d_{R}}`

        where :math:`d_{R} = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_{1}`:

        .. math::

           \\mathbf{T} = t_{1}\\mathbf{a}_{1} + t_{2}\\mathbf{a}_{2}

        """
        return self._t1

    @classmethod
    def compute_t1(cls, n=int, m=int):
        """Compute :math:`t_{1} = \\frac{2m + n}{d_{R}}`

        where :math:`d_{R} = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_{1}`:

        .. math::

           \\mathbf{T} = t_{1}\\mathbf{a}_{1} + t_{2}\\mathbf{a}_{2}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.

        Returns
        -------
        int
            :math:`t_{1}`

        """
        dR = Nanotube.compute_dR(n=n, m=m)
        return int((2 * m + n) / dR)

    @property
    def t2(self):
        """:math:`t_{2} = -\\frac{2n + m}{d_{R}}`

        where :math:`d_{R} = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_{2}`:

        .. math::

           \\mathbf{T} = t_{1}\\mathbf{a}_{1} + t_{2}\\mathbf{a}_{2}

        """
        return self._t2

    @classmethod
    def compute_t2(cls, n=int, m=int):
        """Compute :math:`t_{2} = -\\frac{2n + m}{d_{R}}`

        where :math:`d_{R} = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_{2}`:

        .. math::

           \\mathbf{T} = t_{1}\\mathbf{a}_{1} + t_{2}\\mathbf{a}_{2}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.

        Returns
        -------
        int
            :math:`t_{2}`

        """
        dR = Nanotube.compute_dR(n=n, m=m)
        return -int((2 * n + m) / dR)

    @property
    def N(self):
        """Number of hexagons per nanotube unit cell :math:`N`:

        .. math::

           N = \\frac{4(n^2 + m^2 + nm)}{d_{R}}

        """
        return self._N

    @classmethod
    def compute_N(cls, n=int, m=int):
        """Compute :math:`N = \\frac{2(n^2+m^2+nm)}{d_{R}}`.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.

        Returns
        -------
        int
            Number of hexagons per nanotube unit cell:
            :math:`N = \\frac{2(n^2+m^2+nm)}{d_{R}}`.

        """
        dR = Nanotube.compute_dR(n=n, m=m)
        return int(2 * (n**2 + m**2 + n * m) / dR)

    @property
    def Natoms(self):
        """Number of atoms per nanotube **unit cell** :math:`2N`.

        .. math::

           N_{\\mathrm{atoms}} = 2N = \\frac{4(n^2 + m^2 + nm)}{d_{R}}

        """
        return self._Natoms

    @classmethod
    def compute_Natoms(cls, n=int, m=int):
        """Compute :math:`N_{\mathrm{atoms/cell}} = 2N`.

        .. math::

           N_{\\mathrm{atoms}} = 2N = \\frac{4(n^2 + m^2 + nm)}{d_{R}}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.

        Returns
        -------
        int
            Number of atoms per nanotube *unit cell*:
            N_{\\mathrm{atoms}} = 2N = \\frac{4(n^2 + m^2 + nm)}{d_{R}}


        """
        N = Nanotube.compute_N(n=n, m=m)
        return 2 * N

    @property
    def R(self):
        """Symmetry vector :math:`\\mathbf{R} = (p, q)`.

        .. math::

           \\mathbf{R} = p\\mathbf{a}_{1} + q\\mathbf{a}_{2}

        """
        return self._R

    @classmethod
    def compute_R(cls, n=int, m=int, bond=None, units='angstrom',
                  with_units=False, length=False, magnitude=True):
        """Compute symmetry vector :math:`\\mathbf{R} = (p, q)`

        .. math::

           \\mathbf{R} = p\\mathbf{a}_{1} + q\\mathbf{a}_{2}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.
        length : bool, optional
            if `True`, return the length of R
        magnitude : bool, optional
            if `True`, return the length of R without units

        Returns
        -------
        (p, q) : tuple
            2-tuple of ints which are the components of R vector
        float
            length of R if `length` is `True`

        """
        t1 = Nanotube.compute_t1(n=n, m=m)
        t2 = Nanotube.compute_t2(n=n, m=m)
        N = Nanotube.compute_N(n=n, m=m)

        p = None
        q = None

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

            if with_units and isinstance(bond, float):
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

    @property
    def M(self):
        """The number of :math:`\\mathbf{T}` in :math:`N\\mathbf{R}`"""
        return self._M

    @classmethod
    def compute_M(cls, n=int, m=int):
        """Compute :math:`M = mp - nq`

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.

        Returns
        -------
        int
            :math:`M = mp - nq`

        """
        p, q = Nanotube.compute_R(n=n, m=m)
        return m * p - n * q

    @property
    def Ch(self):
        """Nanotube circumference :math:`\\mathbf{C}_{h} = (n, m)`"""
        return self._Ch

    @classmethod
    def compute_Ch(cls, n=int, m=int, bond=None, with_units=False,
                   units='angstrom', magnitude=True):
        """Compute the nanotube circumference.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.
        bond : float, optional
            distance between nearest neighbor atoms.
            Must be in units of **Angstroms**.

        Returns
        -------
        float
            nanotube circumference in Angstroms

        """
        if bond is None:
            bond = CCbond

        if with_units and isinstance(bond, float):
            bond = Qty(bond, units)

        if magnitude and with_units:
            try:
                return bond.magnitude * np.sqrt(3 * (n**2 + m**2 + n * m))
            except AttributeError:
                return bond * np.sqrt(3 * (n**2 + m**2 + n * m))
        else:
            return bond * np.sqrt(3 * (n**2 + m**2 + n * m))

    @property
    def dt(self):
        """Nanotube diameter :math:`d_{t} = \\frac{|\\mathbf{C}_{h}|}{\\pi}`"""
        return self._dt

    @classmethod
    def compute_tube_diameter(cls, n=int, m=int, bond=None, with_units=False,
                              units='angstrom', magnitude=True):
        """Alias for :meth:`Nanotube.compute_dt`"""
        return Nanotube.compute_dt(n=n, m=m, bond=bond, units=units,
                                   with_units=with_units,
                                   magnitude=magnitude)

    @classmethod
    def compute_dt(cls, n=int, m=int, bond=None, with_units=False,
                   units='angstrom', magnitude=True):
        """Compute nanotube diameter :math:`d_{t}`

        .. math::

           d_{t} = \\frac{|\\mathbf{C}_{h}|}{\\pi}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} =
            (n, m)`.
        bond : float, optional
            distance between nearest neighbor atoms.
            Must be in units of **Angstroms**.

        Returns
        -------
        float
            nanotube diameter in Angstroms

        """
        Ch = Nanotube.compute_Ch(n, m, bond=bond, with_units=with_units,
                                 units=units, magnitude=magnitude)
        return Ch / np.pi

    @property
    def rt(self):
        """Nanotube radius :math:`r_{t} = \\frac{|\\mathbf{C}_{h}|}{2\\pi}`"""
        return self._rt

    @classmethod
    def compute_tube_radius(cls, n=int, m=int, bond=None, with_units=False,
                            units='angstrom', magnitude=True):
        """Alias for :meth:`Nanotube.compute_rt`"""
        return Nanotube.compute_rt(n, m, bond=bond, with_units=with_units,
                                   units=units, magnitude=magnitude)

    @classmethod
    def compute_rt(cls, n=int, m=int, bond=None, with_units=False,
                   units='angstrom', magnitude=True):
        """Compute nanotube radius :math:`r_{t}`

        .. math::

           r_{t} = \\frac{|\\mathbf{C}_{h}|}{2\\pi}


        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.
        bond : float, optional
            distance between nearest neighbor atoms.
            Must be in units of **Angstroms**.

        Returns
        -------
        float
            nanotube radius in Angstroms

        """
        Ch = Nanotube.compute_Ch(n=n, m=m, bond=bond, with_units=with_units,
                                 units=units, magnitude=magnitude)
        return Ch / (2 * np.pi)

    @property
    def d(self):
        """:math:`d=\\gcd{(n, m)}`"""
        return self._d

    @classmethod
    def compute_d(cls, n=int, m=int):
        """Compute :math:`d=\\gcd{(n, m)}`

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} +
            m\\mathbf{a}_{2} = (n, m)`.

        Returns
        -------
        int
            greatest common divisor of :math:`n` and :math:`m`

        """
        return gcd(n, m)

    @property
    def dR(self):
        """:math:`d_R=\\gcd{(2n + m, 2m + n)}`"""
        return self._dR

    @classmethod
    def compute_dR(cls, n=int, m=int):
        """Compute :math:`d_R=\\gcd{(2n + m, 2m + n)}`

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} +
            m\\mathbf{a}_{2} = (n, m)`.

        Returns
        -------
        int
            greatest common divisor of :math:`2n+m` and :math:`2m+n`

        """
        return gcd(2 * m + n, 2 * n + m)

    @property
    def chiral_angle(self):
        """Chiral angle :math:`\\theta_{c}`.

        .. math::

           \\theta_{c} = \\tan^{-1}\\left({\\frac{\\sqrt{3} m}{2n + m}}\\right)

        """
        return self._chiral_angle

    @classmethod
    def compute_chiral_angle(cls, n=int, m=int, rad2deg=True):
        """Compute chiral angle :math:`\\theta_{c}`

        .. math::

           \\theta_{c} = \\tan^{-1}\\left({\\frac{\\sqrt{3} m}{2n + m}}\\right)

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} +
            m\\mathbf{a}_{2} = (n, m)`.
        rad2deg : bool, optional
            If `True`, return angle in degrees

        Returns
        -------
        float
            chiral angle :math:`\\theta_{c}` in degrees.

        """
        theta = np.arctan(np.sqrt(3) * m / (2 * n + m))
        #return np.arccos((2*n + m) / (2 * np.sqrt(n**2 + m**2 + n*m)))
        if rad2deg:
            return np.degrees(theta)
        else:
            return theta

    @property
    def T(self):
        """Unit cell length :math:`|\\mathbf{T}|`.

        .. math::

           |\\mathbf{T}| = \\frac{\\sqrt{3} |\\mathbf{C}_{h}|}{d_{R}}

        """
        return self._T

    @classmethod
    def compute_T(cls, n=None, m=None, bond=None, with_units=False,
                  units='angstrom', length=True, magnitude=True):
        """Compute unit cell length :math:`|\\mathbf{T}|`

        .. math::

           |\\mathbf{T}| = \\frac{\\sqrt{3} |\\mathbf{C}_{h}|}{d_{R}}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} +
            m\\mathbf{a}_{2} = (n, m)`.
        bond : float, optional
            distance between nearest neighbor atoms.
            Must be in units of **Angstroms**.

        Returns
        -------
        float
            length of unit cell in Angstroms

        """

        if length:
            if bond is None:
                bond = CCbond

            if with_units and isinstance(bond, float):
                bond = Qty(bond, units)

            Ch = Nanotube.compute_Ch(n=n, m=m, bond=bond,
                                     with_units=with_units, units=units,
                                     magnitude=magnitude)
            dR = Nanotube.compute_dR(n=n, m=m)

            return np.sqrt(3) * Ch / dR
        else:
            t1 = Nanotube.compute_t1(n=n, m=m)
            t2 = Nanotube.compute_t2(n=n, m=m)

            return (t1, t2)

    @property
    def Natoms_per_tube(self):
        """Number of atoms per nanotube :math:`N_{\\mathrm{atoms/tube}}`."""
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
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.
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
    def nx(self):
        """Number of nanotube unit cells along the :math:`x`-axis."""
        return int(self._nx)

    @nx.setter
    def nx(self, value=int):
        """Set :math:`n_x`"""
        self._nx = value

    @property
    def ny(self):
        """Number of nanotube unit cells along the :math:`y`-axis."""
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
    def nz(self, value=float):
        """Set number of unit cells along the :math:`z`-axis."""
        if self._assume_integer_unit_cells:
            self._nz = int(value)
        else:
            self._nz = value
        self.compute_tube_params()

    @property
    def Ntubes(self):
        """Number of nanotubes."""
        return int(self._Ntubes)

    @property
    def Lx(self):
        """Nanotube :math:`L_{\\mathrm{x}}` in **nanometers**."""
        return self._Lx

    @property
    def Ly(self):
        """Nanotube :math:`L_{\\mathrm{y}}` in **nanometers**."""
        return self._Ly

    @property
    def Lz(self):
        """Nanotube length :math:`L_{\\mathrm{tube}}` in **nanometers**."""
        return self._Lz

    @Lz.setter
    def Lz(self, value=float):
        """Set tube length"""
        self._Lz = value

    @classmethod
    def compute_Lz(cls, n=int, m=int, nz=float, bond=None, with_units=False,
                   units='angstrom', magnitude=True):
        """Compute :math:`L_{\\mathrm{tube}}`.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.
        nz : {int, float}
            Number of nanotube unit cells
        bond : float, optional
            bond length

        Returns
        -------
        float
            :math:`L_{\\mathrm{tube}}` in **nanometers**

        """
        T = Nanotube.compute_T(n=n, m=m, bond=bond, with_units=with_units,
                               units=units, length=True, magnitude=False)
        Lz = nz * T
        if with_units:
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
    def tube_length(self):
        """Nanotube length :math:`L_{\\mathrm{tube}}` in **nanometers**."""
        return self.Lz

    @classmethod
    def compute_tube_length(cls, n=int, m=int, nz=float, bond=None,
                            with_units=False, units='angstrom',
                            magnitude=True):
        """Compute :math:`L_{\\mathrm{tube}}`.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.
        nz : {int, float}
            Number of nanotube unit cells
        bond : float, optional
            bond length

        Returns
        -------
        float
            :math:`L_{\\mathrm{tube}}` in **nanometers**

        """
        return Nanotube.compute_Lz(n=n, m=m, nz=nz, bond=bond,
                                   with_units=with_units, units=units,
                                   magnitude=magnitude)

    @property
    def tube_mass(self):
        """Nanotube mass in grams."""
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
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.
        nz : {int, float}
            Number of nanotube unit cells

        Returns
        -------
        float

        """
        Natoms_per_tube = \
            Nanotube.compute_Natoms_per_tube(n=n, m=m, nz=nz)

        if element1 is None:
            element1 = 'C'
        if element2 is None:
            element2 = 'C'

        atom1 = Atom(element1)
        atom2 = Atom(element2)

        mass = Natoms_per_tube * Atoms([atom1, atom2]).m
        if with_units:
            mass = Qty(mass, units)

        if magnitude and with_units:
            try:
                return mass.magnitude
            except AttributeError:
                return mass
        else:
            return mass

    @property
    def electronic_type(self):
        """Nanotube electronic type.

        .. versionadded:: 0.2.7

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
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.

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

    @classmethod
    def compute_symmetry_chiral_angle(cls, n=int, m=int, rad2deg=True):
        """Compute "chiral angle" of symmetry vector :math:`\\mathbf{R}`

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} +
            m\\mathbf{a}_{2} = (n, m)`.
        rad2deg : bool, optional
            If `True`, return angle in degrees

        Returns
        -------
        float
            chiral angle :math:`\\theta_{R}`

        """
        p, q = Nanotube.compute_R(n=n, m=m)
        theta = np.arctan((np.sqrt(3) * q) / (2 * p + q))
        if rad2deg:
            return np.degrees(theta)
        else:
            return theta

    @classmethod
    def compute_tau(cls, n=int, m=int, bond=None, with_units=False,
                    units='angstrom', magnitude=True):
        M = Nanotube.compute_M(n=n, m=m)
        N = Nanotube.compute_N(n=n, m=m)
        T = Nanotube.compute_T(n=n, m=m, bond=bond, with_units=with_units,
                               units=units, length=True, magnitude=magnitude)
        return M * T / N

    @classmethod
    def compute_psi(cls, n=int, m=int):
        N = Nanotube.compute_N(n=n, m=m)
        return 2 * np.pi / N


class NanotubeBundle(Nanotube):
    u"""Class for creating interactive Nanotube **Bundles**.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        Bond length between nearest neighbor atoms.
        Must be in units of **Angstroms**. Default value is
        the carbon-carbon bond length in graphite:
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.421` \u212b ([SoFaCNTs]_)
    Lx, Ly, Lz : float, optional
        Length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.

        .. versionadded:: 0.2.5

    tube_length : float, optional
        Length of nanotube in units of **nanometers**.

        .. deprecated:: 0.2.5
           Use ``Lz`` instead

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If ``True``, then
        non integer :math:`n_z` cells are permitted.

        .. versionadded:: 0.2.6

    verbose : bool, optional
        Verbose output.

    """
    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond, tube_length=None,
                 vdw_spacing=3.4, bundle_packing=None, bundle_geometry=None,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False, with_units=False,
                 verbose=False):

        if tube_length is not None and Lz is None:
            Lz = tube_length

        super(NanotubeBundle, self).__init__(
            n=n, m=m, nx=nx, ny=ny, nz=nz,
            element1=element1, element2=element2,
            bond=bond, Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz,
            with_units=with_units, verbose=verbose)

        self._vdw_spacing = vdw_spacing
        self._bundle_packing = bundle_packing
        self._bundle_geometry = bundle_geometry

        self._Natoms_per_bundle = None

        self._bundle_mass = None
        self._bundle_density = None

        self.compute_bundle_params()

    def compute_bundle_params(self, d_vdw=None):
        """Compute bundle params."""

        super(NanotubeBundle, self).compute_tube_params()
        self._Ntubes = self.compute_Ntubes(nx=self._nx,
                                           ny=self._ny)
        self._bundle_mass = \
            self.compute_bundle_mass(n=self._n, m=self._m,
                                     nx=self._nx, ny=self._ny, nz=self._nz,
                                     with_units=self._with_units)
        self._Natoms_per_bundle = self._Ntubes * self._Natoms_per_tube
        if d_vdw is None:
            d_vdw = self._vdw_spacing

        if self._with_units and isinstance(d_vdw, float):
            d_vdw = Qty(d_vdw, 'angstrom')

        self._bundle_density = \
            self.compute_bundle_density(n=self._n, m=self._m,
                                        d_vdw=d_vdw, bond=self._bond,
                                        with_units=self._with_units)

    @property
    def Ntubes(self):
        """Number of nanotubes."""
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
        """Compute number of nanotubes.

        Parameters
        ----------
        nx, ny : int, optional
            Number of repeat unit cells in the x,y directions

        Returns
        -------
        int

        """
        return int(nx * ny)

    @property
    def Natoms_per_bundle(self):
        """Number of atoms in bundle.

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
        """Compute number of atoms per bundle.

        .. versionadded:: 0.2.5

        Parameters
        ----------
        n, m : int
        nz : {int, float}
        Ntubes : {None, int}, optional
        nx, ny : {None, int}, optional

        Returns
        -------
        int

        """
        Natoms_per_tube = \
            Nanotube.compute_Natoms_per_tube(n=n, m=m, nz=nz)
        if Ntubes is None and (nx is None or ny is None):
            raise NanotubeError("Ntubes or both nx and ny cells must be set")
        elif Ntubes is None and nx is not None and ny is not None:
            Ntubes = nx * ny
        Natoms_per_bundle = Ntubes * Natoms_per_tube

        return Natoms_per_bundle

    @property
    def bundle_mass(self):
        """Bundle mass in grams."""
        return self._bundle_mass

    @classmethod
    def compute_bundle_mass(cls, n=int, m=int, nx=int, ny=int, nz=None,
                            with_units=False):
        """Bundle mass in grams."""
        Ntubes = \
            NanotubeBundle.compute_Ntubes(nx=nx, ny=ny)
        tube_mass = Nanotube.compute_tube_mass(n=n, m=m, nz=nz,
                                               with_units=with_units)
        return Ntubes * tube_mass

    @property
    def bundle_density(self):
        """Bundle density in :math:`g/cm^3`."""
        return self._bundle_density

    @classmethod
    def compute_bundle_density(cls, n=int, m=int, d_vdw=None, bond=None,
                               element1=None, element2=None,
                               with_units=False):
        """Compute bundle mass.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.
        d_vdw : int
            van der Waals distance between nearest-neighbor tubes
        bond : float, optional
            Bond length.

        Returns
        -------
        float
            :math:`L_{\\mathrm{tube}}` in **nanometers**

        """
        if bond is None:
            bond = CCbond

        if isinstance(bond, float):
            bond = Qty(bond, 'angstroms')

        if d_vdw is None:
            if n == m:
                d_vdw = 3.38
            elif (m == 0) or (n == 0):
                d_vdw = 3.41
            else:
                d_vdw = 3.39

        if isinstance(d_vdw, float):
            d_vdw = Qty(d_vdw, 'angstroms')

        if element1 is None:
            element1 = 'C'
        if element2 is None:
            element2 = 'C'

        atom1 = Atom(element1)
        atom2 = Atom(element2)

        atom_mass = Qty(Atoms([atom1, atom2]).m, 'grams')

        bundle_density = 8 * np.pi**2 * atom_mass * \
            np.sqrt(n**2 + m**2 + n*m) / \
            (9 * np.sqrt(3) * (bond.to('cm'))**3 *
                (np.sqrt(n**2 + m**2 + n*m) +
                    np.pi * d_vdw / (np.sqrt(3) * bond))**2)
        if with_units:
            return bundle_density
        else:
            return bundle_density.magnitude
