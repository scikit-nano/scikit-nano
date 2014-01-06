# -*- coding: utf-8 -*-
"""
=============================================================
Nanotube structure tools (:mod:`sknano.nanogen.nanotube`)
=============================================================

.. currentmodule:: sknano.nanogen.nanotube

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

import copy
#import itertools
import warnings
warnings.filterwarnings('ignore')  # to suppress the Pint UnicodeWarning

from fractions import gcd
from collections import OrderedDict
from pint import UnitRegistry
ureg = UnitRegistry()
Qty = ureg.Quantity

import numpy as np

from pkshared.tools.arrayfuncs import rotation_matrix
from pkshared.tools.strfuncs import plural_word_check
from pkshared.tools.refdata import CCbond

from ..chemistry import Atom, Atoms
from ..structure_io import DATAWriter, XYZWriter, default_structure_format, \
    supported_structure_formats

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
           'Nanotube', 'NanotubeGenerator', 'NanotubeGeneratorError',
           'NanotubeBundle', 'NanotubeBundleGenerator',
           'MWNTGenerator']


class NanotubeGeneratorError(Exception):
    """Base class for NanotubeGenerator exceptions."""
    pass


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
    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1,
                 element1='C', element2='C',
                 bond=CCbond, tube_length=None,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False,
                 verbose=False):

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

        if isinstance(bond, float):
            self._bond = Qty(float(bond), 'angstroms')

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

        self._a = np.sqrt(3) * self._bond

        self._a1 = np.zeros(2, dtype=float)
        self._a2 = np.zeros(2, dtype=float)

        self._a1[0] = self._a2[0] = np.sqrt(3) / 2 * self._a.magnitude
        self._a1[1] = 1 / 2 * self._a.magnitude
        self._a2[1] = -self._a1[1]

        self._a1 = Qty(self._a1, 'angstrom')
        self._a2 = Qty(self._a2, 'angstrom')

        self._b1 = np.zeros(2, dtype=float)
        self._b2 = np.zeros(2, dtype=float)

        self._b1[0] = self._b2[0] = \
            1 / np.sqrt(3) * 2 * np.pi / self._a.magnitude
        self._b1[1] = 2 * np.pi / self._a.magnitude
        self._b2[1] = -self._b1[1]

        self._b1 = Qty(self._b1, '1/angstrom')
        self._b2 = Qty(self._b2, '1/angstrom')

        self._Ntubes = 1
        self._nx = int(nx)
        self._ny = int(ny)

        self._Lx = Lx
        self._Ly = Ly

        if Lz is not None:
            self._Lz = Qty(float(Lz), 'nanometers')
            self._T = self.compute_T(n=self._n, m=self._m, bond=self._bond,
                                     magnitude=False)
            if self._assume_integer_unit_cells:
                self._nz = \
                    int(np.ceil(self._Lz.to('angstroms').magnitude /
                                self._T.magnitude))
            else:
                self._nz = \
                    self._Lz.to('angstroms').magnitude / self._T.magnitude
        else:
            self._nz = int(nz)

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
        self._Ch = self.compute_Ch(n=self._n, m=self._m, bond=self._bond)
        self._T = self.compute_T(n=self._n, m=self._m, bond=self._bond)
        self._dt = self.compute_dt(n=self._n, m=self._m, bond=self._bond)
        self._rt = self.compute_rt(n=self._n, m=self._m, bond=self._bond)
        self._chiral_angle = self.compute_chiral_angle(n=self._n, m=self._m)
        self._M = self.compute_M(n=self._n, m=self._m)

        # Compute physical properties
        self._Lz = \
            self.compute_Lz(n=self._n, m=self._m, nz=self._nz, bond=self._bond)
        self._tube_mass = \
            self.compute_tube_mass(n=self._n, m=self._m, nz=self._nz)
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
                  length=False, magnitude=True):
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
            if ``True``, return the length of R
        magnitude : bool, optional
            if ``True``, return the length of R without units

        Returns
        -------
        (p, q) : tuple
            2-tuple of ints which are the components of R vector
        float
            length of R if ``length`` is ``True``

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

            if isinstance(bond, float):
                bond = Qty(bond, units)

            if magnitude:
                return bond.magnitude * np.sqrt(3 * (p**2 + q**2 + p * q))
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
    def compute_Ch(cls, n=int, m=int, bond=None, units='angstrom',
                   magnitude=True):
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

        if isinstance(bond, float):
            bond = Qty(bond, units)

        if magnitude:
            return bond.magnitude * np.sqrt(3 * (n**2 + m**2 + n * m))
        else:
            return bond * np.sqrt(3 * (n**2 + m**2 + n * m))

    @property
    def dt(self):
        """Nanotube diameter :math:`d_{t} = \\frac{|\\mathbf{C}_{h}|}{\\pi}`"""
        return self._dt

    @classmethod
    def compute_tube_diameter(cls, n=int, m=int, bond=None,
                              units='angstrom', magnitude=True):
        """Alias for :meth:`Nanotube.compute_dt`"""
        return Nanotube.compute_dt(n, m, bond=bond, units=units,
                                   magnitude=magnitude)

    @classmethod
    def compute_dt(cls, n=int, m=int, bond=None, units='angstrom',
                   magnitude=True):
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
        Ch = Nanotube.compute_Ch(n, m, bond=bond, units=units,
                                 magnitude=magnitude)
        return Ch / np.pi

    @property
    def rt(self):
        """Nanotube radius :math:`r_{t} = \\frac{|\\mathbf{C}_{h}|}{2\\pi}`"""
        return self._rt

    @classmethod
    def compute_tube_radius(cls, n=int, m=int, bond=None, units='angstrom',
                            magnitude=True):
        """Alias for :meth:`Nanotube.compute_rt`"""
        return Nanotube.compute_rt(n, m, bond=bond, units=units,
                                   magnitude=magnitude)

    @classmethod
    def compute_rt(cls, n=int, m=int, bond=None, units='angstrom',
                   magnitude=True):
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
        Ch = Nanotube.compute_Ch(n, m, bond=bond, units=units,
                                 magnitude=magnitude)
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
    def compute_chiral_angle(cls, n=int, m=int):
        """Compute chiral angle :math:`\\theta_{c}`

        .. math::

           \\theta_{c} = \\tan^{-1}\\left({\\frac{\\sqrt{3} m}{2n + m}}\\right)

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} +
            m\\mathbf{a}_{2} = (n, m)`.

        Returns
        -------
        float
            chiral angle :math:`\\theta_{c}` in degrees.

        """
        #return np.arccos((2*n + m) / (2 * np.sqrt(n**2 + m**2 + n*m)))
        return np.degrees(np.arctan(np.sqrt(3) * m / (2 * n + m)))

    @property
    def T(self):
        """Unit cell length :math:`|\\mathbf{T}|`.

        .. math::

           |\\mathbf{T}| = \\frac{\\sqrt{3} |\\mathbf{C}_{h}|}{d_{R}}

        """
        return self._T

    @classmethod
    def compute_T(cls, n=None, m=None, bond=None,
                  units='angstrom', magnitude=True):
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

        Raises
        ------
        NanotubeError
            if the parameters not valid.

        """
        if bond is None:
            bond = CCbond

        if isinstance(bond, float):
            bond = Qty(bond, units)

        Ch = Nanotube.compute_Ch(n=n, m=m, bond=bond, units=units,
                                 magnitude=magnitude)
        dR = Nanotube.compute_dR(n=n, m=m)
        return np.sqrt(3) * Ch / dR

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
    def compute_Lz(cls, n=int, m=int, nz=float, bond=None, units='angstrom',
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
        T = Nanotube.compute_T(n=n, m=m, bond=bond, units=units,
                               magnitude=False)
        Lz = nz * T
        Lz.ito('nanometer')
        if magnitude:
            return Lz.magnitude
        else:
            return Lz

    @property
    def tube_length(self):
        """Nanotube length :math:`L_{\\mathrm{tube}}` in **nanometers**."""
        return self.Lz

    @classmethod
    def compute_tube_length(cls, n=int, m=int, nz=float, bond=None,
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
        return Nanotube.compute_Lz(n=n, m=m, nz=nz, bond=bond,
                                   units=units, magnitude=magnitude)

    @property
    def tube_mass(self):
        """Nanotube mass in grams."""
        return self._tube_mass

    @classmethod
    def compute_tube_mass(cls, n=int, m=int, nz=float,
                          element1=None, element2=None,
                          units='grams', magnitude=True):
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

        mass = Qty(Natoms_per_tube * Atoms([atom1, atom2]).m, units)
        if magnitude:
            return mass.magnitude
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
                 Lx=None, Ly=None, Lz=None, fix_Lz=False, verbose=False):

        if tube_length is not None and Lz is None:
            Lz = tube_length

        super(NanotubeBundle, self).__init__(
            n=n, m=m, nx=ny, ny=ny, nz=nz,
            element1=element1, element2=element2,
            bond=bond, Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz,
            verbose=verbose)

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
        self._bundle_mass = self.compute_bundle_mass(n=self._n, m=self._m,
                                                     nx=self._nx,
                                                     ny=self._ny,
                                                     nz=self._nz)
        self._Natoms_per_bundle = self._Ntubes * self._Natoms_per_tube
        if d_vdw is None:
            d_vdw = self._vdw_spacing

        if isinstance(d_vdw, float):
            d_vdw = Qty(d_vdw, 'angstrom')

        self._bundle_density = \
            self.compute_bundle_density(n=self._n, m=self._m,
                                        d_vdw=d_vdw, bond=self._bond)

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
    def compute_bundle_mass(cls, n=int, m=int, nx=int, ny=int, nz=None):
        """Bundle mass in grams."""
        Ntubes = \
            NanotubeBundle.compute_Ntubes(nx=nx, ny=ny)
        tube_mass = Nanotube.compute_tube_mass(n=n, m=m, nz=nz)
        return Ntubes * tube_mass

    @property
    def bundle_density(self):
        """Bundle density in :math:`g/cm^3`."""
        return self._bundle_density

    @classmethod
    def compute_bundle_density(cls, n=int, m=int, d_vdw=None, bond=None,
                               element1=None, element2=None):
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
        return bundle_density


class NanotubeGenerator(Nanotube):
    u"""Class for generating nanotube structures.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nz : int, optional
        Number of repeat unit cells in the :math:`z` direction, along
        the *length* of the nanotube.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    Lz : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the ``nz`` value.

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

    autogen : bool, optional
        if ``True``, automatically call
        :py:meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by :py:meth:`~NanotubeGenerator.generate_structure_data`.
    verbose : bool, optional
        if ``True``, show verbose output

    Examples
    --------
    First, load the :py:class:`~sknano.nanogen.NanoGenerator` class.

    >>> from sknano.nanogen import NanotubeGenerator

    Now let's generate a :math:`\\mathbf{C}_{\\mathrm{h}} = (10, 5)`
    SWCNT unit cell.

    >>> nt = NanotubeGenerator(n=10, m=5)
    >>> nt.save_data(fname='10,5_unit_cell.xyz')

    The rendered structure looks like (orhographic view):

    .. image:: /images/10,5_unit_cell_orthographic_view.png

    and the perspective view:

    .. image:: /images/10,5_unit_cell_perspective_view.png

    """

    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1,
                 element1='C', element2='C',
                 bond=CCbond, Lx=None, Ly=None, Lz=None,
                 tube_length=None, fix_Lz=False,
                 autogen=True, verbose=False):

        if tube_length is not None and Lz is None:
            Lz = tube_length

        super(NanotubeGenerator, self).__init__(
            n=n, m=m, nx=nx, ny=ny, nz=nz,
            element1=element1, element2=element2,
            bond=bond, Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz,
            verbose=verbose)

        self._fname = None
        self.unit_cell = None
        self.structure_atoms = None

        if autogen:
            self.generate_unit_cell()
            self.generate_structure_data()

    @property
    def fname(self):
        """Structure file name."""
        return self._fname

    def generate_unit_cell(self):
        """Generate the unit cell."""
        n = self._n
        m = self._m
        t1 = self._t1
        t2 = self._t2
        Ch = self._Ch
        rt = self._rt
        N = self._N
        dR = self._dR

        e1 = self._element1
        e2 = self._element2

        self.unit_cell = Atoms()

        for q in xrange(t2, m + 1):
            for p in xrange(0, t1 + n + 1):
                M = m * p - n * q

                g_atom1 = Atom(e1)
                g_atom1.x = Ch * (q * t1 - p * t2) / N
                g_atom1.y = np.sqrt(3) * Ch * M / (N * dR)

                g_atom2 = Atom(e2)
                g_atom2.x = g_atom1.x + Ch * (n + m) / (N * dR)
                g_atom2.y = \
                    g_atom1.y - np.sqrt(3) * Ch * (n - m) / (3 * N * dR)

                phi1 = g_atom1.x / rt
                phi2 = g_atom2.x / rt

                if (g_atom1.x >= 0 and (q * t1 - p * t2) < N
                        and g_atom1.y >= 0 and (M < N)):

                    nt_atom1 = Atom(e1)

                    nt_atom1.x = rt * np.cos(phi1)
                    nt_atom1.y = rt * np.sin(phi1)
                    nt_atom1.z = g_atom1.y
                    self.unit_cell.append(nt_atom1)

                    nt_atom2 = Atom(e2)
                    nt_atom2.x = rt * np.cos(phi2)
                    nt_atom2.y = rt * np.sin(phi2)
                    nt_atom2.z = g_atom2.y
                    self.unit_cell.append(nt_atom2)

    def generate_structure_data(self):
        """Generate structure data."""
        self.structure_atoms = []
        for nz in xrange(int(np.ceil(self._nz))):
            dr = np.array([0.0, 0.0, nz * self.T])
            for uc_atom in self.unit_cell.atoms:
                nt_atom = Atom(uc_atom.symbol)
                nt_atom.r = uc_atom.r + dr
                self.structure_atoms.append(nt_atom)

    def save_data(self, fname=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  center_CM=True):
        """Save structure data.

        Parameters
        ----------
        fname : {None, str}, optional
            file name string
        structure_format : {None, str}, optional
            chemical file format of saved structure data. Must be one of:

                - xyz
                - data

            If ``None``, then guess based on ``fname`` file extension or
            default to ``xyz`` format.
        rotation_angle : {None, float}, optional
            Angle of rotation
        rot_axis : {'x', 'y', 'z'}, optional
            Rotation axis
        deg2rad : bool, optional
            Convert ``rotation_angle`` from degrees to radians.
        center_CM : bool, optional
            Center center-of-mass on origin.

        """
        if (fname is None and structure_format not in
                supported_structure_formats) or \
                (fname is not None and not
                    fname.endswith(supported_structure_formats) and
                    structure_format not in supported_structure_formats):
            structure_format = default_structure_format

        if fname is None:
            chirality = '{}{}r'.format('{}'.format(self._n).zfill(2),
                                       '{}'.format(self._m).zfill(2))
            if self._assume_integer_unit_cells:
                nz = ''.join(('{}'.format(self._nz),
                              plural_word_check('cell', self._nz)))
            else:
                nz = ''.join(('{:.2f}'.format(self._nz),
                              plural_word_check('cell', self._nz)))
            fname_wordlist = (chirality, nz)
            fname = '_'.join(fname_wordlist)
            fname += '.' + structure_format
        else:
            if fname.endswith(supported_structure_formats) and \
                    structure_format is None:
                for ext in supported_structure_formats:
                    if fname.endswith(ext):
                        structure_format = ext
                        break
            else:
                if structure_format is None or \
                        structure_format not in supported_structure_formats:
                    structure_format = default_structure_format

        #structure_atoms = list(itertools.chain(*self.structure_atoms))
        structure_atoms = None
        if isinstance(self.structure_atoms, list):
            structure_atoms = Atoms(self.structure_atoms)
        elif isinstance(self.structure_atoms, Atoms):
            structure_atoms = self.structure_atoms

        if center_CM:
            structure_atoms.center_CM()

        if self._L0 is not None and self._fix_Lz:
            structure_atoms.clip_bounds(abs_limit=(10 * self._L0 + 0.5) / 2,
                                        r_indices=[2])

        if rotation_angle is not None:
            R_matrix = rotation_matrix(rotation_angle,
                                       rot_axis=rot_axis,
                                       deg2rad=deg2rad)
            structure_atoms.rotate(R_matrix)

        if structure_format == 'data':
            DATAWriter.write(fname=fname, atoms=structure_atoms)
        else:
            XYZWriter.write(fname=fname, atoms=structure_atoms)

        self._fname = fname


class NanotubeBundleGenerator(NanotubeGenerator, NanotubeBundle):
    u"""Class for generating nanotube bundles.

    .. versionadded:: 0.2.4

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    vdw_spacing : float, optional
        van der Waals distance between nearest neighbor tubes

        .. versionadded:: 0.2.5

    bundle_packing : {None, 'hexagonal', 'cubic'}, optional
        close packing arrangement of bundles

        .. versionadded:: 0.2.5

    bundle_geometry : {None, 'triangle', 'hexagon', 'square', 'rectangle',
                       'rhombus', 'rhomboid'}, optional

        .. versionadded:: 0.2.5

    Lx, Ly, Lz : float, optional
        length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.

        .. versionadded:: 0.2.5

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If ``True``, then
        non integer :math:`n_z` cells are permitted.

        .. versionadded:: 0.2.6

    autogen : bool, optional
        if ``True``, automatically call
        :py:meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by :py:meth:`~NanotubeGenerator.generate_structure_data`.
    verbose : bool, optional
        if ``True``, show verbose output

    Examples
    --------

    Using the :py:class:`NanotubeBundleGenerator` class, you can
    generate cubic close packed (ccp) or hexagonal close packed
    bundles arrangements. In general, specifying **ccp** bundling will
    generate rectangular bundles (square bundles if :math:`n_x = n_y`)
    and specifying **hcp** bundling will generate *rhomboidal* bundles
    (*i.e.* bundles arranged within a rhomboid) (rhombuses if
    :math:`n_x = n_y`). However, you can also enforce a specific
    *bundle geometry* which will try and reshape the bundle arrangement so
    that it "fits inside" the boundaries of a specified geometric shape.
    This allows you to generate **hcp** bundles that are trianglar,
    hexagonal, or rectangular in *shape*, as some of the examples below
    illustrate.

    To start, let's generate an hcp bundle of
    :math:`C_{\\mathrm{h}} = (10, 5)` SWCNTs and cell count
    :math:`n_x=10, n_y=3, n_z=5`.

    >>> from sknano.nanogen import NanotubeBundleGenerator
    >>> SWCNTbundle = NanotubeBundleGenerator(n=10, m=5, nx=10,
    ...                                       ny=3, nz=5)
    >>> SWCNTbundle.save_data()

    The rendered structure looks like:

    .. image:: /images/1005_hcp_10cellsx3cellsx5cells-001.png

    Now let's generate a nice hexagon bundle, 3 tubes wide, with
    :math:`C_{\\mathrm{h}} = (6, 5)`.

    >>> SWCNTbundle = NanotubeBundleGenerator(n=6, m=5, nz=5,
    ...                                       bundle_geometry='hexagon')
    >>> SWCNTbundle.save_data()

    which looks like:

    .. image:: /images/0605_hcp_7tube_hexagon-001.png

    Remember, all of the :py:meth:`~NanotubeBundleGenerator.save_data`
    methods allow you to rotate the structure data before saving:

    >>> SWCNTbundle.save_data(fname='0605_hcp_7tube_hexagon_rot-30deg.xyz',
    ...                       rot_axis='z', rotation_angle=30)

    .. image:: /images/0605_hcp_7tube_hexagon_rot-30deg-001.png

    Now, just because we can, let's make a big ass hexagon bundle with
    :math:`C_{\\mathrm{h}} = (10, 0)`.

    >>> BIGASSHEXABUN = NanotubeBundleGenerator(n=10, m=0, nx=25,
    ...                                         ny=25, nz=1,
    ...                                         bundle_geometry='hexagon')

    You're looking at 469 :math:`(10, 0)` unit cells! That's
    :math:`N_{\\mathrm{atoms}} = 18760`.

    .. image:: /images/1000_hcp_469tube_hexagon-001.png

    """

    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond,
                 vdw_spacing=3.4, bundle_packing=None, bundle_geometry=None,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False,
                 autogen=True, verbose=False):

        super(NanotubeBundleGenerator, self).__init__(
            n=n, m=m, nx=nx, ny=ny, nz=nz,
            element1=element1, element2=element2,
            Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz,
            bond=bond, autogen=False,
            verbose=verbose)

        self.compute_bundle_params()

        self._r1 = np.zeros(3)
        self._r1[0] = Nanotube.compute_dt(n=n, m=m, bond=bond) + \
            vdw_spacing

        self._r2 = np.zeros(3)

        if bundle_packing is None and \
                bundle_geometry in ('square', 'rectangle'):
            bundle_packing = 'cubic'
        elif bundle_packing is None:
            bundle_packing = 'hexagonal'
        elif (bundle_packing == 'cubic' and bundle_geometry not in
                (None, 'square', 'rectangle')) or \
                (bundle_packing == 'hexagonal' and bundle_geometry not in
                    (None, 'triangle', 'hexagon', 'rhombus', 'rhomboid')):
            bundle_geometry = None

        if bundle_packing == 'cubic':
            self._r2[1] = self._r1[0]
        else:
            self._r2[0] = self._r1[0] * np.cos(np.pi / 3)
            self._r2[1] = self._r1[0] * np.sin(np.pi / 3)

        self._bundle_packing = bundle_packing
        self._bundle_geometry = bundle_geometry

        if autogen:
            super(NanotubeBundleGenerator, self).generate_unit_cell()
            self.generate_structure_data()

    def generate_structure_data(self):
        """Generate structure data."""
        super(NanotubeBundleGenerator, self).generate_structure_data()

        self._Ntubes = 0

        swcnt0 = copy.deepcopy(self.structure_atoms)
        self.structure_atoms = Atoms()
        if self._bundle_geometry == 'hexagon':
            nrows = max(self._nx, self._ny, 3)
            if nrows % 2 != 1:
                nrows += 1

            ntubes_per_end_rows = int((nrows + 1) / 2)

            row = 0
            ntubes_per_row = nrows
            while ntubes_per_row >= ntubes_per_end_rows:
                if row == 0:
                    for n in xrange(ntubes_per_row):
                        swcnt = Atoms(atoms=swcnt0, deepcopy=True)
                        swcnt.center_CM()
                        dr = n * self._r1
                        swcnt.translate(dr)
                        self.structure_atoms.extend(swcnt.atoms)
                        self._Ntubes += 1
                else:
                    for nx in xrange(ntubes_per_row):
                        for ny in (-row, row):
                            swcnt = Atoms(atoms=swcnt0, deepcopy=True)
                            swcnt.center_CM()
                            dy = np.zeros(3)
                            dy[0] = abs(ny) * self._r2[0]
                            dy[1] = ny * self._r2[1]
                            dr = nx * self._r1 + dy
                            swcnt.translate(dr)
                            self.structure_atoms.extend(swcnt.atoms)
                            self._Ntubes += 1
                row += 1
                ntubes_per_row = nrows - row
        else:
            for nx in xrange(self._nx):
                for ny in xrange(self._ny):
                    swcnt = Atoms(atoms=swcnt0, deepcopy=True)
                    swcnt.center_CM()
                    dr = nx * self._r1 + ny * self._r2
                    swcnt.translate(dr)
                    self.structure_atoms.extend(swcnt.atoms)
                    self._Ntubes += 1
        self._Natoms_per_bundle = \
            self.compute_Natoms_per_bundle(n=self._n, m=self._m,
                                           nz=self._nz,
                                           Ntubes=self._Ntubes)

    def save_data(self, fname=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  center_CM=True):
        """Save structure data.

        Parameters
        ----------
        fname : {None, str}, optional
            file name string
        structure_format : {None, str}, optional
            chemical file format of saved structure data. Must be one of:

                - xyz
                - data

            If ``None``, then guess based on ``fname`` file extension or
            default to ``xyz`` format.
        rotation_angle : {None, float}, optional
            Angle of rotation
        rot_axis : {'x', 'y', 'z'}, optional
            Rotation axis
        deg2rad : bool, optional
            Convert ``rotation_angle`` from degrees to radians.
        center_CM : bool, optional
            Center center-of-mass on origin.

        """
        if (fname is None and structure_format not in
                supported_structure_formats) or \
                (fname is not None and not
                    fname.endswith(supported_structure_formats) and
                    structure_format not in supported_structure_formats):
            structure_format = default_structure_format

        if fname is None:

            chirality = '{}{}'.format('{}'.format(self._n).zfill(2),
                                      '{}'.format(self._m).zfill(2))
            packing = '{}cp'.format(self._bundle_packing[0])
            #Ntubes = ''.join(('{}'.format(self._Ntubes),
            #                  plural_word_check('tube', self._Ntubes)))
            Ntube = '{}tube'.format(self._Ntubes)

            fname_wordlist = None
            if self._bundle_geometry is None:
                nx = ''.join(('{}'.format(self._nx),
                             plural_word_check('cell', self._nx)))
                ny = ''.join(('{}'.format(self._ny),
                             plural_word_check('cell', self._ny)))
                if self._assume_integer_unit_cells:
                    nz = ''.join(('{}'.format(self._nz),
                                  plural_word_check('cell', self._nz)))
                else:
                    nz = ''.join(('{:.2f}'.format(self._nz),
                                  plural_word_check('cell', self._nz)))
                cells = 'x'.join((nx, ny, nz))
                fname_wordlist = (chirality, packing, cells)
            else:
                fname_wordlist = \
                    (chirality, packing, Ntube, self._bundle_geometry)

            fname = '_'.join(fname_wordlist)
            fname += '.' + structure_format

        else:
            if fname.endswith(supported_structure_formats) and \
                    structure_format is None:
                for ext in supported_structure_formats:
                    if fname.endswith(ext):
                        structure_format = ext
                        break
            else:
                if structure_format is None or \
                        structure_format not in supported_structure_formats:
                    structure_format = default_structure_format

        super(NanotubeBundleGenerator, self).save_data(
            fname=fname, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            deg2rad=deg2rad, center_CM=center_CM)


class MWNTGenerator(NanotubeGenerator, NanotubeBundle):
    u"""Class for generating multi-walled nanotubes.

    .. versionadded:: 0.2.7

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    vdw_spacing : float, optional
        van der Waals distance between nearest neighbor tubes
    bundle_packing : {None, 'hexagonal', 'cubic'}, optional
        close packing arrangement of bundles
    bundle_geometry : {None, 'triangle', 'hexagon', 'square', 'rectangle',
                       'rhombus', 'rhomboid'}, optional
    Lx, Ly, Lz : float, optional
        length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.
    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If ``True``, then
        non integer :math:`n_z` cells are permitted.
    autogen : bool, optional
        if ``True``, automatically call
        :py:meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by :py:meth:`~NanotubeGenerator.generate_structure_data`.
    verbose : bool, optional
        if ``True``, show verbose output

    Examples
    --------

    >>> from sknano.nanogen import MWNTGenerator
    >>> mwnt = MWNTGenerator(n=40, m=40, max_shells=5, Lz=1.0, fix_Lz=True)
    >>> mwnt.save_data()

    .. image:: /images/5shell_mwnt_4040_outer_Ch_1cellx1cellx4.06cells-001.png

    >>> mwntbundle = MWNTGenerator(n=40, m=40, max_shells=5, Lz=1.0,
    ...                            fix_Lz=True, bundle_geometry='hexagon')
    >>> mwntbundle.save_data()

    .. image:: /images/5shell_mwnt_4040_outer_Ch_hcp_7tube_hexagon-001.png

    """
    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond,
                 vdw_spacing=3.4, bundle_packing=None, bundle_geometry=None,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False,
                 max_shells=None, min_shell_diameter=0.0,
                 shell_spacing=3.4, inner_shell_Ch_type=None,
                 autogen=True, verbose=False):

        super(MWNTGenerator, self).__init__(
            n=n, m=m, nx=nx, ny=ny, nz=nz,
            element1=element1, element2=element2,
            Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz,
            bond=bond, autogen=False,
            verbose=verbose)

        self._Lzmin = np.inf

        self._max_shells = max_shells
        if max_shells is None:
            self._max_shells = np.inf

        self._min_shell_diameter = min_shell_diameter
        self._shell_spacing = shell_spacing
        self._inner_shell_Ch_type = inner_shell_Ch_type

        self._Nshells_per_tube = 1
        self._Natoms_per_tube = 0

        self.compute_bundle_params()

        self._r1 = np.zeros(3)
        self._r1[0] = Nanotube.compute_dt(n=n, m=m, bond=bond) + \
            vdw_spacing

        self._r2 = np.zeros(3)

        if bundle_packing is None and \
                bundle_geometry in ('square', 'rectangle'):
            bundle_packing = 'cubic'
        elif bundle_packing is None:
            bundle_packing = 'hexagonal'
        elif (bundle_packing == 'cubic' and bundle_geometry not in
                (None, 'square', 'rectangle')) or \
                (bundle_packing == 'hexagonal' and bundle_geometry not in
                    (None, 'triangle', 'hexagon', 'rhombus', 'rhomboid')):
            bundle_geometry = None

        if bundle_packing == 'cubic':
            self._r2[1] = self._r1[0]
        else:
            self._r2[0] = self._r1[0] * np.cos(np.pi / 3)
            self._r2[1] = self._r1[0] * np.sin(np.pi / 3)

        self._bundle_packing = bundle_packing
        self._bundle_geometry = bundle_geometry

        if autogen:
            super(MWNTGenerator, self).generate_unit_cell()
            self.generate_structure_data()

    def generate_unit_cell(self, n=int, m=int):
        """Generate the unit cell."""
        t1 = Nanotube.compute_t1(n=n, m=m)
        t2 = Nanotube.compute_t2(n=n, m=m)
        Ch = Nanotube.compute_Ch(n=n, m=m)
        rt = Nanotube.compute_rt(n=n, m=m)
        N = Nanotube.compute_N(n=n, m=m)
        dR = Nanotube.compute_dR(n=n, m=m)

        e1 = self._element1
        e2 = self._element2

        unit_cell = Atoms()

        for q in xrange(t2, m + 1):
            for p in xrange(0, t1 + n + 1):
                M = m * p - n * q

                g_atom1 = Atom(e1)
                g_atom1.x = Ch * (q * t1 - p * t2) / N
                g_atom1.y = np.sqrt(3) * Ch * M / (N * dR)

                g_atom2 = Atom(e2)
                g_atom2.x = g_atom1.x + Ch * (n + m) / (N * dR)
                g_atom2.y = \
                    g_atom1.y - np.sqrt(3) * Ch * (n - m) / (3 * N * dR)

                phi1 = g_atom1.x / rt
                phi2 = g_atom2.x / rt

                if (g_atom1.x >= 0 and (q * t1 - p * t2) < N
                        and g_atom1.y >= 0 and (M < N)):

                    nt_atom1 = Atom(e1)

                    nt_atom1.x = rt * np.cos(phi1)
                    nt_atom1.y = rt * np.sin(phi1)
                    nt_atom1.z = g_atom1.y
                    unit_cell.append(nt_atom1)

                    nt_atom2 = Atom(e2)
                    nt_atom2.x = rt * np.cos(phi2)
                    nt_atom2.y = rt * np.sin(phi2)
                    nt_atom2.z = g_atom2.y
                    unit_cell.append(nt_atom2)

        return unit_cell

    def generate_structure_data(self):
        """Generate structure data."""
        super(MWNTGenerator, self).generate_structure_data()

        self._Ntubes = 0

        dt = []
        Ch = []
        for n in xrange(0, 201):
            for m in xrange(0, 201):
                dt.append(Nanotube.compute_dt(n=n, m=m))
                Ch.append((n, m))
        dt = np.asarray(dt)
        Ch = np.asarray(Ch)

        swnt0 = copy.deepcopy(self.structure_atoms)
        mwnt0 = Atoms(atoms=swnt0, deepcopy=True)
        self._Lzmin = min(self._Lzmin, self._Lz)
        mwnt0.center_CM()

        next_dt = self.dt - 2 * self._shell_spacing
        while next_dt >= self._min_shell_diameter and \
                self._Nshells_per_tube < self._max_shells:
            # get chiral indices for next_dt
            next_Ch_candidates = []
            delta_dt = 0.05
            while len(next_Ch_candidates) == 0:
                if self._inner_shell_Ch_type in ('AC', 'armchair'):
                    next_Ch_candidates = \
                        Ch[np.where(
                            np.logical_and(np.abs(dt - next_dt) <= delta_dt,
                                           Ch[:,0] == Ch[:,1]))]
                elif self._inner_shell_Ch_type in ('ZZ', 'zigzag'):
                    next_Ch_candidates = \
                        Ch[np.where(
                            np.logical_and(np.abs(dt - next_dt) <= delta_dt,
                                           np.logical_or(Ch[:,0] == 0,
                                                         Ch[:,1] == 0)))]
                elif self._inner_shell_Ch_type == 'achiral':
                    next_Ch_candidates = \
                        Ch[np.where(
                            np.logical_and(np.abs(dt - next_dt) <= delta_dt,
                                           np.logical_or(
                                               Ch[:,0] == Ch[:,1],
                                               np.logical_or(
                                                   Ch[:,0] == 0,
                                                   Ch[:,1] == 0))))]
                elif self._inner_shell_Ch_type == 'chiral':
                    next_Ch_candidates = \
                        Ch[np.where(
                            np.logical_and(np.abs(dt - next_dt) <= delta_dt,
                                           np.logical_and(
                                               Ch[:,0] != Ch[:,1],
                                               np.logical_and(
                                                   Ch[:,0] != 0,
                                                   Ch[:,1] != 1))))]
                else:
                    next_Ch_candidates = \
                        Ch[np.where(np.abs(dt - next_dt) <= delta_dt)]

                next_dt -= delta_dt

                #delta_dt += 0.05

            n, m = \
                next_Ch_candidates[
                    np.random.choice(np.arange(len(next_Ch_candidates)))]
            # generate unit cell for new chiral indices
            unit_cell = self.generate_unit_cell(n=n, m=m)
            if self._verbose:
                print('next_dt: {:.4f}'.format(next_dt))
                print('n, m = {}, {}'.format(n, m))
                print('unit_cell.Natoms: {}\n'.format(unit_cell.Natoms))
            T = Nanotube.compute_T(n=n, m=m)
            Lz = Nanotube.compute_Lz(n=n, m=m, nz=self._nz)
            self._Lzmin = min(self._Lzmin, Lz)
            shell_atoms = Atoms()
            for nz in xrange(int(np.ceil(self._nz))):
                dr = np.array([0.0, 0.0, nz * T])
                for uc_atom in unit_cell.atoms:
                    nt_atom = Atom(uc_atom.symbol)
                    nt_atom.r = uc_atom.r + dr
                    shell_atoms.append(nt_atom)
            shell_atoms.center_CM()
            mwnt0.extend(shell_atoms.atoms)
            next_dt -= 2 * self._shell_spacing
            self._Nshells_per_tube += 1

        if self._L0 is not None and self._fix_Lz:
            mwnt0.clip_bounds(abs_limit=(10 * self._L0 + 0.5) / 2,
                              r_indices=[2])
        else:
            mwnt0.clip_bounds(abs_limit=(10 * self._Lzmin + 0.5) / 2,
                              r_indices=[2])

        self._Natoms_per_tube = mwnt0.Natoms

        self.structure_atoms = Atoms()

        if self._bundle_geometry == 'hexagon':
            nrows = max(self._nx, self._ny, 3)
            if nrows % 2 != 1:
                nrows += 1

            ntubes_per_end_rows = int((nrows + 1) / 2)

            row = 0
            ntubes_per_row = nrows
            while ntubes_per_row >= ntubes_per_end_rows:
                if row == 0:
                    for n in xrange(ntubes_per_row):
                        mwnt = Atoms(atoms=mwnt0.atoms, deepcopy=True)
                        mwnt.center_CM()
                        dr = n * self._r1
                        mwnt.translate(dr)
                        self.structure_atoms.extend(mwnt.atoms)
                        self._Ntubes += 1
                else:
                    for nx in xrange(ntubes_per_row):
                        for ny in (-row, row):
                            mwnt = Atoms(atoms=mwnt0.atoms, deepcopy=True)
                            mwnt.center_CM()
                            dy = np.zeros(3)
                            dy[0] = abs(ny) * self._r2[0]
                            dy[1] = ny * self._r2[1]
                            dr = nx * self._r1 + dy
                            mwnt.translate(dr)
                            self.structure_atoms.extend(mwnt.atoms)
                            self._Ntubes += 1
                row += 1
                ntubes_per_row = nrows - row
        else:
            for nx in xrange(self._nx):
                for ny in xrange(self._ny):
                    mwnt = Atoms(atoms=mwnt0.atoms, deepcopy=True)
                    mwnt.center_CM()
                    dr = nx * self._r1 + ny * self._r2
                    mwnt.translate(dr)
                    self.structure_atoms.extend(mwnt.atoms)
                    self._Ntubes += 1
        self._Natoms_per_bundle = self._Ntubes * self._Natoms_per_tube

        if self._verbose:
            print('Ntubes: {}'.format(self._Ntubes))
            print('Nshells_per_tube: {}'.format(self._Nshells_per_tube))
            print('Natoms_per_tube: {}'.format(self._Natoms_per_tube))
            print('Natoms_per_bundle: {}'.format(self._Natoms_per_bundle))

    def save_data(self, fname=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  center_CM=True):
        """Save structure data.

        Parameters
        ----------
        fname : {None, str}, optional
            file name string
        structure_format : {None, str}, optional
            chemical file format of saved structure data. Must be one of:

                - xyz
                - data

            If ``None``, then guess based on ``fname`` file extension or
            default to ``xyz`` format.
        rotation_angle : {None, float}, optional
            Angle of rotation
        rot_axis : {'x', 'y', 'z'}, optional
            Rotation axis
        deg2rad : bool, optional
            Convert ``rotation_angle`` from degrees to radians.
        center_CM : bool, optional
            Center center-of-mass on origin.

        """
        if (fname is None and structure_format not in
                supported_structure_formats) or \
                (fname is not None and not
                    fname.endswith(supported_structure_formats) and
                    structure_format not in supported_structure_formats):
            structure_format = default_structure_format

        if fname is None:

            Nshells = '{}shell_mwnt'.format(self._Nshells_per_tube)

            chirality = '{}{}_outer_Ch'.format('{}'.format(self._n).zfill(2),
                                               '{}'.format(self._m).zfill(2))
            packing = '{}cp'.format(self._bundle_packing[0])

            Ntube = '{}tube'.format(self._Ntubes)

            fname_wordlist = None
            if self._bundle_geometry is None:
                nx = ''.join(('{}'.format(self._nx),
                             plural_word_check('cell', self._nx)))
                ny = ''.join(('{}'.format(self._ny),
                             plural_word_check('cell', self._ny)))
                if self._assume_integer_unit_cells:
                    nz = ''.join(('{}'.format(self._nz),
                                  plural_word_check('cell', self._nz)))
                else:
                    nz = ''.join(('{:.2f}'.format(self._nz),
                                  plural_word_check('cell', self._nz)))
                cells = 'x'.join((nx, ny, nz))

                if self._nx == 1 and self._ny == 1:
                    fname_wordlist = (Nshells, chirality, cells)
                else:
                    fname_wordlist = (Nshells, chirality, packing, cells)
            else:
                fname_wordlist = \
                    (Nshells, chirality, packing, Ntube, self._bundle_geometry)

            fname = '_'.join(fname_wordlist)
            fname += '.' + structure_format

        else:
            if fname.endswith(supported_structure_formats) and \
                    structure_format is None:
                for ext in supported_structure_formats:
                    if fname.endswith(ext):
                        structure_format = ext
                        break
            else:
                if structure_format is None or \
                        structure_format not in supported_structure_formats:
                    structure_format = default_structure_format

        super(MWNTGenerator, self).save_data(
            fname=fname, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            deg2rad=deg2rad, center_CM=center_CM)
