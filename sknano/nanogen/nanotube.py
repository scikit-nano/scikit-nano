# -*- coding: utf-8 -*-
"""
=============================================================
Nanotube structure tools (:mod:`sknano.nanogen.nanotube`)
=============================================================

.. currentmodule:: sknano.nanogen.nanotube

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

from fractions import gcd
from math import pi
from collections import OrderedDict
#import itertools

import numpy as np

from pksci.chemistry import Atom, Atoms
from pkshared.tools.arrayfuncs import rotation_matrix
from pkshared.tools.strfuncs import plural_word_check
from pkshared.tools.refdata import ccbond

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
           'Nanotube', 'NanotubeGenerator',
           'NanotubeBundle', 'NanotubeBundleGenerator',
           'SWNTGenerator', 'MWNTGenerator']


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
    nzcells : int, optional
        Number of repeat unit cells in the :math:`z` direction, along
        the *length* of the nanotube.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        bond length between nearest neighbor atoms.
        Must be in units of **Angstroms**. Default value is
        the carbon-carbon bond length in graphite:
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.421` \u212b ([SoFaCNTs]_)
    tube_length : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the ``nzcells`` value.
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

    def __init__(self, n=int, m=int, nzcells=1, element1='C', element2='C',
                 bond=ccbond, tube_length=None, verbose=False):

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

        self._n = int(n)
        self._m = int(m)
        self._element1 = element1
        self._element2 = element2

        self._bond = float(bond)
        self._verbose = verbose

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
        self._tube_length = None
        self._Natoms_per_tube = None

        if tube_length is not None:
            self._tube_length = float(tube_length)
            self._nzcells = int(np.ceil(10 * self._tube_length / self._T))
        else:
            self._nzcells = int(nzcells)

        for k, v in self.__dict__.iteritems():
            p = k.strip('_')
            if p in self._params.keys():
                self._params[p]['units'] = param_units.get(p)
                self._params[p]['strfmt'] = param_strfmt.get(p)
                if param_symbols.get(p) is not None:
                    self._params[p]['var'] = param_symbols[p]
                else:
                    self._params[p]['var'] = p

        self._compute_tube_params()

    def _compute_tube_params(self):
        self._d = self.compute_d(n=self._n, m=self._m)
        self._dR = self.compute_dR(n=self._n, m=self._m)
        self._t1 = self.compute_t1(n=self._n, m=self._m)
        self._t2 = self.compute_t2(n=self._n, m=self._m)

        #Compute geometric properties
        self._Ch = self.compute_Ch(n=self._n, m=self._m, bond=self._bond)
        self._T = self.compute_T(n=self._n, m=self._m, bond=self._bond)
        self._dt = self.compute_dt(n=self._n, m=self._m, bond=self._bond)
        self._rt = self.compute_rt(n=self._n, m=self._m, bond=self._bond)
        self._chiral_angle = self.compute_chiral_angle(n=self._n, m=self._m)
        self._tube_length = \
            self.compute_tube_length(n=self._n,
                                     m=self._m,
                                     nzcells=self._nzcells,
                                     bond=self._bond)

        self._M = self.compute_M(n=self._n, m=self._m)

        # Compute symmetry properties
        self._R = self.compute_R(n=self._n, m=self._m)

        #Compute atomistic properties
        self._N = self.compute_N(n=self._n, m=self._m)
        self._Natoms = self.compute_Natoms(n=self._n, m=self._m)
        self._Natoms_per_tube = \
            self.compute_Natoms_per_tube(n=self._n,
                                         m=self._m,
                                         nzcells=self._nzcells)

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

    @property
    def n(self):
        """Chiral index :math:`n`"""
        return self._n

    @n.setter
    def n(self, value):
        self._n = int(value)
        self._compute_tube_params()

    @property
    def m(self):
        """Chiral index :math:`m`"""
        return self._m

    @m.setter
    def m(self, value):
        self._m = int(value)
        self._compute_tube_params()

    @property
    def bond(self):
        """Bond length in **Angstroms**."""
        return self._bond

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
        """Number of hexagons per nanotube unit cell :math:`N`."""
        return self._N

    @classmethod
    def compute_N(cls, n=int, m=int):
        """Compute :math:`N = \\frac{2(n^2+m^2+nm)}{d_{R}}`.
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.

        Returns
        -------
        int
            Number of hexagons per nanotube unit cell :math:`N`.

        """
        dR = Nanotube.compute_dR(n=n, m=m)
        return int(2 * (n**2 + m**2 + n * m) / dR)

    @property
    def Natoms(self):
        """Number of atoms per nanotube **unit cell** :math:`2N`."""

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
            Number of atoms per nanotube **unit cell** :math:`2N`.

        """
        N = Nanotube.compute_N(n=n, m=m)
        return 2 * N

    @property
    def R(self):
        """Symmetry vector :math:`\\mathbf{R} = (p, q)`"""
        return self._R

    @classmethod
    def compute_R(cls, n=int, m=int, bond=None, magnitude=False):
        """Compute symmetry vector :math:`\\mathbf{R} = (p, q)`

        .. math::

           \\mathbf{R} = p\\mathbf{a}_{1} + q\\mathbf{a}_{2}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.
        magnitude : bool, optional
            if ``True``, return the magnitude of R

        Returns
        -------
        (p, q) : tuple
            2-tuple of ints which are the components of R vector
        float
            magnitude of R if ``magnitude`` is `True`

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

        if magnitude:
            if bond is None:
                bond = ccbond
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
    def compute_Ch(cls, n=int, m=int, bond=None):
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
            bond = ccbond

        return bond * np.sqrt(3 * (n**2 + m**2 + n * m))

    @property
    def dt(self):
        """Nanotube diameter :math:`d_{t} = \\frac{|\\mathbf{C}_{h}|}{\\pi}`"""
        return self._dt

    @classmethod
    def compute_tube_diameter(cls, n=int, m=int, bond=None):
        """Alias for :meth:`Nanotube.compute_dt`"""
        return Nanotube.compute_dt(n, m, bond)

    @classmethod
    def compute_dt(cls, n=int, m=int, bond=None):
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
        Ch = Nanotube.compute_Ch(n, m, bond)
        return Ch / pi

    @property
    def rt(self):
        """Nanotube radius :math:`r_{t} = \\frac{|\\mathbf{C}_{h}|}{2\\pi}`"""
        return self._rt

    @classmethod
    def compute_tube_radius(cls, n=int, m=int, bond=None):
        """Alias for :meth:`Nanotube.compute_rt`"""
        return Nanotube.compute_rt(n, m, bond)

    @classmethod
    def compute_rt(cls, n=int, m=int, bond=None):
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
        Ch = Nanotube.compute_Ch(n, m, bond)
        return Ch / (2 * pi)

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

           \\theta_{c} = \\atan{\\frac{\\sqrt{3} m}{2n + m}}

        """
        return self._chiral_angle

    @classmethod
    def compute_chiral_angle(cls, n=int, m=int):
        """Compute chiral angle :math:`\\theta_{c}`

        .. math::

           \\theta_{c} = \\atan{\\frac{\\sqrt{3} m}{2n + m}}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} +
            m\\mathbf{a}_{2} = (n, m)`.

        """
        #return np.arccos((2*n + m) / (2 * np.sqrt(n**2 + m**2 + n*m)))
        return np.degrees(np.arctan(np.sqrt(3) * m / (2 * n + m)))

    @property
    def T(self):
        """Unit cell length.

        :math:`|\\mathbf{T}| = \\frac{\\sqrt{3} |\\mathbf{C}_{h}|}{d_{R}}`

        """
        return self._T

    @classmethod
    def compute_T(cls, n=None, m=None, Ch=None, dR=None, bond=None):
        """Compute unit cell length :math:`|\\mathbf{T}|`

        .. math::

           |\\mathbf{T}| = \\frac{\\sqrt{3} |\\mathbf{C}_{h}|}{d_{R}}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} +
            m\\mathbf{a}_{2} = (n, m)`.
        Ch : float, optional
            nanotube circumference in Angstroms
        dR : int, optional
            greatest common divisor of :math:`2n + m` :math:`2m + n`
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
            bond = ccbond
        if n is not None and m is not None:
            Ch = Nanotube.compute_Ch(n=n, m=m, bond=bond)
            dR = Nanotube.compute_dR(n=n, m=m)
            return np.sqrt(3) * Ch / dR
        elif Ch is not None and dR is not None:
            return np.sqrt(3) * Ch / dR
        else:
            raise NanotubeError("Invalid parameters")

    @property
    def Natoms_per_tube(self):
        """Number of atoms per nanotube :math:`N_{\\mathrm{atoms/tube}}`."""
        return self._Natoms_per_tube

    @classmethod
    def compute_Natoms_per_tube(cls, n=int, m=int, nzcells=int):
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
        nzcells : int
            Number of nanotube unit cells

        Returns
        -------
        int
            :math:`N_{\\mathrm{atoms/tube}}`
        """
        Natoms = Nanotube.compute_Natoms(n=n, m=m)
        return int(Natoms * nzcells)

    @property
    def nzcells(self):
        """Number of nanotube unit cells along the :math:`z`-axis."""
        return self._nzcells

    @property
    def tube_length(self):
        """Nanotube length :math:`L_{\\mathrm{tube}}` in **nanometers**."""
        return self._tube_length

    @classmethod
    def compute_tube_length(cls, n=int, m=int, nzcells=int, bond=None):
        """Compute :math:`L_{\\mathrm{tube}}`.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.
        nzcells : int
            Number of nanotube unit cells
        bond : float, optional
            bond length

        Returns
        -------
        float
            :math:`L_{\\mathrm{tube}}` in **nanometers**

        """
        T = Nanotube.compute_T(n=n, m=m, bond=bond)
        return nzcells * T / 10.

    def _ncell_check(self, nzcells=None):
        """Check unit cell count.

        Parameters
        ----------
        nzcells : int
            the number of unit cells

        Raises
        ------
        CellError
            if unit cell count is invalid.

        """
        #if (nzcells <= 0) or (int(nzcells) != nzcells):
        #    raise CellError('unit cell count must be positive integer')
        if (nzcells <= 0):
            raise CellError('unit cell count must be positive')

    def _tube_length_check(self, tube_length=None):
        """Check tube length.

        Parameters
        ----------
        tube_length : float
            the length of the nanotube in nanometers

        Raises
        ------
        TubeLengthError
            if tube_length is invalid.

        """
        if (float(tube_length) <= 0):
            raise TubeLengthError('nanotube length must be positive')


class NanotubeBundle(Nanotube):
    u"""Class for creating interactive Nanotube *Bundle* objects.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nxcells, nycells, nzcells : int, optional
        Number of repeat unit cells in the x,y,z directions
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        bond length between nearest neighbor atoms.
        Must be in units of **Angstroms**. Default value is
        the carbon-carbon bond length in graphite:
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.421` \u212b ([SoFaCNTs]_)
    tube_length : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the ``nzcells`` value.
    verbose : bool, optional
        verbose output

    """
    def __init__(self, n=int, m=int, nxcells=1, nycells=1, nzcells=1,
                 element1='C', element2='C', bond=ccbond, tube_length=None,
                 verbose=False):

        super(NanotubeBundle, self).__init__(
            n=n, m=m, nzcells=nzcells, element1=element1, element2=element2,
            bond=bond, tube_length=tube_length, verbose=verbose)

        self._nxcells = int(nxcells)
        self._nycells = int(nycells)
        self._Ntubes = self._nxcells * self._nycells

    @property
    def nxcells(self):
        """Number of nanotube unit cells along the :math:`x`-axis."""
        return int(self._nxcells)

    @property
    def nycells(self):
        """Number of nanotube unit cells along the :math:`y`-axis."""
        return int(self._nycells)

    @property
    def Ntubes(self):
        """Number of nanotubes."""
        return int(self._Ntubes)


class NanotubeGenerator(Nanotube):
    u"""Class for generating nanotube structures.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nzcells : int, optional
        Number of repeat unit cells in the :math:`z` direction, along
        the *length* of the nanotube.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    tube_length : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the ``nzcells`` value.
    autogen : bool, optional
        if ``True``, automatically call
        :py:meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by
        :py:meth:`~NanotubeGenerator.generate_structure_data`.
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

    def __init__(self, n=int, m=int, nzcells=1, element1='C', element2='C',
                 bond=ccbond, tube_length=None, autogen=True, verbose=False):

        super(NanotubeGenerator, self).__init__(
            n=n, m=m, nzcells=nzcells, element1=element1, element2=element2,
            bond=bond, tube_length=tube_length, verbose=verbose)

        self.unit_cell = None
        self.structure_atoms = None

        if autogen:
            self.generate_unit_cell()
            self.generate_structure_data()

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
        for nz in xrange(self._nzcells):
            dr = np.array([0.0, 0.0, nz * self.T])
            for uc_atom in self.unit_cell.atomlist:
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
            Center center-of-mass on on origin.

        """
        if (fname is None and structure_format not in
                supported_structure_formats) or \
                (fname is not None and not
                    fname.endswith(supported_structure_formats) and
                    structure_format not in supported_structure_formats):
            structure_format = default_structure_format

        #structure_atoms = list(itertools.chain(*self.structure_atoms))
        structure_atoms = Atoms(self.structure_atoms)
        if center_CM:
            structure_atoms.center_CM(r_indices=[2])
        if rotation_angle is not None:
            R_matrix = rotation_matrix(rotation_angle,
                                       rot_axis=rot_axis,
                                       deg2rad=deg2rad)
            structure_atoms.rotate(R_matrix)
        if fname is None:
            chirality = '{}{}r'.format('{}'.format(self._n).zfill(2),
                                       '{}'.format(self._m).zfill(2))
            nzcells = ''.join(('{}'.format(self._nzcells),
                               plural_word_check('cell', self._nzcells)))
            fname_wordlist = (chirality, nzcells)
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
        if structure_format == 'data':
            DATAWriter.write(fname=fname, atoms=structure_atoms)
        else:
            XYZWriter.write(fname=fname, atoms=structure_atoms)


class NanotubeBundleGenerator(NanotubeGenerator):
    u"""Class for generating nanotube bundles.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nxcells, nycells, nzcells : int, optional
        Number of repeat unit cells in the x,y,z directions
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    tube_length : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the ``nzcells`` value.
    autogen : bool, optional
        if ``True``, automatically call
        :py:meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by
        :py:meth:`~NanotubeGenerator.generate_structure_data`.
    verbose : bool, optional
        if ``True``, show verbose output

    Examples
    --------

    """

    def __init__(self, n=int, m=int, nxcells=1, nycells=1, nzcells=1,
                 element1='C', element2='C', bond=ccbond, tube_length=None,
                 autogen=True, verbose=False):

        super(NanotubeBundleGenerator, self).__init__(
            n=n, m=m, nxcells=nxcells, nycells=nycells, nzcells=nzcells,
            element1=element1, element2=element2, bond=bond,
            tube_length=tube_length, verbose=verbose)

        if autogen:
            super(NanotubeBundleGenerator, self).generate_unit_cell()
            self.generate_structure_data()

    def generate_structure_data(self):
        """Generate structure data."""
        super(NanotubeBundleGenerator, self).generate_structure()
        for nx in xrange(self._nxcells):
            for ny in xrange(self._nycells):
                pass


class SWNTGenerator(NanotubeGenerator):
    """Class for generating single-walled nanotubes."""
    pass


class MWNTGenerator(NanotubeGenerator):
    """Class for generating multi-walled nanotubes."""
    pass
