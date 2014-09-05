# -*- coding: utf-8 -*-
"""
==============================================================================
SWNT structure class (:mod:`sknano.structures._swnt`)
==============================================================================

.. currentmodule:: sknano.structures._swnt

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numbers
import numpy as np

from ._base import StructureBase
from ._compute_funcs import *

__all__ = ['SWNT', 'Nanotube']


class SWNT(StructureBase):
    u"""Class for creating interactive SWNT objects.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nz : int, optional
        Number of repeat unit cells along the :math:`z` dimension
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b
    Lz : float, optional
        Length of the nanotube in **nanometers**.
        Overrides the :math:`n_z` cell values.

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

    Examples
    --------

    >>> from sknano.structures import SWNT

    Create a SWNT with :math:`\\mathbf{C}_{h} = (10, 10)` chirality.

    >>> nt = SWNT(n=10, m=10, verbose=True)
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
    def __init__(self, n=10, m=0, nz=1, tube_length=None, Lz=None,
                 fix_Lz=False, **kwargs):

        super(SWNT, self).__init__(**kwargs)

        if tube_length is not None and Lz is None:
            Lz = tube_length

        # add each parameter in the order I want them to appear in
        # verbose output mode
        self._params = ['n', 'm', 't1', 't2', 'd', 'dR', 'N',
                        'R', 'chiral_angle', 'Ch', 'T', 'dt', 'rt',
                        'electronic_type']

        self.n = n
        self.m = m

        if Lz is not None:
            self.nz = 10 * float(Lz) / self.T
        else:
            self.nz = nz

        self._fix_Lz = fix_Lz
        self._assume_integer_unit_cells = True
        if fix_Lz:
            self._assume_integer_unit_cells = False

        if self._assume_integer_unit_cells:
            self.nz = int(np.ceil(self.nz))

        self._L0 = Lz  # store initial value of Lz

    def __str__(self):
        """Return nice string representation of `SWNT`."""
        return repr(self)

    def __repr__(self):
        """Return canonical string representation of `SWNT`."""
        strrep = "SWNT(n={!r}, m={!r}, element1={!r}, element2={!r}, bond={!r}"
        if self._fix_Lz:
            strrep += ", Lz={!r}, fix_Lz={!r})"
            return strrep.format(self.n, self.m, self.element1, self.element2,
                                 self.bond, self.Lz, self._fix_Lz)
        else:
            strrep += ")"
            return strrep.format(self.n, self.m, self.element1,
                                 self.element2, self.bond)

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
        if not (isinstance(value, numbers.Real) or value >= 0):
            raise TypeError('Expected an integer.')
        self._n = int(value)

    @n.deleter
    def n(self):
        del self._n

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
        if not (isinstance(value, numbers.Real) or value >= 0):
            raise TypeError('Expected an integer.')
        self._m = int(value)

    @m.deleter
    def m(self):
        del self._m

    @property
    def d(self):
        """:math:`d=\\gcd{(n, m)}`

        :math:`d` is the **G**\ reatest **C**\ ommon **D**\ ivisor of
        :math:`n` and :math:`m`.

        """
        return compute_d(self.n, self.m)

    @property
    def dR(self):
        """:math:`d_R=\\gcd{(2n + m, 2m + n)}`

        :math:`d_R` is the **G**\ reatest **C**\ ommon **D**\ ivisor of
        :math:`2n + m` and :math:`2m + n`.

        """
        return compute_dR(self.n, self.m)

    @property
    def N(self):
        """Number of graphene hexagons in nanotube *unit cell*.

        .. math::

           N = \\frac{4(n^2 + m^2 + nm)}{d_R}

        """
        return compute_N(self.n, self.m)

    @property
    def t1(self):
        """:math:`t_{1} = \\frac{2m + n}{d_{R}}`

        where :math:`d_R = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_1`:

        .. math::

           \\mathbf{T} = t_1\\mathbf{a}_{1} + t_2\\mathbf{a}_2

        """
        return compute_t1(self.n, self.m)

    @property
    def t2(self):
        """:math:`t_2 = -\\frac{2n + m}{d_R}`

        where :math:`d_R = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_2`:

        .. math::

           \\mathbf{T} = t_1\\mathbf{a}_1 + t_2\\mathbf{a}_2

        """
        return compute_t2(self.n, self.m)

    @property
    def Ch(self):
        u"""SWNT circumference :math:`|\\mathbf{C}_h|` in **\u212b**"""
        return compute_Ch(self.n, self.m, bond=self.bond)

    @property
    def chiral_angle(self):
        """Chiral angle :math:`\\theta_c` in **degrees**.

        .. math::

           \\theta_c = \\tan^{-1}\\left(\\frac{\\sqrt{3} m}{2n + m}\\right)

        """
        return compute_chiral_angle(self.n, self.m)

    @property
    def T(self):
        u"""Length of nanotube unit cell :math:`|\\mathbf{T}|` in \u212b.

        .. math::

           |\\mathbf{T}| = \\frac{\\sqrt{3} |\\mathbf{C}_{h}|}{d_{R}}

        """
        return compute_T(self.n, self.m, bond=self.bond, length=True)

    @property
    def dt(self):
        u"""Nanotube diameter :math:`d_t = \\frac{|\\mathbf{C}_h|}{\\pi}`
        in \u212b."""
        return compute_dt(self.n, self.m, bond=self.bond)

    @property
    def rt(self):
        u"""Nanotube radius :math:`r_t = \\frac{|\\mathbf{C}_h|}{2\\pi}`
        in \u212b."""
        return compute_rt(self.n, self.m, bond=self.bond)

    @property
    def M(self):
        """:math:`M = np - nq`

        :math:`M` is the number of multiples of the translation vector
        :math:`\\mathbf{T}` in the vector :math:`N\\mathbf{R}`.

        """
        return compute_M(self.n, self.m)

    @property
    def R(self):
        """Symmetry vector :math:`\\mathbf{R} = (p, q)`.

        .. math::

           \\mathbf{R} = p\\mathbf{a}_1 + q\\mathbf{a}_2

        """
        return compute_R(self.n, self.m, bond=self.bond, length=False)

    @property
    def nz(self):
        """Number of nanotube unit cells along the :math:`z`-axis."""
        return self._nz

    @nz.setter
    def nz(self, value):
        """Set number of nanotube unit cells along the :math:`z`-axis."""
        if not (isinstance(value, numbers.Real) or value > 0):
            raise TypeError('Expected a real, positive number.')
        if self._assume_integer_unit_cells:
            self._nz = int(value)
        else:
            self._nz = value

    @nz.deleter
    def nz(self):
        del self._nz

    @property
    def electronic_type(self):
        """SWNT electronic type.

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
        return compute_electronic_type(self.n, self.m)

    @property
    def Natoms(self):
        """Number of atoms in nanotube.

        .. versionchanged:: 0.3.0

           **Returns total number of atoms per nanotube.**
           Use :attr:`~SWNT.Natoms_per_unit_cell` to get the number of
           atoms per unit cell.

        .. math::

           N_{\\mathrm{atoms}} = 2N\\times n_z =
           \\frac{4(n^2 + m^2 + nm)}{d_R}\\times n_z

        where :math:`N` is the number of graphene hexagons mapped to the
        nanotube unit cell and :math:`n_z` is the number of unit cells.

        """
        return compute_Natoms(self.n, self.m, self.nz)

    @property
    def Natoms_per_tube(self):
        """Number of atoms in nanotube :math:`N_{\\mathrm{atoms/tube}}`."""
        return self.Natoms

    @property
    def Natoms_per_unit_cell(self):
        """Number of atoms in nanotube unit cell.

        .. math::

           N_{\\mathrm{atoms}} = 2N = \\frac{4(n^2 + m^2 + nm)}{d_R}

        where :math:`N` is the number of graphene hexagons mapped to the
        nanotube unit cell.

        """
        return compute_Natoms_per_unit_cell(self.n, self.m)

    @property
    def unit_cell_mass(self):
        """Unit cell mass in atomic mass units."""
        return compute_unit_cell_mass(self.n, self.m,
                                      element1=self.element1,
                                      element2=self.element2)

    @property
    def linear_mass_density(self):
        """Linear mass density of nanotube in g/nm."""
        return compute_linear_mass_density(self.n, self.m, bond=self.bond,
                                           element1=self.element1,
                                           element2=self.element2)

    @property
    def Ntubes(self):
        """Number of nanotubes."""
        return 1

    @property
    def Lz(self):
        """SWNT length :math:`L_z = L_{\\mathrm{tube}}` in
        **nanometers**."""
        return self.nz * self.T

    @property
    def tube_length(self):
        """Alias for :attr:`SWNT.Lz`"""
        return self.Lz

    @property
    def tube_mass(self):
        """SWNT mass in **grams**."""
        return compute_tube_mass(self.n, self.m, self.nz,
                                 element1=self.element1,
                                 element2=self.element2)

Nanotube = SWNT
