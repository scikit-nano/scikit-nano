# -*- coding: utf-8 -*-
"""
==============================================================================
Mixin structure classes (:mod:`sknano.structures._mixins`)
==============================================================================

.. currentmodule:: sknano.structures._mixins

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
from six.moves import range
__docformat__ = 'restructuredtext en'

import numbers

import numpy as np

from sknano.core.math import Vector
from ._compute_funcs import *
from ._extras import get_Ch_type

__all__ = ['MWNTMixin', 'NanotubeMixin', 'NanotubeBundleMixin',
           'UnrolledSWNTMixin']


class NanotubeMixin(object):
    """Mixin class for nanotube classes."""
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
        """SWNT circumference :math:`|\\mathbf{C}_h|` in **\u212b**"""
        return compute_Ch(self.n, self.m, bond=self.bond)

    @property
    def dt(self):
        """Nanotube diameter :math:`d_t = \\frac{|\\mathbf{C}_h|}{\\pi}` \
        in \u212b."""
        return compute_dt(self.n, self.m, bond=self.bond)

    @property
    def rt(self):
        """Nanotube radius :math:`r_t = \\frac{|\\mathbf{C}_h|}{2\\pi}` \
        in \u212b."""
        return compute_rt(self.n, self.m, bond=self.bond)

    @property
    def chiral_angle(self):
        """Chiral angle :math:`\\theta_c` in **degrees**.

        .. math::

           \\theta_c = \\tan^{-1}\\left(\\frac{\\sqrt{3} m}{2n + m}\\right)

        """
        return compute_chiral_angle(self.n, self.m)

    @property
    def chiral_type(self):
        return get_Ch_type((self.n, self.m))

    @property
    def T(self):
        """Length of nanotube unit cell :math:`|\\mathbf{T}|` in \u212b.

        .. math::

           |\\mathbf{T}| = \\frac{\\sqrt{3} |\\mathbf{C}_{h}|}{d_{R}}

        """
        return compute_T(self.n, self.m, bond=self.bond, length=True)

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
        self._nz = value
        if self._assert_integer_nz:
            self._nz = int(np.ceil(value))

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
    def Lz(self):
        """SWNT length :math:`L_z = L_{\\mathrm{tube}}` in
        **nanometers**."""
        return self.nz * self.T

    @property
    def fix_Lz(self):
        return self._fix_Lz

    @fix_Lz.setter
    def fix_Lz(self, value):
        if not isinstance(value, bool):
            raise TypeError('Expected `True` or `False`')
        self._fix_Lz = value
        self._assert_integer_nz = True
        if self.fix_Lz:
            self._assert_integer_nz = False

    @property
    def unit_cell_symmetry_params(self):
        psi, tau = compute_symmetry_operation(self.n, self.m, bond=self.bond)
        aCh = compute_chiral_angle(self.n, self.m, rad2deg=False)
        dpsi = self.bond * np.cos(np.pi / 6 - aCh) / self.rt
        dtau = self.bond * np.sin(np.pi / 6 - aCh)

        return psi, tau, dpsi, dtau


class MWNTMixin(object):
    """Mixin class for MWNTs."""
    pass


class NanotubeBundleMixin(object):
    """Mixin class for nanotube bundles."""

    def compute_bundle_params(self):
        """Compute/update nanotube bundle parameters."""

        self.r1.x = self.dt + self.vdw_spacing
        if self.bundle_packing is None and \
                self.bundle_geometry in ('square', 'rectangle'):
            self._bundle_packing = 'ccp'
        elif self.bundle_packing is None and \
                self.bundle_geometry in ('triangle', 'hexagon'):
            self._bundle_packing = 'hcp'

        if self.bundle_packing in ('cubic', 'ccp'):
            self.r2.y = self.r1.x
        else:
            self.r2.x = self.r1.x * np.cos(2 * np.pi / 3)
            self.r2.y = self.r1.x * np.sin(2 * np.pi / 3)
            if self.bundle_packing is None:
                self._bundle_packing = 'hcp'

        if self.bundle_geometry == 'hexagon':
            nrows = max(self.nx, self.ny, 3)
            if nrows % 2 != 1:
                nrows += 1

            ntubes_per_end_rows = int((nrows + 1) / 2)

            row = 0
            ntubes_per_row = nrows
            while ntubes_per_row >= ntubes_per_end_rows:
                if row == 0:
                    for n in range(ntubes_per_row):
                        dr = n * self.r1
                        self.bundle_coords.append(dr)
                else:
                    for nx in range(ntubes_per_row):
                        for ny in (-row, row):
                            dr = Vector()
                            dr.x = abs(ny * self.r2.x)
                            dr.y = ny * self.r2.y
                            dr = nx * self.r1 + dr
                            self.bundle_coords.append(dr)
                row += 1
                ntubes_per_row = nrows - row

        elif self.bundle_geometry == 'rectangle':
            Lx = 10 * self.Lx
            for nx in range(self.nx):
                for ny in range(self.ny):
                    dr = nx * self.r1 + ny * self.r2
                    while dr.x < 0:
                        dr.x += Lx
                    self.bundle_coords.append(dr)

        elif self.bundle_geometry == 'square':
            pass
        elif self.bundle_geometry == 'triangle':
            pass
        else:
            for nx in range(self.nx):
                for ny in range(self.ny):
                    dr = nx * self.r1 + ny * self.r2
                    self.bundle_coords.append(dr)

    @property
    def nx(self):
        """Number of nanotubes along the :math:`x`-axis."""
        return self._nx

    @nx.setter
    def nx(self, value):
        """Set :math:`n_x`"""
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a positive integer.')
        self._nx = int(value)

    @nx.deleter
    def nx(self):
        del self._nx

    @property
    def ny(self):
        """Number of nanotubes along the :math:`y`-axis."""
        return self._ny

    @ny.setter
    def ny(self, value):
        """Set :math:`n_y`"""
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a positive integer.')
        self._ny = int(value)

    @ny.deleter
    def ny(self):
        del self._ny

    @property
    def Lx(self):
        return self.nx * (self.dt + self.vdw_spacing) / 10

    @property
    def Ly(self):
        return self.ny * (self.dt + self.vdw_spacing) / 10

    @property
    def bundle_packing(self):
        return self._bundle_packing

    @bundle_packing.setter
    def bundle_packing(self, value):
        if value not in ('ccp', 'hcp'):
            raise ValueError('Expected value to be `hcp` or `ccp`')
        self._bundle_packing = value
        self.compute_bundle_params()

    @bundle_packing.deleter
    def bundle_packing(self):
        del self._bundle_packing

    @property
    def bundle_mass(self):
        return self.Ntubes * self.tube_mass

    @property
    def Natoms_per_bundle(self):
        return self.Ntubes * self.Natoms_per_tube

    @property
    def Ntubes(self):
        return len(self.bundle_coords)


class UnrolledSWNTMixin(object):
    """Mixin class for unrolled nanotubes."""

    @property
    def Lx(self):
        return self.nx * self.Ch / 10

    @property
    def Ly(self):
        return self.Nlayers * self.layer_spacing / 10

    @property
    def nx(self):
        """Number of unit cells along the :math:`x`-axis."""
        return self._nx

    @nx.setter
    def nx(self, value):
        """Set :math:`n_x`"""
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a real positive number.')
        self._nx = int(value)

    @property
    def Nlayers(self):
        """Number of layers."""
        return self._Nlayers

    @Nlayers.setter
    def Nlayers(self, value):
        """Set :attr:`Nlayers`."""
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a real positive number.')
        self._Nlayers = int(value)

    @Nlayers.deleter
    def Nlayers(self):
        del self._Nlayers

    @property
    def fix_Lx(self):
        return self._fix_Lx

    @fix_Lx.setter
    def fix_Lx(self, value):
        if not isinstance(value, bool):
            raise TypeError('Expected `True` or `False`')
        self._fix_Lx = value
        self._assert_integer_nx = True
        if self.fix_Lx:
            self._assert_integer_nx = False
