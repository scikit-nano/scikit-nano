# -*- coding: utf-8 -*-
"""
==============================================================================
Mixin structure classes (:mod:`sknano.structures._mixins`)
==============================================================================

.. currentmodule:: sknano.structures._mixins

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core.atoms import Atom
from sknano.core.refdata import CCbond, dVDW, grams_per_Da
from ._swnt import Nanotube

__all__ = ['NanotubeBundleMixin', 'MWNTMixin', 'UnrolledSWNTMixin']


class UnrolledSWNTMixin(object):
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

    @classmethod
    def compute_Lx(cls, n=int, m=int, nx=int, bond=None, **kwargs):
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
        dt = Nanotube.compute_Ch(n=n, m=m, bond=bond, **kwargs)
        return nx * dt / 10

    @classmethod
    def compute_Ly(cls, n=int, m=int, ny=int, bond=None, **kwargs):
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
        return ny * dVDW / 10


class MWNTMixin(object):

    @property
    def max_shells(self):
        return self._max_shells

    @max_shells.setter
    def max_shells(self, value):
        self._max_shells = value

    @property
    def min_shells(self):
        return self._min_shells

    @min_shells.setter
    def min_shells(self, value):
        self._min_shells = value

    @property
    def Natoms_per_tube(self):
        """Number of atoms in nanotube :math:`N_{\\mathrm{atoms/tube}}`."""
        return self._Natoms_per_tube

    @property
    def Nshells_per_tube(self):
        """Number of shells in MWNT :math:`N_{\\mathrm{shells}}`."""
        return self._Nshells_per_tube


class NanotubeBundleMixin(object):

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

    @property
    def r1(self):
        """Bundle lattice vector :math:`\\mathbf{r}_1`."""
        return self._r1

    @property
    def r2(self):
        """Bundle lattice vector :math:`\\mathbf{r}_2`."""
        return self._r2

    @property
    def bundle_geometry(self):
        """Bundle geometry."""
        return self._bundle_geometry

    @property
    def bundle_packing(self):
        """Bundle packing arrangement."""
        return self._bundle_packing

    @property
    def bundle_density(self):
        """Nanotube bundle mass density :math:`\\rho_{\\mathrm{bundle}}` in
        :math:`\\mathrm{g/cm^3}`."""
        return self._bundle_density

    @property
    def bundle_mass(self):
        """Nanotube bundle mass :math:`M_{\\mathrm{bundle}}` in **grams**."""
        return self._bundle_mass

    @classmethod
    def compute_Lx(cls, n=int, m=int, nx=int, bond=None, **kwargs):
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
        dt = Nanotube.compute_dt(n=n, m=m, bond=bond, **kwargs)
        return nx * (dt + dVDW) / 10

    @classmethod
    def compute_Ly(cls, n=int, m=int, ny=int, bond=None, **kwargs):
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
        dt = Nanotube.compute_dt(n=n, m=m, bond=bond, **kwargs)
        return ny * (dt + dVDW) / 10

    def compute_bundle_params(self, d_vdw=None):
        """Compute/update nanotube bundle parameters."""
        print('In compute_bundle_params\n' +
              'called by type(self): {}'.format(type(self)))

        self._r1.x = \
            Nanotube.compute_dt(n=self.n, m=self.m, bond=self.bond) + \
            self._vdw_spacing
        if self._bundle_packing is None and \
                self._bundle_geometry in ('square', 'rectangle'):
            self._bundle_packing = 'ccp'
        elif self._bundle_packing is None and \
                self._bundle_geometry in ('triangle', 'hexagon'):
            self._bundle_packing = 'hcp'

        if self._bundle_packing in ('cubic', 'ccp'):
            self._r2.y = self._r1.x
        else:
            self._r2.x = self._r1.x * np.cos(2 * np.pi / 3)
            self._r2.y = self._r1.x * np.sin(2 * np.pi / 3)
            if self._bundle_packing is None:
                self._bundle_packing = 'hcp'

        self._Lx = \
            self.compute_Lx(n=self._n, m=self._m, nx=self._nx, bond=self._bond)

        self._Ly = \
            self.compute_Ly(n=self._n, m=self._m, ny=self._ny, bond=self._bond)

        self._Ntubes = self.compute_Ntubes(nx=self._nx, ny=self._ny)
        self._bundle_mass = \
            self.compute_bundle_mass(n=self._n, m=self._m,
                                     nx=self._nx, ny=self._ny, nz=self._nz)
        self._Natoms_per_bundle = self._Ntubes * self._Natoms_per_tube
        if d_vdw is None:
            d_vdw = self._vdw_spacing

        self._bundle_density = \
            self.compute_bundle_density(n=self.n, m=self.m, d_vdw=d_vdw,
                                        bond=self.bond, element1=self.element1,
                                        element2=self.element2)

    @classmethod
    def compute_bundle_density(cls, n=int, m=int, d_vdw=None, bond=None,
                               element1=None, element2=None):
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

        if d_vdw is None:
            if n == m:
                d_vdw = 3.38
            elif (m == 0) or (n == 0):
                d_vdw = 3.41
            else:
                d_vdw = 3.39

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

        # there are 1.6605e-24 grams / Da and 1e-8 cm / angstrom
        bundle_density *= grams_per_Da / (1e-8)**3
        return bundle_density

    @classmethod
    def compute_bundle_mass(cls, n=int, m=int, nx=int, ny=int, nz=float,
                            element1=None, element2=None):
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
            :class:`~sknano.core.atoms.Atoms` 1 and 2

        Returns
        -------
        float
            Nanotube bundle mass :math:`M_{\\mathrm{bundle}}` in **grams**.

        """
        Ntubes = int(nx * ny)
        tube_mass = Nanotube.compute_tube_mass(n=n, m=m, nz=nz,
                                               element1=element1,
                                               element2=element2)
        return Ntubes * tube_mass

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
