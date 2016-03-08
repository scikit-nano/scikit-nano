# -*- coding: utf-8 -*-
"""
==============================================================================
Nanotube bundle base class (:mod:`sknano.core.structures._nanotube_bundle`)
==============================================================================

.. currentmodule:: sknano.core.structures._nanotube_bundle

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numbers
import re

import numpy as np

from sknano.core.atoms import Atom, vdw_radius_from_basis
from sknano.core.refdata import aCC, grams_per_Da
from sknano.core.math import Vector
from ._extras import get_chiral_indices

__all__ = ['compute_bundle_density', 'NanotubeBundleMixin',
           'NanotubeBase']


def compute_bundle_density(*Ch, r_vdw=None, bond=None,
                           element1=None, element2=None):
    """Compute nanotube bundle mass density \
    :math:`\\rho_{\\mathrm{bundle}}(n, m)` in :math:`\\mathrm{g/cm^3}`.

    .. math::

        \\rho_{\\mathrm{bundle}}(n, m) = \\frac{8\\pi^2 m_{\\mathrm{C}}
        \\sqrt{n^2 + m^2 + nm}}{9\\sqrt{3}a_{\\mathrm{CC}}^3 \\times
        \\left(\\sqrt{n^2 + m^2 + nm} +
        \\frac{\\pi d_{\\mathrm{vdW}}}{\\sqrt{3}a_{\\mathrm{CC}}}\\right)^2}

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    r_vdw : int
        van der Waals radius of nanotube atoms
    bond : float, optional
        Bond length.

    Returns
    -------
    float
        :math:`\\rho_{\\mathrm{bundle}}` in units of
        :math:`\\mathrm{\\frac{g}{cm^3}}`

    """
    n, m, _ = get_chiral_indices(*Ch)

    if bond is None:
        bond = aCC

    if element1 is None:
        element1 = 'C'
    if element2 is None:
        element2 = 'C'

    if r_vdw is None:
        r_vdw = vdw_radius_from_basis(element1, element2)

    if element1 == element2:
        bundle_density = 8 * np.pi ** 2 * Atom(element1).mass * \
            np.sqrt(n ** 2 + m ** 2 + n * m) / \
            (9 * np.sqrt(3) * bond ** 3 *
                (np.sqrt(n ** 2 + m ** 2 + n * m) +
                    2 * np.pi * r_vdw / (np.sqrt(3) * bond)) ** 2)
    else:
        bundle_density = 0

    # there are 1.6605e-24 grams / Da and 1e-8 cm / angstrom
    bundle_density *= grams_per_Da / (1e-8) ** 3
    return bundle_density


class NanotubeBundleMixin:
    """Mixin class for nanotube bundles."""
    @property
    def bundle_density(self):
        """Compute the nanotube bundle density."""
        return compute_bundle_density(self.n, self.m, r_vdw=self.vdw_radius,
                                      bond=self.bond, basis=self.basis)

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
        """Axis-aligned length along the `x`-axis in **Angstroms**.

        Calculated as:

        .. math::

           L_x = n_x * (d_t + 2 r_{\\mathrm{vdW}})

        """
        return self.nx * (self.dt + 2 * self.vdw_radius)

    @property
    def Ly(self):
        """Axis-aligned length along the `y`-axis in **Angstroms**.

        Calculated as:

        .. math::

           L_y = n_y * (d_t + 2 r_{\\mathrm{vdW}})

        """
        return self.ny * (self.dt + 2 * self.vdw_radius)

    @property
    def bundle_geometry(self):
        """Bundle geometry."""
        return self._bundle_geometry

    @bundle_geometry.setter
    def bundle_geometry(self, value):
        if value is not None and value not in self._bundle_geometries:
            print('Unrecognized `bundle_geometry`: {!r}'.format(value))
            value = None
        self._bundle_geometry = value

    @property
    def bundle_packing(self):
        """Bundle packing."""
        return self._bundle_packing

    @bundle_packing.setter
    def bundle_packing(self, value):
        if value is None and \
                self.bundle_geometry in ('square', 'rectangle'):
            value = 'ccp'
        elif value is None and \
                self.bundle_geometry in ('triangle', 'hexagon'):
            value = 'hcp'

        if value is not None and value not in ('ccp', 'hcp'):
            raise ValueError('Expected value to be `hcp` or `ccp`')

        self._bundle_packing = value
        # self.generate_bundle_coords()

    @bundle_packing.deleter
    def bundle_packing(self):
        del self._bundle_packing

    @property
    def bundle_mass(self):
        """Bundle mass."""
        return self.Ntubes * self.tube_mass

    @property
    def Natoms(self):
        """Number of atoms in nanotube bundle.

           **Returns total number of atoms in nanotube bundle.**
           Use :attr:`~NanotubeBundleMixin.Natoms_per_tube` to
           get a list of the number of atoms in each nanotube in
           the bundle.

        """
        if self.is_bundle:
            return sum(self.Natoms_list)
        return super().Natoms

    @property
    def Natoms_per_bundle(self):
        """Alias for :attr:`~NanotubeBundleMixin.Natoms`."""
        return self.Natoms

    @property
    def Natoms_list(self):
        """:class:`~python:list` of `Natoms` per nanotube in bundle."""
        if self.is_bundle:
            return [nanotube.Natoms for nanotube in self.bundle_list]
        return super().Natoms_list

    @property
    def Ntubes(self):
        """Number of nanotubes in bundle."""
        try:
            return len(self.bundle_coords)
        except AttributeError:
            return 1

    @property
    def Natoms_per_tube(self):
        """Alias for :attr:`~NanotubeBundleMixin.Natoms_list`."""
        try:
            val = self.Natoms_list[:]
            return val if len(val) > 1 else val[0]
        except AttributeError:
            return super().Natoms_per_tube

    def init_bundle_parameters(self):
        """Initialize bundle attributes."""
        self.bundle_list = []
        self.generate_bundle_coords()
        fmtstr = super().fmtstr
        match = re.search('[nL]z=', fmtstr)
        if match:
            self.fmtstr = fmtstr[:match.start()] + \
                "nx={nx!r}, ny={ny!r}, " + fmtstr[match.start():] + \
                ", bundle_packing={bundle_packing!r}, " + \
                "bundle_geometry={bundle_geometry!r}"

    def generate_bundle_coords(self):
        """Generate coordinates of bundle tubes."""
        self.r1 = Vector()
        self.r2 = Vector()
        self.bundle_coords = []

        self.r1.x = self.dt + 2 * self.vdw_radius
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
            Lx = self.Lx
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


class NanotubeBase(NanotubeBundleMixin):
    """Nanotube bundle structure base class."""

    _bundle_geometries = ['square', 'rectangle', 'hexagon']

    def __init__(self, *args, basis=['C', 'C'], bond=aCC, nx=1, ny=1,
                 bundle_packing=None, bundle_geometry=None, **kwargs):

        super().__init__(*args, basis=basis, bond=bond, **kwargs)

        self.nx = nx
        self.ny = ny
        self.bundle_geometry = bundle_geometry
        self.bundle_packing = bundle_packing

        self.is_bundle = False
        if nx != 1 or ny != 1 or bundle_geometry is not None:
            self.is_bundle = True
            self.init_bundle_parameters()

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attrdict = super().todict()
        if self.is_bundle:
            attrdict.update(dict(nx=self.nx, ny=self.ny,
                                 bundle_packing=self.bundle_packing,
                                 bundle_geometry=self.bundle_geometry))
        return attrdict
