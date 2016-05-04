# -*- coding: utf-8 -*-
"""
==============================================================================
Nanotube bundle classes (:mod:`sknano.core.structures.nanotube_bundle`)
==============================================================================

.. currentmodule:: sknano.core.structures.nanotube_bundle

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numbers
import re

import numpy as np

from sknano.core import deprecated, deprecated_kwargs
from sknano.core.atoms import Atom, vdw_radius_from_basis
from sknano.core.refdata import aCC, grams_per_Da
from sknano.core.math import Vector
from .base import NanoStructureBase
from .extras import get_chiral_indices

__all__ = ['compute_bundle_density', 'NanotubeBundleMixin',
           'NanotubeBundleBase']


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
    @deprecated(since='0.4.0', alternative='n1', obj_type='attribute')
    def nx(self):
        """Number of nanotubes along the :math:`x`-axis."""
        return self.n1

    @nx.setter
    def nx(self, value):
        self.n1 = value

    @property
    @deprecated(since='0.4.0', alternative='n2', obj_type='attribute')
    def ny(self):
        """Number of nanotubes along the :math:`x`-axis."""
        return self.n2

    @ny.setter
    def ny(self, value):
        self.n2 = value

    @property
    @deprecated(since='0.4.0', alternative='lattice.a', obj_type='attribute')
    def Lx(self):
        """Axis-aligned length along the `x`-axis in **Angstroms**."""
        return self.lattice.a

    @property
    @deprecated(since='0.4.0', alternative='lattice.b', obj_type='attribute')
    def Ly(self):
        """Axis-aligned length along the `y`-axis in **Angstroms**."""
        return self.lattice.b

    @property
    def n1(self):
        """Number of nanotubes along the :math:`x`-axis."""
        return self._n1

    @n1.setter
    def n1(self, value):
        """Set :math:`n_x`"""
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a positive integer.')
        self._n1 = int(value)

    @property
    def n2(self):
        """Number of nanotubes along the :math:`y`-axis."""
        return self._n2

    @n2.setter
    def n2(self, value):
        """Set :math:`n_y`"""
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a positive integer.')
        self._n2 = int(value)

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
    def mass(self):
        """Bundle mass."""
        return self.Ntubes * super().mass

    @property
    def bundle_mass(self):
        """An alias for :attr:`~NanotubeBundleMixin.mass`."""
        return self.mass

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
        # try:
        #     val = self.Natoms_list[:]
        #     return val if len(val) > 1 else val[0]
        # except AttributeError:
        #     return super().Natoms_per_tube
        return super().Natoms_per_tube

    def init_bundle_parameters(self):
        """Initialize bundle attributes."""
        self.bundle_list = []
        self.generate_bundle_coords()
        fmtstr = super().fmtstr
        match = re.search('(n3|L)=', fmtstr)
        if match:
            self.fmtstr = fmtstr[:match.start()] + \
                "n1={n1!r}, n2={n2!r}, " + fmtstr[match.start():] + \
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
            nrows = max(self.n1, self.n2, 3)
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
                    for n1 in range(ntubes_per_row):
                        for n2 in (-row, row):
                            dr = Vector()
                            dr.x = abs(n2 * self.r2.x)
                            dr.y = n2 * self.r2.y
                            dr = n1 * self.r1 + dr
                            self.bundle_coords.append(dr)
                row += 1
                ntubes_per_row = nrows - row

        elif self.bundle_geometry == 'rectangle':
            a = self.lattice.a
            for n1 in range(self.n1):
                for n2 in range(self.n2):
                    dr = n1 * self.r1 + n2 * self.r2
                    while dr.x < 0:
                        dr.x += a
                    self.bundle_coords.append(dr)

        elif self.bundle_geometry == 'square':
            pass
        elif self.bundle_geometry == 'triangle':
            pass
        else:
            for n1 in range(self.n1):
                for n2 in range(self.n2):
                    dr = n1 * self.r1 + n2 * self.r2
                    self.bundle_coords.append(dr)


class NanotubeBundleBase(NanotubeBundleMixin, NanoStructureBase):
    """Nanotube bundle structure base class.

    Parameters
    ----------
    n1, n2 : :class:`~python:int`
    bundle_packing : {None, 'hcp', 'ccp'}, optional
    bundle_geometry : {None, 'hexagon', 'rectangle'}, optional

    """

    _bundle_geometries = ['square', 'rectangle', 'hexagon']

    @deprecated_kwargs({'nx': 'n1', 'ny': 'n2'}, since='0.4.0')
    def __init__(self, *args, n1=1, n2=1, bundle_packing=None,
                 bundle_geometry=None, **kwargs):
        super().__init__(*args, **kwargs)

        [setattr(atom, 'mol', 1) for atom in self.crystal_cell.basis]

        self.n1 = n1
        self.n2 = n2
        self.bundle_geometry = bundle_geometry
        self.bundle_packing = bundle_packing

        self.is_bundle = False
        if n1 != 1 or n2 != 1 or bundle_geometry is not None:
            self.is_bundle = True
            self.init_bundle_parameters()
            self.__generate_unit_cell()

    def generate_unit_cell(self):
        """Generate Nanotube unit cell."""
        self.scaling_matrix = [self.n1, self.n2, 1]

    __generate_unit_cell = generate_unit_cell

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attrdict = super().todict()
        if self.is_bundle:
            attrdict.update(dict(n1=self.n1, n2=self.n2,
                                 bundle_packing=self.bundle_packing,
                                 bundle_geometry=self.bundle_geometry))
        return attrdict
