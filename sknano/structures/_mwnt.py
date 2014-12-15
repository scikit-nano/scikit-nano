# -*- coding: utf-8 -*-
"""
==============================================================================
MWNT structure class (:mod:`sknano.structures._mwnt`)
==============================================================================

.. currentmodule:: sknano.structures._mwnt

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
from six.moves import range
__docformat__ = 'restructuredtext en'

import numbers

import numpy as np

from sknano.core.refdata import dVDW

from ._base import StructureBase
from ._swnt import SWNT
from ._compute_funcs import compute_dt

__all__ = ['MWNT']


class MWNT(StructureBase):
    """MWNT structure class."""
    def __init__(self, Ch=None, Nwalls=None, Lz=None, fix_Lz=False, nz=None,
                 add_inner_shells=False, add_outer_shells=True,
                 max_shells=None, max_shell_diameter=np.inf,
                 min_shells=None, min_shell_diameter=0.0,
                 new_shell_type=None, shell_spacing=dVDW, **kwargs):

        super(MWNT, self).__init__(**kwargs)

        self.Ch = Ch
        self.shells = []

        if Ch is None or not isinstance(Ch, list):
            self.Ch = []

            self._add_inner_shells = add_inner_shells
            self._add_outer_shells = add_outer_shells
            self._starting_shell_position = 'outer'

            if max_shells is None:
                max_shells = 10
            self.max_shells = max_shells
            self.max_shell_diameter = max_shell_diameter

            if min_shells is None:
                min_shells = 2
            self.min_shells = min_shells
            self.min_shell_diameter = min_shell_diameter

            self.new_shell_type = new_shell_type
            self.shell_spacing = shell_spacing

            Ch = []
            dt = []
            for n in range(0, 501):
                for m in range(0, 501):
                    if (n <= 2 and m <= 2):
                        continue
                    else:
                        Ch.append((n, m))
                        dt.append(compute_dt(n, m, bond=self.bond))
            Ch = np.asarray(Ch)
            dt = np.asarray(dt)

            self.max_shell_diameter = min(self.max_shell_diameter, dt.max())
            self.min_shell_diameter = max(self.min_shell_diameter, dt.min())

            delta_dt = -2 * self.shell_spacing
            max_dt_diff = 0.05

            if self._add_outer_shells:
                delta_dt = -delta_dt
                self._starting_shell_position = 'inner'

            next_dt = dt[np.where(dt >= self.min_shell_diameter)][0] + delta_dt
            while len(self.Ch) < self.max_shells and \
                    next_dt <= self.max_shell_diameter and \
                    next_dt >= self.min_shell_diameter:

                # get chiral indices for next_dt
                next_Ch_candidates = []
                while len(next_Ch_candidates) == 0 and \
                        next_dt <= self.max_shell_diameter and \
                        next_dt >= self.min_shell_diameter:

                    if self.new_shell_type in ('AC', 'armchair'):
                        next_Ch_candidates = \
                            Ch[np.where(
                                np.logical_and(
                                    np.abs(dt - next_dt) <= max_dt_diff,
                                    Ch[:,0] == Ch[:,1]))]
                    elif self.new_shell_type in ('ZZ', 'zigzag'):
                        next_Ch_candidates = \
                            Ch[np.where(
                                np.logical_and(
                                    np.abs(dt - next_dt) <= max_dt_diff,
                                    np.logical_or(Ch[:,0] == 0,
                                                  Ch[:,1] == 0)))]
                    elif self.new_shell_type == 'achiral':
                        next_Ch_candidates = \
                            Ch[np.where(
                                np.logical_and(
                                    np.abs(dt - next_dt) <= max_dt_diff,
                                    np.logical_or(Ch[:,0] == Ch[:,1],
                                                  np.logical_or(
                                                      Ch[:,0] == 0,
                                                      Ch[:,1] == 0))))]
                    elif self.new_shell_type == 'chiral':
                        next_Ch_candidates = \
                            Ch[np.where(
                                np.logical_and(
                                    np.abs(dt - next_dt) <= max_dt_diff,
                                    np.logical_and(Ch[:,0] != Ch[:,1],
                                                   np.logical_and(
                                                       Ch[:,0] != 0,
                                                       Ch[:,1] != 1))))]
                    else:
                        next_Ch_candidates = \
                            Ch[np.where(np.abs(dt - next_dt) <= max_dt_diff)]

                    if self._add_outer_shells:
                        next_dt += max_dt_diff
                    else:
                        next_dt -= max_dt_diff

                if len(next_Ch_candidates) > 0:
                    n, m = next_Ch_candidates[
                        np.random.choice(np.arange(len(next_Ch_candidates)))]
                    self.Ch.append((n, m))
                    next_dt += delta_dt
                else:
                    break

        for Ch in self.Ch:
            self.shells.append(SWNT(n=Ch[0], m=Ch[-1], Lz=Lz, fix_Lz=fix_Lz,
                                    nz=nz, **kwargs))

        self.fix_Lz = fix_Lz
        if nz is not None:
            self.nz = nz
        else:
            self.nz = 1

        if self.verbose:
            print(self.shells)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        """Return canonical string representation of `MWNT`."""
        strrep = "MWNT(Ch={!r}, Nwalls={!r}, Lz={!r}, fix_Lz={!r}, " + \
            "nz={!r}, element1={!r}, element2={!r}, bond={!r})"
        return strrep.format(self.Ch, self.Nwalls, self.Lz, self.fix_Lz,
                             self.nz, self.element1, self.element2, self.bond)

    @property
    def chiral_types(self):
        """List of chiral types for each `MWNT` shell."""
        return [swnt.chiral_type for swnt in self.shells]

    @property
    def chiral_set(self):
        """Set of all chiral types in `MWNT`."""
        return set(self.chiral_types)

    @property
    def nz(self):
        """Number of nanotube unit cells along the :math:`z`-axis."""
        if len(self.chiral_set) == 1:
            return self._nz
        return [swnt.nz for swnt in self.shells]

    @nz.setter
    def nz(self, value):
        """Set number of nanotube unit cells along the :math:`z`-axis."""
        if not (isinstance(value, numbers.Real) or value > 0):
            raise TypeError('Expected a real, positive number.')
        self._nz = value
        if self._assert_integer_nz:
            self._nz = int(np.ceil(value))

    @property
    def Lz(self):
        """MWNT length :math:`L_z = L_{\\mathrm{tube}}` in **nanometers**."""
        if len(self.chiral_set) == 1:
            return self.nz * self.T
        return np.asarray([swnt.Lz for swnt in self.shells]).max()

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
    def Nwalls(self):
        return len(self.shells)

    @property
    def T(self):
        """Length of `MWNT` unit cell :math:`|\\mathbf{T}|` in \u212b."""
        if len(self.chiral_set) == 1:
            return self.shells[-1].T
        return np.asarray([swnt.T for swnt in self.shells]).max()

    @property
    def dt(self):
        """`MWNT` diameter :math:`d_t=\\frac{|\\mathbf{C}_h|}{\\pi}` \
        in \u212b."""
        return np.asarray([swnt.dt for swnt in self.shells]).max()

    @property
    def rt(self):
        """`MWNT` radius :math:`r_t=\\frac{|\\mathbf{C}_h|}{2\\pi}` \
        in \u212b."""
        return np.asarray([swnt.rt for swnt in self.shells]).max()

    @property
    def Natoms(self):
        """Number of atoms in nanotube.

        .. versionchanged:: 0.3.0

           **Returns total number of atoms per nanotube.**
           Use :attr:`~MWNT.Natoms_per_unit_cell` to get the number of
           atoms per unit cell.

        .. math::

           N_{\\mathrm{atoms}} = 2N\\times n_z =
           \\frac{4(n^2 + m^2 + nm)}{d_R}\\times n_z

        where :math:`N` is the number of graphene hexagons mapped to the
        nanotube unit cell and :math:`n_z` is the number of unit cells.

        """
        return np.asarray([swnt.Natoms for swnt in self.shells]).sum()

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
        return np.asarray([swnt.Natoms_per_unit_cell for swnt in
                           self.shells]).sum()

    @property
    def unit_cell_mass(self):
        """Unit cell mass in atomic mass units."""
        return np.asarray([swnt.unit_cell_mass for swnt in self.shells]).sum()

    @property
    def Ntubes(self):
        """Number of `MWNT`."""
        return 1

    @property
    def Nshells(self):
        return self.Nwalls

    @property
    def tube_mass(self):
        """MWNT mass in **grams**."""
        return np.asarray([swnt.tube_mass for swnt in self.shells]).sum()

    def todict(self):
        return dict(Ch=self.Ch, Nwalls=self.Nwalls, nz=self.nz, Lz=self.Lz,
                    fix_Lz=self.fix_Lz, element1=self.element1,
                    element2=self.element2, bond=self.bond)
