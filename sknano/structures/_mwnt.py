# -*- coding: utf-8 -*-
"""
==============================================================================
MWNT structure class (:mod:`sknano.structures._mwnt`)
==============================================================================

.. currentmodule:: sknano.structures._mwnt

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numbers
import numpy as np

from sknano.core.refdata import dVDW

from ._base import StructureBase
from ._swnt import SWNT
from ._compute_funcs import compute_dt

__all__ = ['MWNT']


class MWNT(StructureBase):
    def __init__(self, Ch=None, Nwalls=None, Lz=None, fix_Lz=False,
                 add_inner_shells=False, add_outer_shells=True,
                 max_shells=None, max_shell_diameter=np.inf,
                 min_shells=None, min_shell_diameter=0.0,
                 new_shell_type=None, shell_spacing=dVDW, **kwargs):

        super(MWNT, self).__init__(**kwargs)

        self.shells = []

        self.Ch = Ch
        if Ch is None or not isinstance(Ch, list):
            self.Ch = []
            self._L0 = Lz  # store initial value of Lz

            self._fix_Lz = fix_Lz
            self._assume_integer_unit_cells = True
            if fix_Lz:
                self._assume_integer_unit_cells = False

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
            for n in xrange(0, 501):
                for m in xrange(0, 501):
                    if (n <= 2 and m <= 2):
                        continue
                    else:
                        Ch.append((n, m))
                        dt.append(compute_dt(n, m, bond=self.bond))
            Ch = np.asarray(Ch)
            dt = np.asarray(dt)

            self.max_shell_diameter = min(self.max_shell_diameter, dt.max())
            self.min_shell_diameter = max(self.min_shell_diameter, dt.min())

            #Lzmin = self.Lz

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
            print(self.shells)
            self.shells.append(SWNT(Ch[0], Ch[-1]))

    def __repr__(self):
        """Return canonical string representation of `SWNT`."""
        retstr = "MWNT(Ch={!r}, ".format(self.Ch) + \
            "element1={!r}, element2={!r}, bond={!r}".format(
                self.element1, self.element2, self.bond)
        if self._fix_Lz:
            retstr += ", Lz={!r}, fix_Lz={!r})".format(self.Lz, self._fix_Lz)
        else:
            retstr += ")"

        return retstr

    @property
    def Nwalls(self):
        return len(self.shells)

    @property
    def dt(self):
        u"""`MWNT` diameter :math:`d_t = \\frac{|\\mathbf{C}_h|}{\\pi}`
        in \u212b."""
        return self.shells[-1].dt

    @property
    def rt(self):
        u"""`MWNT` radius :math:`r_t = \\frac{|\\mathbf{C}_h|}{2\\pi}`
        in \u212b."""
        return self.shells[-1].rt

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
        """SWNT mass in **grams**."""
        return np.asarray([swnt.tube_mass for swnt in self.shells]).sum()

    @property
    def max_shells(self):
        return self._max_shells

    @max_shells.setter
    def max_shells(self, value):
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a real, positive integer.')
        self._max_shells = int(value)

    @property
    def min_shells(self):
        return self._min_shells

    @min_shells.setter
    def min_shells(self, value):
        if not (isinstance(value, numbers.Number) or value >= 0):
            raise TypeError('Expected a real integer.')
        self._min_shells = int(value)
