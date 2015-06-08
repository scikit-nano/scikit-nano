# -*- coding: utf-8 -*-
"""
==============================================================================
MWNT structure class (:mod:`sknano.structures._mwnt`)
==============================================================================

.. currentmodule:: sknano.structures._mwnt

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
from builtins import range
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core.crystallography import Crystal3DLattice, UnitCell
from sknano.core.refdata import aCC, dVDW

from ._base import StructureBase
from ._swnt import SWNT
from ._compute_funcs import compute_dt, compute_T
from ._extras import generate_Ch_list
from ._mixins import MWNTMixin

__all__ = ['MWNT']


class MWNT(MWNTMixin, StructureBase):
    """MWNT structure class."""
    def __init__(self, Ch_list=None, Nwalls=None, Lz=None,
                 min_shell_diameter=None, max_shell_diameter=None,
                 max_shells=None, chiral_types=None, shell_spacing=dVDW,
                 basis=['C', 'C'], bond=aCC, **kwargs):
        if Ch_list is None and 'Ch' in kwargs:
            Ch_list = kwargs['Ch']
            del kwargs['Ch']

        super().__init__(bond=bond, **kwargs)

        if Ch_list is None or not isinstance(Ch_list, list):

            if Nwalls is not None:
                max_shells = Nwalls

            if max_shells is None:
                max_shells = 10

            if max_shell_diameter is None:
                max_shell_diameter = np.inf

            if min_shell_diameter is None:
                min_shell_diameter = 5.0

            delta_dt = 2 * shell_spacing

            imax = 100

            self._Ch_pool = \
                np.asarray(generate_Ch_list(imax=imax,
                                            chiral_type=chiral_types))
            self._dt_pool = np.asarray([compute_dt(_Ch, bond=self.bond) for _Ch
                                       in self._Ch_pool])

            dt_mask = np.logical_and(self._dt_pool >= min_shell_diameter,
                                     self._dt_pool <= max_shell_diameter)

            self._Ch_pool = self._Ch_pool[dt_mask]
            self._dt_pool = self._dt_pool[dt_mask]

            if max_shell_diameter < np.inf:
                dt_list = []
                dt = self._dt_pool.min()
                while dt <= max_shell_diameter:
                    dt_list.append(dt)
                    dt += delta_dt
            else:
                dt_list = [self._dt_pool.min() + i * delta_dt
                           for i in range(max_shells)]

            dt_masks = [self.generate_dt_mask(_dt) for _dt in dt_list]

            Ch_list = [tuple(self._Ch_pool[_mask][
                np.random.choice(list(range(len(self._Ch_pool[_mask]))))].tolist())
                for _mask in dt_masks]

        self.Ch_list = Ch_list[:]
        self.min_shell_diameter = min_shell_diameter
        self.max_shell_diameter = max_shell_diameter
        self.max_shells = max_shells

        if Lz is None:
            Lz = 1.0
        self.L0 = Lz

        self.shells = \
            [SWNT(Ch, Lz=Lz, fix_Lz=True, basis=basis, bond=bond, **kwargs)
             for Ch in self.Ch_list]

        a = compute_dt(self.Ch_list[-1], bond) + dVDW
        c = compute_T(self.Ch_list[-1], bond, length=True)
        lattice = Crystal3DLattice.hexagonal(a, c)

        self.unit_cell = UnitCell(lattice, basis)

        if self.verbose:
            print(self.shells)

        self.fmtstr = "Ch_list={Ch_list!r}, Lz={Lz!r}, bond={bond!r}, " + \
            "basis={basis!r}, min_shell_diameter={min_shell_diameter!r}, " + \
            "max_shell_diameter={max_shell_diameter!r}, " + \
            "max_shells={max_shells!r}, chiral_types={chiral_types!r}, " + \
            "shell_spacing={shell_spacing!r}"

    def todict(self):
        """Return :class:`~python:dict` of `MWNT` attributes."""
        return dict(Ch_list=self.Ch_list, Lz=self.Lz, bond=self.bond,
                    basis=[self.element1, self.element2],
                    min_shell_diameter=self.min_shell_diameter,
                    max_shell_diameter=self.max_shell_diameter,
                    max_shells=self.max_shells,
                    chiral_types=self.chiral_types,
                    shell_spacing=self.shell_spacing)
