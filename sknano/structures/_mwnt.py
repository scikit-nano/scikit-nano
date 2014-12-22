# -*- coding: utf-8 -*-
"""
==============================================================================
MWNT structure class (:mod:`sknano.structures._mwnt`)
==============================================================================

.. currentmodule:: sknano.structures._mwnt

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core.refdata import dVDW

from ._base import StructureBase
from ._swnt import SWNT
from ._compute_funcs import compute_dt
from ._extras import generate_Ch_list
from ._mixins import MWNTMixin

__all__ = ['MWNT']


class MWNT(MWNTMixin, StructureBase):
    """MWNT structure class."""
    def __init__(self, Ch_list=None, Nwalls=None, Lz=None,
                 min_shell_diameter=None, max_shell_diameter=None,
                 max_shells=None, chiral_types=None, shell_spacing=dVDW,
                 **kwargs):
        if Ch_list is None and 'Ch' in kwargs:
            Ch_list = kwargs['Ch']
            del kwargs['Ch']

        super(MWNT, self).__init__(**kwargs)

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
            #self._Ch_pool = \
            #    np.asarray([(n, m) for n in range(imax) for m in range(imax)
            #                if not n == m == 0])
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
                np.random.choice(range(len(self._Ch_pool[_mask])))].tolist())
                for _mask in dt_masks]

        self.Ch_list = Ch_list[:]

        if Lz is None:
            Lz = 1.0
        self.L0 = Lz

        self.shells = \
            [SWNT(n=_Ch[0], m=_Ch[-1], Lz=Lz, fix_Lz=True, **kwargs)
             for _Ch in self.Ch_list]

        #if Lz is not None:
        #    self.nz = 10 * float(Lz) / min([shell.T for shell in self.shells])
        #elif nz is not None:

        if self.verbose:
            print(self.shells)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        """Return canonical string representation of `MWNT`."""
        strrep = "MWNT(Ch_list={!r}, Nwalls={!r}, Lz={!r}, " + \
            "element1={!r}, element2={!r}, bond={!r})"
        return strrep.format(self.Ch_list, self.Nwalls, self.Lz,
                             self.element1, self.element2, self.bond)

    def todict(self):
        """Return :class:`~python:dict` of `MWNT` attributes."""
        return dict(Ch_list=self.Ch_list, Nwalls=self.Nwalls, Lz=self.Lz,
                    element1=self.element1, element2=self.element2,
                    bond=self.bond)
