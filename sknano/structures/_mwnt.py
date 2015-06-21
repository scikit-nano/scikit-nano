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
    """MWNT structure class.

    Parameters
    ----------
    Ch_list : :class:`python:list`, optional
        (:attr:`~SWNT.n`, :attr:`~SWNT.m`) for each `SWNT` wall in `MWNT`.
    Nwalls : int, optional
        Number of `SWNT` walls in `MWNT`.
    Lz : float, optional
        `MWNT` length in **nanometers**.
    min_shell_diameter : float, optional
        Minimum `MWNT` wall diameter, in units of **Angstroms**.
    max_shell_diameter : float, optional
        Maximum `MWNT` wall diameter, in units of **Angstroms**.
    max_shells : int, optional
        Maximum number of `MWNT` walls.
    chiral_types : {None, 'armchair', 'zigzag', 'achiral', 'chiral'}, optional
        If `None`, the :attr:`~SWNT.chiral_type` of each `MWNT` walls
        will be random and determined by the set of randomly selected
        chiral indices (:attr:`~SWNT.n`, :attr:`~SWNT.m`).
    shell_spacing : float, optional
        Inter-wall spacing in units of **Angstroms**.
        Default value is the van der Waals interaction distance of 3.35
        Angstroms.
    basis : {:class:`python:list`}, optional
        List of :class:`python:str`\ s of element symbols or atomic number
        of the two atom basis (default: ['C', 'C'])

        .. versionadded:: 0.3.10

    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.Atom` 1 and 2

        .. deprecated:: 0.3.10
           Use `basis` instead

    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms, in units of **Angstroms**.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.generators import MWNT

    """
    def __init__(self, Ch_list=None, Nwalls=None, Lz=None,
                 min_shell_diameter=None, max_shell_diameter=None,
                 max_shells=None, chiral_types=None, shell_spacing=dVDW,
                 basis=['C', 'C'], bond=aCC, **kwargs):
        if Ch_list is None and 'Ch' in kwargs:
            Ch_list = kwargs['Ch']
            del kwargs['Ch']

        super().__init__(basis=basis, bond=bond, **kwargs)

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
        self.shell_spacing = shell_spacing

        if Lz is None:
            Lz = 1.0
        self.L0 = Lz

        self.shells = \
            [SWNT(Ch, Lz=Lz, fix_Lz=True, basis=self.basis, bond=self.bond,
                  **kwargs)
             for Ch in self.Ch_list]

        a = compute_dt(self.Ch_list[-1], bond=bond) + dVDW
        c = compute_T(self.Ch_list[-1], bond=bond, length=True)
        lattice = Crystal3DLattice.hexagonal(a, c)

        self.unit_cell = UnitCell(lattice, basis=self.basis)

        if self.verbose:
            print(self.shells)

        self.fmtstr = "Ch_list={Ch_list!r}, Lz={Lz!r}, bond={bond!r}, " + \
            "basis={basis!r}, min_shell_diameter={min_shell_diameter!r}, " + \
            "max_shell_diameter={max_shell_diameter!r}, " + \
            "max_shells={max_shells!r}, chiral_types={chiral_types!r}, " + \
            "shell_spacing={shell_spacing!r}"

    def todict(self):
        """Return :class:`~python:dict` of `MWNT` attributes."""
        return dict(Ch_list=self.Ch_list, Lz=self.Lz,
                    basis=self.basis, bond=self.bond,
                    min_shell_diameter=self.min_shell_diameter,
                    max_shell_diameter=self.max_shell_diameter,
                    max_shells=self.max_shells,
                    chiral_types=self.chiral_types,
                    shell_spacing=self.shell_spacing)
