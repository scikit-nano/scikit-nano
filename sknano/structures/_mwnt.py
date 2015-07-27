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

# import numpy as np

from sknano.core.crystallography import Crystal3DLattice, UnitCell
from sknano.core.refdata import aCC, dVDW

from ._base import StructureBase
from ._compute_funcs import compute_dt, compute_T
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
    min_wall_diameter : float, optional
        Minimum `MWNT` wall diameter, in units of **Angstroms**.
    max_wall_diameter : float, optional
        Maximum `MWNT` wall diameter, in units of **Angstroms**.
    max_walls : int, optional
        Maximum number of `MWNT` walls.
    chiral_types : {None, 'armchair', 'zigzag', 'achiral', 'chiral'}, optional
        If `None`, the :attr:`~SWNT.chiral_type` of each `MWNT` walls
        will be random and determined by the set of randomly selected
        chiral indices (:attr:`~SWNT.n`, :attr:`~SWNT.m`).
    wall_spacing : float, optional
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
                 min_wall_diameter=None, max_wall_diameter=None,
                 max_walls=None, chiral_types=None, wall_spacing=dVDW,
                 basis=['C', 'C'], bond=aCC, **kwargs):
        if Ch_list is None and 'Ch' in kwargs:
            Ch_list = kwargs['Ch']
            del kwargs['Ch']

        super().__init__(basis=basis, bond=bond, **kwargs)

        if Ch_list is None or not isinstance(Ch_list, list):
            Ch_list = \
                self.generate_Ch_list(Nwalls=Nwalls, max_walls=max_walls,
                                      min_wall_diameter=min_wall_diameter,
                                      max_wall_diameter=max_wall_diameter,
                                      chiral_types=chiral_types,
                                      wall_spacing=wall_spacing)

        self.Ch_list = Ch_list[:]
        self._min_wall_diameter = min_wall_diameter
        self._max_wall_diameter = max_wall_diameter
        self._max_walls = max_walls
        self._wall_spacing = wall_spacing

        if Lz is None:
            Lz = 1.0
        self.Lz = Lz

        a = compute_dt(self.Ch_list[-1], bond=bond) + dVDW
        c = compute_T(self.Ch_list[-1], bond=bond, length=True)

        self.unit_cell = UnitCell(lattice=Crystal3DLattice.hexagonal(a, c),
                                  basis=self.basis)

        if self.verbose:
            print(self.walls)

        self.fmtstr = "Ch_list={Ch_list!r}, Lz={Lz!r}, bond={bond!r}, " + \
            "basis={basis!r}, min_wall_diameter={min_wall_diameter!r}, " + \
            "max_wall_diameter={max_wall_diameter!r}, " + \
            "max_walls={max_walls!r}, chiral_types={chiral_types!r}, " + \
            "wall_spacing={wall_spacing!r}"

    def todict(self):
        """Return :class:`~python:dict` of `MWNT` attributes."""
        return dict(Ch_list=self.Ch_list, Lz=self.Lz,
                    basis=self.basis, bond=self.bond,
                    min_wall_diameter=self.min_wall_diameter,
                    max_wall_diameter=self.max_wall_diameter,
                    max_walls=self.max_walls,
                    chiral_types=self.chiral_types,
                    wall_spacing=self.wall_spacing)
