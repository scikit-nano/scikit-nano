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

# from sknano.core.crystallography import Crystal3DLattice, UnitCell
from sknano.core.refdata import aCC, element_data

from ._base import NanoStructureBase, r_CC_vdw
from ._swnt import SWNT, compute_dt  # , compute_T
from ._extras import generate_Ch_list

__all__ = ['MWNTMixin', 'MWNT']


class MWNTMixin:
    """Mixin class for MWNTs."""
    @property
    def Ch_list(self):
        return self._Ch_list

    @Ch_list.setter
    def Ch_list(self, value):
        if not isinstance(value, list):
            raise TypeError('Expected a list')
        self._Ch_list = value[:]

    @property
    def chiral_types(self):
        """List of chiral types for each `MWNT` wall."""
        return [swnt.chiral_type for swnt in self.walls]

    @chiral_types.setter
    def chiral_types(self, value):
        if not isinstance(value, list):
            raise TypeError('Expected a list')
        self.update_Ch_list(chiral_types=value)

    @property
    def chiral_set(self):
        """Set of all chiral types in `MWNT`."""
        return set(self.chiral_types)

    @property
    def dt(self):
        """`MWNT` wall diameters :math:`d_t=\\frac{|\\mathbf{C}_h|}{\\pi}` \
        in \u212b."""
        return self.walls[-1].dt

    @property
    def rt(self):
        """`MWNT` wall radii :math:`r_t=\\frac{|\\mathbf{C}_h|}{2\\pi}` \
        in \u212b."""
        return self.walls[-1].rt

    @property
    def Natoms(self):
        """Number of atoms in `MWNT`.

           **Returns total number of atoms in `MWNT`.**
           Use :attr:`~MWNT.Natoms_per_wall` to get a list of the number of
           atoms in each `MWNT` wall.

        .. math::

           N_{\\mathrm{atoms}} = \\sum_{\\mathrm{walls}}

        """
        return np.asarray(self.Natoms_per_wall).sum()

    @property
    def Natoms_per_tube(self):
        """Number of atoms in `MWNT`."""
        return self.Natoms

    @property
    def Ntubes(self):
        """Number of `MWNT`\ s."""
        return 1

    @property
    def Nwalls(self):
        """Number of `MWNT` walls."""
        return len(self.Ch_list)

    @Nwalls.setter
    def Nwalls(self, value):
        self.update_Ch_list(Nwalls=value)

    @property
    def min_wall_diameter(self):
        return self._min_wall_diameter

    @min_wall_diameter.setter
    def min_wall_diameter(self, value):
        self._min_wall_diameter = value
        self.update_Ch_list()

    @property
    def max_wall_diameter(self):
        return self._max_wall_diameter

    @max_wall_diameter.setter
    def max_wall_diameter(self, value):
        self._max_wall_diameter = value
        self.update_Ch_list()

    @property
    def max_walls(self):
        return self._max_walls

    @max_walls.setter
    def max_walls(self, value):
        self._max_walls = value

    @property
    def wall_spacing(self):
        return self._wall_spacing

    @wall_spacing.setter
    def wall_spacing(self, value):
        self._wall_spacing = value
        self.update_Ch_list()

    @property
    def tube_mass(self):
        """MWNT mass in **grams**."""
        return np.asarray([swnt.tube_mass for swnt in self.walls]).sum()

    # @property
    # def Lz(self):
    #     return self._Lz

    # @Lz.setter
    # def Lz(self, value):
    #     self._Lz = value

    @property
    def Natoms_per_wall(self):
        """Alias for :attr:`MWNT.Natoms_list`"""
        return self.Natoms_list

    @property
    def Natoms_list(self):
        """List of `MWNT` `SWNT` wall's number of atoms \
        :attr:`~SWNT.Natoms`."""
        return [swnt.Natoms for swnt in self.walls]

    @property
    def nz_list(self):
        """Number of nanotube unit cells along the :math:`z`-axis."""
        return [swnt.nz for swnt in self.walls]

    @property
    def Lz_list(self):
        """MWNT length :math:`L_z = L_{\\mathrm{tube}}` in **nanometers**."""
        return [swnt.Lz for swnt in self.walls]

    @property
    def T_list(self):
        """Length of `MWNT` unit cell :math:`|\\mathbf{T}|` in \u212b."""
        return [swnt.T for swnt in self.walls]

    @property
    def dt_list(self):
        """List of `MWNT` `SWNT` wall diameters :attr:`~SWNT.dt` \
        :math:`d_t=\\frac{|\\mathbf{C}_h|}{\\pi}` in \u212b."""
        return [swnt.dt for swnt in self.walls]

    @property
    def rt_list(self):
        """List of `MWNT` `SWNT` wall radii :attr:`~SWNT.rt` \
        :math:`r_t=\\frac{|\\mathbf{C}_h|}{2\\pi}` in \u212b."""
        return [swnt.rt for swnt in self.walls]

    @property
    def wall_diameters(self):
        """Alias for :attr:`MWNTMixin.dt_list`."""
        return self.dt_list

    @property
    def wall_radii(self):
        """Alias for :attr:`MWNTMixin.rt_list`."""
        return self.rt_list

    @property
    def walls(self):
        """List of `MWNT` `SWNT` wall structures."""
        return [SWNT(Ch, Lz=self.Lz, fix_Lz=True, basis=self.basis,
                     bond=self.bond) for Ch in self.Ch_list]

    def get_wall(self, Ch):
        """Return the :class:`~sknano.structures.SWNT` structure with \
            chirality `Ch`.

        """
        return SWNT(Ch, Lz=self.Lz, fix_Lz=True, basis=self.basis,
                    bond=self.bond) if Ch in self.Ch_list else None

    def generate_dt_mask(self, dt, max_dt_diff=0.5):
        """Generate boolean mask array.

        Parameters
        ----------
        dt : float
        max_dt_diff : float, optional

        Returns
        -------
        dt_mask : :class:`~numpy:numpy.ndarray`

        """
        dt_mask = np.abs(self._dt_pool - dt) <= max_dt_diff
        while not np.any(dt_mask):
            max_dt_diff += max_dt_diff
            dt_mask = np.abs(self._dt_pool - dt) <= max_dt_diff
        return dt_mask

    def generate_Ch_list(self, Nwalls=None, max_walls=None,
                         min_wall_diameter=None, max_wall_diameter=None,
                         chiral_types=None, wall_spacing=None):
        if Nwalls is not None:
            max_walls = Nwalls

        if max_walls is None:
            max_walls = 10

        if max_wall_diameter is None:
            max_wall_diameter = np.inf

        if min_wall_diameter is None:
            min_wall_diameter = 5.0

        if wall_spacing is None:
            wall_spacing = 2 * element_data['C']['VanDerWaalsRadius']

        delta_dt = 2 * wall_spacing

        imax = 100

        self._Ch_pool = \
            np.asarray(generate_Ch_list(imax=imax,
                                        chiral_types=chiral_types))
        self._dt_pool = np.asarray([compute_dt(_Ch, bond=self.bond) for _Ch
                                   in self._Ch_pool])

        dt_mask = np.logical_and(self._dt_pool >= min_wall_diameter,
                                 self._dt_pool <= max_wall_diameter)

        self._Ch_pool = self._Ch_pool[dt_mask]
        self._dt_pool = self._dt_pool[dt_mask]

        if max_wall_diameter < np.inf:
            dt_list = []
            dt = self._dt_pool.min()
            while dt <= max_wall_diameter and len(dt_list) < max_walls:
                dt_list.append(dt)
                dt += delta_dt
        else:
            dt_list = [self._dt_pool.min() + i * delta_dt
                       for i in range(max_walls)]

        dt_masks = [self.generate_dt_mask(_dt) for _dt in dt_list]

        return [tuple(self._Ch_pool[_mask][np.random.choice(
            list(range(len(self._Ch_pool[_mask]))))].tolist())
            for _mask in dt_masks]

    def update_Ch_list(self, Nwalls=None, min_wall_diameter=None,
                       max_wall_diameter=None, wall_spacing=None,
                       chiral_types=None):
        if Nwalls is None:
            Nwalls = self.Nwalls
        if min_wall_diameter is None:
            min_wall_diameter = self.min_wall_diameter
        if max_wall_diameter is None:
            max_wall_diameter = self.max_wall_diameter
        if wall_spacing is None:
            wall_spacing = self.wall_spacing
        self.Ch_list = \
            self.generate_Ch_list(Nwalls=Nwalls,
                                  min_wall_diameter=min_wall_diameter,
                                  max_wall_diameter=max_wall_diameter,
                                  chiral_types=chiral_types,
                                  wall_spacing=wall_spacing)


class MWNT(MWNTMixin, NanoStructureBase):
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
                 max_walls=None, chiral_types=None, wall_spacing=2 * r_CC_vdw,
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

        self.unit_cell = self.get_wall(self.Ch_list[-1]).unit_cell

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
