# -*- coding: utf-8 -*-
"""
==============================================================================
SWNT structure class (:mod:`sknano.structures._swnt`)
==============================================================================

.. currentmodule:: sknano.structures._swnt

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._base import StructureBase
from ._compute_funcs import compute_dt, compute_rt, compute_tube_mass, \
    compute_linear_mass_density
from ._mixins import NanotubeMixin

__all__ = ['SWNT', 'Nanotube']


class SWNT(NanotubeMixin, StructureBase):
    u"""Class for creating interactive SWNT objects.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nz : int, optional
        Number of repeat unit cells along the :math:`z` dimension
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b
    Lz : float, optional
        Length of the nanotube in **nanometers**.
        Overrides the :math:`n_z` cell values.

        .. versionadded:: 0.2.5

    tube_length : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides `nz` value.

        .. deprecated:: 0.2.5
           Use `Lz` instead

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.

    Examples
    --------

    >>> from sknano.structures import SWNT

    Create a SWNT with :math:`\\mathbf{C}_{h} = (10, 10)` chirality.

    >>> nt = SWNT(n=10, m=10, verbose=True)
    n: 10
    m: 10
    t₁: 1
    t₂: -1
    d: 10
    dR: 30
    N: 20
    R: (1, 0)
    θc: 30.00°
    Ch: 42.63 Å
    T: 2.46 Å
    dt: 13.57 Å
    rt: 6.78 Å

    Change the chirality to :math:`\\mathbf{C}_{h} = (20, 10)`.

    >>> nt.n = 20
    n: 20
    m: 10
    t₁: 4
    t₂: -5
    d: 10
    dR: 10
    N: 140
    R: (1, -1)
    θc: 19.11°
    Ch: 65.12 Å
    T: 11.28 Å
    dt: 20.73 Å
    rt: 10.36 Å

    Change the chirality to :math:`\\mathbf{C}_{h} = (20, 0)`.

    >>> nt.m = 0
    n: 20
    m: 0
    t₁: 1
    t₂: -2
    d: 20
    dR: 20
    N: 40
    R: (1, -1)
    θc: 0.00°
    Ch: 49.22 Å
    T: 4.26 Å
    dt: 15.67 Å
    rt: 7.83 Å

    """
    def __init__(self, n=10, m=0, nz=1, tube_length=None, Lz=None,
                 fix_Lz=False, **kwargs):

        super(SWNT, self).__init__(**kwargs)

        if tube_length is not None and Lz is None:
            Lz = tube_length

        self.L0 = Lz  # store initial value of Lz

        # add each parameter in the order I want them to appear in
        # verbose output mode
        self._params = ['n', 'm', 't1', 't2', 'd', 'dR', 'N',
                        'R', 'chiral_angle', 'Ch', 'T', 'dt', 'rt',
                        'electronic_type']

        self.n = n
        self.m = m

        self.fix_Lz = fix_Lz
        if Lz is not None:
            self.nz = 10 * float(Lz) / self.T
        elif nz is not None:
            self.nz = nz
        else:
            self.nz = 1

    def __str__(self):
        """Return nice string representation of `SWNT`."""
        return repr(self)

    def __repr__(self):
        """Return canonical string representation of `SWNT`."""
        strrep = "SWNT(n={!r}, m={!r}, element1={!r}, element2={!r}, bond={!r}"
        if self.fix_Lz:
            strrep += ", Lz={!r}, fix_Lz={!r})"
            return strrep.format(self.n, self.m, self.element1, self.element2,
                                 self.bond, self.Lz, self.fix_Lz)
        else:
            strrep += ", nz={!r})"
            return strrep.format(self.n, self.m, self.element1, self.element2,
                                 self.bond, self.nz)

    @property
    def Natoms_per_tube(self):
        """Number of atoms in nanotube :math:`N_{\\mathrm{atoms/tube}}`."""
        return self.Natoms

    @property
    def linear_mass_density(self):
        """Linear mass density of nanotube in g/nm."""
        return compute_linear_mass_density(self.n, self.m, bond=self.bond,
                                           element1=self.element1,
                                           element2=self.element2)

    @property
    def Ntubes(self):
        """Number of nanotubes."""
        return 1

    @property
    def tube_length(self):
        """Alias for :attr:`SWNT.Lz`"""
        return self.Lz

    @property
    def tube_mass(self):
        """SWNT mass in **grams**."""
        return compute_tube_mass(self.n, self.m, self.nz,
                                 element1=self.element1,
                                 element2=self.element2)

    def todict(self):
        return dict(n=self.n, m=self.m, nz=self.nz, Lz=self.Lz,
                    fix_Lz=self.fix_Lz, element1=self.element1,
                    element2=self.element2, bond=self.bond)

Nanotube = SWNT
