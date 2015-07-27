# -*- coding: utf-8 -*-
"""
==============================================================================
SWNT structure class (:mod:`sknano.structures._swnt`)
==============================================================================

.. currentmodule:: sknano.structures._swnt

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core.atoms import BasisAtom as Atom
from sknano.core.crystallography import Crystal3DLattice, UnitCell
from sknano.core.refdata import aCC, dVDW
from ._base import StructureBase
from ._compute_funcs import compute_dt, compute_T
from ._extras import attr_strfmt, attr_symbols, attr_units, get_chiral_indices
from ._mixins import SWNTMixin

__all__ = ['SWNT', 'Nanotube']


class SWNT(SWNTMixin, StructureBase):
    """SWNT structure class.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of integers (i.e., *Ch = ((n, m)) or
        2 integers (i.e., *Ch = (n, m) specifying the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nz : :class:`python:int`, optional
        Number of repeat unit cells in the :math:`z` direction, along
        the *length* of the nanotube.
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
        nearest neighbor atoms. Must be in units of **Angstroms**.
    Lz : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the `nz` value.

        .. versionadded:: 0.2.5

    tube_length : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the `nz` value.

        .. deprecated:: 0.2.5
           Use `Lz` instead

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.

        .. versionadded:: 0.2.6

    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.structures import SWNT

    Create a SWNT with :math:`\\mathbf{C}_{h} = (10, 10)` chirality.

    >>> swnt = SWNT((10, 10), verbose=True)
    >>> print(swnt)
    SWNT((10, 10), nz=1)
    n: 10
    m: 10
    t₁: 1
    t₂: -1
    d: 10
    dR: 30
    N: 20
    R: (1, 0)
    θc: 30.00°
    Ch: 42.60 Å
    T: 2.46 Å
    dt: 13.56 Å
    rt: 6.78 Å
    electronic_type: metallic

    Change the chirality to :math:`\\mathbf{C}_{h} = (20, 10)`.

    >>> swnt.n = 20
    >>> print(swnt)
    SWNT((20, 10), nz=1)
    n: 20
    m: 10
    t₁: 4
    t₂: -5
    d: 10
    dR: 10
    N: 140
    R: (1, -1)
    θc: 19.11°
    Ch: 65.07 Å
    T: 11.27 Å
    dt: 20.71 Å
    rt: 10.36 Å
    electronic_type: semiconducting, type 2

    Change the chirality to :math:`\\mathbf{C}_{h} = (20, 0)`.

    >>> swnt.m = 0
    >>> print(swnt)
    SWNT((20, 0), nz=1)
    n: 20
    m: 0
    t₁: 1
    t₂: -2
    d: 20
    dR: 20
    N: 40
    R: (1, -1)
    θc: 0.00°
    Ch: 49.19 Å
    T: 4.26 Å
    dt: 15.66 Å
    rt: 7.83 Å
    electronic_type: semiconducting, type 1

    """
    # add each attribute in the order I want them to appear in
    # verbose output mode
    _structure_attrs = ['n', 'm', 't1', 't2', 'd', 'dR', 'N', 'R',
                        'chiral_angle', 'Ch', 'T', 'dt', 'rt',
                        'electronic_type']

    def __init__(self, *Ch, nz=None, basis=['C', 'C'], bond=aCC,
                 Lz=None, fix_Lz=False, **kwargs):

        n, m, kwargs = get_chiral_indices(*Ch, **kwargs)

        if 'tube_length' in kwargs:
            Lz = kwargs['tube_length']
            del kwargs['tube_length']

        super().__init__(basis=basis, bond=bond, **kwargs)

        a = compute_dt(n, m, bond=bond) + dVDW
        c = compute_T(n, m, bond=bond, length=True)
        self.unit_cell = UnitCell(lattice=Crystal3DLattice.hexagonal(a, c),
                                  basis=self.basis)

        self.L0 = Lz  # store initial value of Lz

        self.n = n
        self.m = m

        self.fix_Lz = fix_Lz
        if Lz is not None:
            self.nz = 10 * float(Lz) / self.T
        elif nz is not None:
            self.nz = nz
        else:
            self.nz = 1

        self.generate_unit_cell()

        fmtstr = "{Ch!r}, "
        if self.fix_Lz:
            fmtstr += "Lz={Lz!r}, fix_Lz=True, "
        else:
            fmtstr += "nz={nz!r}, "

        self.fmtstr = fmtstr + "basis={basis!r}, bond={bond!r}"

    def __str__(self):
        """Return nice string representation of `SWNT`."""
        fmtstr = repr(self)
        if self.verbose:
            fmtstr += '\n'
            for attr in self._structure_attrs:
                var = attr
                if attr in attr_symbols:
                    var = attr_symbols[attr]
                if attr in attr_strfmt:
                    if attr in attr_units:
                        fmtstr += \
                            "{}: {}{}\n".format(
                                var, attr_strfmt[attr].format(
                                    getattr(self, attr)), attr_units[attr])
                    else:
                        fmtstr += "{}: {}\n".format(
                            var, attr_strfmt[attr].format(getattr(self, attr)))
                else:
                    if attr in attr_units:
                        fmtstr += "{}: {}{}\n".format(
                            var, getattr(self, attr), attr_units[attr])
                    else:
                        fmtstr += "{}: {}\n".format(
                            var, getattr(self, attr))

        return fmtstr

    def generate_unit_cell(self):
        """Generate the nanotube unit cell."""
        eps = 0.01

        e1 = self.element1
        e2 = self.element2
        N = self.N
        T = self.T
        rt = self.rt

        psi, tau, dpsi, dtau = self.unit_cell_symmetry_params

        lattice = self.unit_cell.lattice
        self.unit_cell.basis.clear()
        if self.verbose:
            print('dpsi: {}'.format(dpsi))
            print('dtau: {}\n'.format(dtau))

        for i in range(N):
            for j, element in enumerate((e1, e2), start=1):
                theta = i * psi
                h = i * tau

                if j == 2:
                    theta += dpsi
                    h -= dtau

                x = rt * np.cos(theta)
                y = rt * np.sin(theta)
                z = h

                while z > T - eps:
                    z -= T

                if z < 0:
                    z += T

                xs, ys, zs = lattice.cartesian_to_fractional([x, y, z])
                if self.debug:
                    print('i={}: x, y, z = ({:.6f}, {:.6f}, {:.6f})'.format(
                        i, x, y, z))

                    print('xs, ys, zs = ({:.6f}, {:.6f}, {:.6f})'.format(
                        xs, ys, zs))

                # atom = Atom(element, lattice=lattice, x=x, y=y, z=z)
                # print('i={}: x, y, z = ({:.6f}, {:.6f}, {:.6f})'.format(
                #     i, x, y, z))
                atom = Atom(element, lattice=lattice, xs=xs, ys=ys, zs=zs)
                atom.rezero()

                if self.verbose:
                    print('Basis Atom:\n{}'.format(atom))

                self.unit_cell.basis.append(atom)

    def todict(self):
        """Return :class:`~python:dict` of `SWNT` attributes."""
        return dict(Ch=(self.n, self.m), nz=self.nz,
                    bond=self.bond, basis=self.basis,
                    Lz=self.Lz)
Nanotube = SWNT
