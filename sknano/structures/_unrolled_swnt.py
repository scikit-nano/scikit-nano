# -*- coding: utf-8 -*-
"""
==============================================================================
Unrolled SWNT structure class (:mod:`sknano.structures._unrolled_swnt`)
==============================================================================

.. currentmodule:: sknano.structures._unrolled_swnt

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from sknano.core.atoms import BasisAtom as Atom
from sknano.core.crystallography import Crystal3DLattice, UnitCell
from sknano.core.refdata import aCC, dVDW

from ._mixins import NanotubeMixin, UnrolledSWNTMixin
from ._base import StructureBase
from ._compute_funcs import compute_dt, compute_T

__all__ = ['UnrolledSWNT']


class UnrolledSWNT(UnrolledSWNTMixin, NanotubeMixin, StructureBase):
    """Unrolled SWNT structure class."""

    def __init__(self, *Ch, nx=1, nz=1, bond=aCC, basis=['C', 'C'],
                 nlayers=1, layer_spacing=dVDW, stacking_order='AB',
                 Lx=None, fix_Lx=False, Lz=None, fix_Lz=False, **kwargs):

        try:
            n, m = Ch
        except ValueError:
            try:
                n, m = Ch[0]
            except IndexError:
                n = kwargs['n']
                del kwargs['n']
                m = kwargs['m']
                del kwargs['m']

        a = compute_dt(n, m, bond) + dVDW
        c = compute_T(n, m, bond, length=True)
        lattice = Crystal3DLattice.hexagonal(a, c)

        super().__init__(bond=bond, **kwargs)

        self.unit_cell = UnitCell(lattice, basis)

        self.n = n
        self.m = m

        self.fix_Lx = fix_Lx
        if Lx is not None:
            self.nx = 10 * float(Lx) / self.Ch
        elif nx is not None:
            self.nx = nx
        else:
            self.nx = 1

        self.fix_Lz = fix_Lz
        if Lz is not None:
            self.nz = 10 * float(Lz) / self.T
        elif nz is not None:
            self.nz = nz
        else:
            self.nz = 1

        self.nlayers = nlayers
        self.layer_spacing = layer_spacing
        self.stacking_order = stacking_order
        self.generate_unit_cell()

        self.fmtstr = "{Ch!r}, nx={nx!r}, nz={nz!r}, bond={bond!r}, " + \
            "basis={basis!r}, nlayers={nlayers!r}, " + \
            "layer_spacing={layer_spacing!r}, " + \
            "stacking_order={stacking_order!r}"

    def generate_unit_cell(self):
        """Generate the nanotube unit cell."""
        eps = 0.01

        e1 = self.element1
        e2 = self.element2
        N = self.N
        T = self.T
        rt = self.rt

        psi, tau, dpsi, dtau = self.unit_cell_symmetry_params

        self.basis.clear()

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

                x = rt * theta
                z = h

                while z > T - eps:
                    z -= T

                if z < 0:
                    z += T

                if self.debug:
                    print('i={}: x, z = ({:.6f}, {:.6f})'.format(i, x, z))

                atom = Atom(element, x=x, z=z)
                atom.rezero()

                if self.verbose:
                    print('Basis Atom:\n{}'.format(atom))

                self.basis.append(atom)

    def todict(self):
        """Return :class:`~python:dict` of `SWNT` attributes."""
        return dict(Ch=(self.n, self.m), nx=self.nx, nz=self.nz,
                    bond=self.bond, basis=[self.element1, self.element2],
                    nlayers=self.nlayers, layer_spacing=self.layer_spacing,
                    stacking_order=self.stacking_order)
