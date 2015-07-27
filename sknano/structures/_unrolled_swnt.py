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

import numbers

import numpy as np

from sknano.core.atoms import BasisAtom as Atom
from sknano.core.crystallography import Crystal3DLattice, UnitCell
from sknano.core.math import Vector
from sknano.core.refdata import aCC, dVDW

from ._base import StructureBase
from ._compute_funcs import compute_dt, compute_T
from ._extras import get_chiral_indices
from ._mixins import NanotubeMixin, UnrolledSWNTMixin

__all__ = ['UnrolledSWNT']


class UnrolledSWNT(UnrolledSWNTMixin, NanotubeMixin, StructureBase):
    """Unrolled SWNT structure class.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of integers (i.e., *Ch = ((n, m)) or
        2 integers (i.e., *Ch = (n, m) specifying the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nx : :class:`python:int`, optional
        Number of repeat unit cells in the :math:`x` direction, along
        the *unrolled* chiral vector.
    nz : :class:`python:int`, optional
        Number of repeat unit cells in the :math:`z` direction, along
        the *length* of the nanotube.
    basis : {:class:`python:list`}, optional
        List of :class:`python:str`\ s of element symbols or atomic number
        of the two atom basis (default: ['C', 'C'])

        .. versionadded:: 0.3.10

    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.Atom`\ s 1 and 2

        .. deprecated:: 0.3.10
           Use `basis` instead

    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    nlayers : int, optional
        Number of layers (default: 1)
    layer_spacing : float, optional
        Distance between layers in **Angstroms** (default: 3.35).
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of layers.
    layer_rotation_angles : list, optional
        list of rotation angles for each layer in **degrees** if
        `degrees` is `True` (default), otherwise in radians.
        The list length must equal the number of layers.
    layer_rotation_increment : float, optional
        incremental layer rotation angle in **degrees** if
        `degrees` is `True` (default), otherwise in radians.
        Each subsequent layer will
        be rotated by `layer_rotation_increment` relative to the layer
        below it.
    Lx : float, optional
        Length of the unrolled swnt sheet along the chiral vector
        in units of **nanometers**. Overrides the `nx` value.
    fix_Lx : bool, optional
        Generate the unrolled swnt sheet with the length along the
        chiral vector as close to the specified :math:`L_x` as possible.
        If `True`, then non integer :math:`n_x` cells are permitted.
    Lz : float, optional
        Length of the unrolled swnt sheet along the translation vector
        in units of **nanometers**. Overrides the `nz` value.
    fix_Lz : bool, optional
        Generate the unrolled swnt sheet with the length along the
        translation vector as close to the specified :math:`L_z` as possible.
        If `True`, then non integer :math:`n_z` cells are permitted.

    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.structures import UnrolledSWNT

    """

    def __init__(self, *Ch, nx=1, nz=1, basis=['C', 'C'], bond=aCC,
                 nlayers=1, layer_spacing=dVDW, layer_rotation_angles=None,
                 layer_rotation_increment=None, degrees=True,
                 stacking_order='AB', Lx=None, fix_Lx=False,
                 Lz=None, fix_Lz=False, **kwargs):

        n, m, kwargs = get_chiral_indices(*Ch, **kwargs)

        super().__init__(basis=basis, bond=bond, **kwargs)

        a = compute_dt(n, m, bond=bond) + dVDW
        c = compute_T(n, m, bond=bond, length=True)

        self.unit_cell = UnitCell(lattice=Crystal3DLattice.hexagonal(a, c),
                                  basis=self.basis)

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

        if layer_rotation_increment is not None and \
                layer_rotation_angles is None:
            layer_rotation_angles = layer_rotation_increment * \
                np.arange(self.nlayers)
        elif isinstance(layer_rotation_angles, numbers.Number):
            layer_rotation_angles = layer_rotation_angles * \
                np.ones(self.nlayers)
        elif layer_rotation_angles is None or \
                isinstance(layer_rotation_angles, (tuple, list, np.ndarray)) \
                and len(layer_rotation_angles) != self.nlayers:
            layer_rotation_angles = np.zeros(self.nlayers)
            degrees = False

        if degrees:
            layer_rotation_angles = \
                np.radians(np.asarray(layer_rotation_angles)).tolist()

        self.layer_rotation_angles = \
            np.asarray(layer_rotation_angles).tolist()

        self.layer_shift = Vector()
        self.stacking_order = stacking_order
        if nlayers > 1 and stacking_order == 'AB':
            self.layer_shift.x = self.bond

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

                x = rt * theta
                z = h

                while z > T - eps:
                    z -= T

                if z < 0:
                    z += T

                xs, ys, zs = lattice.cartesian_to_fractional([x, 0, z])
                if self.debug:
                    print('i={}: x, z = ({:.6f}, {:.6f})'.format(i, x, z))

                atom = Atom(element, lattice=lattice, xs=xs, ys=ys, zs=zs)
                atom.rezero()

                if self.verbose:
                    print('Basis Atom:\n{}'.format(atom))

                self.unit_cell.basis.append(atom)

    def todict(self):
        """Return :class:`~python:dict` of `SWNT` attributes."""
        return dict(Ch=(self.n, self.m), nx=self.nx, nz=self.nz,
                    bond=self.bond, basis=self.basis,
                    nlayers=self.nlayers, layer_spacing=self.layer_spacing,
                    stacking_order=self.stacking_order)
