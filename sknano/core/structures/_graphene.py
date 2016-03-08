# -*- coding: utf-8 -*-
"""
===============================================================================
Graphene structure class (:mod:`sknano.core.structures._graphene`)
===============================================================================

.. currentmodule:: sknano.core.structures._graphene

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import itertools
import numbers

import numpy as np

# from sknano.core.atoms import Atom
from sknano.core.crystallography import Crystal3DLattice, UnitCell
from sknano.core.math import Vector
from sknano.core.refdata import aCC  # , grams_per_Da
from ._base import NanoStructureBase, r_CC_vdw

__all__ = ['GraphenePrimitiveCell', 'GrapheneConventionalCell',
           'GrapheneMixin', 'Graphene', 'GrapheneBase',
           'PrimitiveCellGraphene', 'HexagonalGraphene',
           'HexagonalCellGraphene', 'ConventionalCellGraphene',
           'RectangularGraphene', 'RectangularCellGraphene',
           'GrapheneNanoribbon']


class GraphenePrimitiveCell(UnitCell):
    """Primitive graphene unit cell with 2 atom basis.

    Parameters
    ----------
    bond : :class:`~python:float`, optional
    a : :class:`~python:float`, optional
    gamma : {60, 120}, optional
    basis : {:class:`~python:list`, :class:`~sknano.core.atoms.BasisAtoms`}, \
    optional
    coords : {:class:`~python:list`}, optional
    cartesian : {:class:`~python:bool`}, optional

    """
    def __init__(self, bond=aCC, a=np.sqrt(3) * aCC, c=2 * r_CC_vdw,
                 gamma=60, basis=['C', 'C'], coords=[[0, 0, 0], [aCC, 0, 0]],
                 cartesian=True):
        lattice = Crystal3DLattice(a=a, b=a, c=c, alpha=90., beta=90.,
                                   gamma=gamma)
        lattice.rotate(angle=-np.pi / 6, axis='z')
        super().__init__(lattice, basis, coords, cartesian)


class GrapheneConventionalCell(UnitCell):
    """Conventional (rectangular) graphene unit cell with 4 atom basis.

    Parameters
    ----------
    bond : :class:`~python:float`, optional
    a, b : :class:`~python:float`, optional
    basis : {:class:`~python:list`, :class:`~sknano.core.atoms.BasisAtoms`}, \
    optional
    coords : {:class:`~python:list`}, optional
    cartesian : {:class:`~python:bool`}, optional

    """
    def __init__(self, bond=aCC, a=3 * aCC, b=np.sqrt(3) * aCC, c=2 * r_CC_vdw,
                 basis=4 * ['C'],
                 coords=[[0, 0, 0], [aCC, 0, 0],
                         [3 / 2 * aCC, np.sqrt(3) / 2 * aCC, 0],
                         [5 / 2 * aCC, np.sqrt(3) / 2 * aCC, 0]],
                 cartesian=True):
        lattice = Crystal3DLattice.orthorhombic(a, b, c)
        super().__init__(lattice, basis, coords, cartesian)


class GrapheneMixin:
    """Mixin class for graphene structure classes."""
    @property
    def n1(self):
        """Number of unit cells along :attr:`Crystal3DLattice.a1`."""
        try:
            return self._n1
        except AttributeError:
            return int(np.ceil(self.l1 / self.unit_cell.a1.length))

    @n1.setter
    def n1(self, value):
        self._n1 = int(value)

    @property
    def n2(self):
        """Number of unit cells along :attr:`Crystal3DLattice.a2`."""
        try:
            return self._n2
        except AttributeError:
            return int(np.ceil(self.l2 / self.unit_cell.a2.length))

    @n2.setter
    def n2(self, value):
        self._n2 = int(value)

    @property
    def r1(self):
        """Vector :attr:`GrapheneMixin.n1` :math:`\\times` \
            :attr:`Crystal3DLattice.a1`."""
        return self.n1 * self.unit_cell.a1

    @property
    def r2(self):
        """Vector :attr:`GrapheneMixin.n2` :math:`\\times` \
            :attr:`Crystal3DLattice.a2`."""
        return self.n2 * self.unit_cell.a2

    @property
    def area(self):
        """Total area of graphene supercell."""
        return np.abs(self.r1[:2].cross(self.r2[:2]))

    @property
    def N(self):
        """Number of graphene unit cells.

        .. math::

           N = \\frac{A_{\\mathrm{sheet}}}{A_{\\mathrm{cell}}}

        """
        return int(self.area /
                   np.abs(self.unit_cell.a1[:2].cross(self.unit_cell.a2[:2])))

    @property
    def Natoms(self):
        """Total number of atoms."""
        return self.nlayers * self.Natoms_per_layer

    @property
    def Natoms_per_layer(self):
        """Number of atoms per layer."""
        return self.N * self.Natoms_per_unit_cell

    @property
    def Natoms_per_unit_cell(self):
        """Number of atoms per unit cell."""
        return self.unit_cell.basis.Natoms


class GrapheneBase(GrapheneMixin, NanoStructureBase):
    """Graphene base structure class.

    .. versionadded:: 0.3.11

    Parameters
    ----------
    basis : {:class:`python:list`}, optional
        List of :class:`python:str`\ s of element symbols or atomic number
        of the two atom basis (default: ['C', 'C'])
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    nlayers : int, optional
        Number of graphene layers (default: 1)
    layer_spacing : float, optional
        Distance between layers in **Angstroms** (default: 3.4).
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers.
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
    verbose : bool, optional
        verbose output
    """

    def __init__(self, basis=['C', 'C'], bond=aCC, nlayers=1,
                 layer_spacing=2 * r_CC_vdw, layer_rotation_angles=None,
                 layer_rotation_increment=None, stacking_order='AB',
                 degrees=True, **kwargs):

        if 'deg2rad' in kwargs:
            degrees = kwargs['deg2rad']
            del kwargs['deg2rad']

        if 'layer_rotation_angle' in kwargs and layer_rotation_angles is None:
                layer_rotation_angles = kwargs['layer_rotation_angle']
                del kwargs['layer_rotation_angle']

        super().__init__(basis=basis, bond=bond, **kwargs)

        self.layers = []
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

        self.fmtstr = 'bond={bond!r}, basis={basis!r}, ' + \
            'nlayers={nlayers!r}, layer_spacing={layer_spacing!r}, ' + \
            'layer_rotation_angles={layer_rotation_angles!r}, ' + \
            'stacking_order={stacking_order!r}'

    def __str__(self):
        strrep = repr(self)
        strrep += '\n\n'
        strrep += 'area: {:.2f} \u212b^2\n'.format(self.area)
        strrep += 'N atoms/unit cell: {:d}\n'.format(self.Natoms_per_unit_cell)
        strrep += 'N unit cells: {:d}\n'.format(self.N)
        strrep += 'N atoms/layer: {:d}\n'.format(self.Natoms_per_layer)
        strrep += 'N layers: {:d}\n'.format(self.nlayers)
        strrep += 'N atoms: {:d}\n'.format(self.Natoms)
        return strrep

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(bond=self.bond, basis=self.basis,
                    nlayers=self.nlayers, layer_spacing=self.layer_spacing,
                    layer_rotation_angles=self.layer_rotation_angles,
                    stacking_order=self.stacking_order)


class PrimitiveCellGraphene(GrapheneBase):
    """Graphene structure class built from a primitive unit cell.

    .. versionadded:: 0.3.11

    Parameters
    ----------
    edge_length : float, optional
        length of graphene edges in **Angstroms**.

        .. versionchanged:: 0.4.0

           Changed units from nanometers to **Angstroms**

    basis : {:class:`python:list`}, optional
        List of :class:`python:str`\ s of element symbols or atomic number
        of the two atom basis (default: ['C', 'C'])
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    nlayers : int, optional
        Number of graphene layers (default: 1)
    layer_spacing : float, optional
        Distance between layers in **Angstroms** (default: 3.4).
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers.
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
    verbose : bool, optional
        verbose output

    """

    def __init__(self, edge_length=None, ncells=None, **kwargs):

        self.edge_length = edge_length
        self.l1 = self.l2 = self.edge_length

        super().__init__(**kwargs)

        self.unit_cell = GraphenePrimitiveCell(bond=self.bond,
                                               basis=self.basis)
        self.crystal_cell.scaling_matrix = [self.n1, self.n2, self.nlayers]
        self.fmtstr = 'edge_length={edge_length!r}, ' + self.fmtstr

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attr_dict = super().todict()
        attr_dict.update(dict(edge_length=self.edge_length))
        return attr_dict

HexagonalGraphene = HexagonalCellGraphene = PrimitiveCellGraphene


class ConventionalCellGraphene(GrapheneBase):
    """Graphene structure class built from a conventional unit cell.

    .. versionadded:: 0.3.11

    Parameters
    ----------
    armchair_edge_length : float, optional
        Length of armchair edge in **Angstroms**

        .. versionchanged:: 0.4.0

           Changed units from nanometers to **Angstroms**

    zigzag_edge_length : float, optional
        Length of zigzag edge in **Angstroms**

        .. versionchanged:: 0.4.0

           Changed units from nanometers to **Angstroms**

    basis : {:class:`python:list`}, optional
        List of :class:`python:str`\ s of element symbols or atomic number
        of the two atom basis (default: ['C', 'C'])
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    nlayers : int, optional
        Number of graphene layers (default: 1)
    layer_spacing : float, optional
        Distance between layers in **Angstroms** (default: 3.4).
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers.
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
    verbose : bool, optional
        verbose output

    """
    def __init__(self, armchair_edge_length=None, zigzag_edge_length=None,
                 n1=None, n2=None, **kwargs):

        if 'length' in kwargs:
            armchair_edge_length = kwargs['length']
            del kwargs['length']

        if 'width' in kwargs:
            zigzag_edge_length = kwargs['width']
            del kwargs['width']

        self.l1 = self.armchair_edge_length = armchair_edge_length
        self.l2 = self.zigzag_edge_length = zigzag_edge_length

        super().__init__(**kwargs)

        self.unit_cell = GrapheneConventionalCell(bond=self.bond,
                                                  basis=2 * self.basis)
        self.crystal_cell.scaling_matrix = [self.n1, self.n2, self.nlayers]

        self.fmtstr = 'armchair_edge_length={armchair_edge_length!r}, ' + \
            'zigzag_edge_length={zigzag_edge_length!r}, ' + self.fmtstr

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attr_dict = super().todict()
        attr_dict.update(dict(armchair_edge_length=self.armchair_edge_length,
                         zigzag_edge_length=self.zigzag_edge_length))
        return attr_dict

RectangularGraphene = RectangularCellGraphene = ConventionalCellGraphene


class GrapheneNanoribbon(ConventionalCellGraphene):
    """Graphene nanoribbon structure class.

    .. versionadded:: 0.3.11

    """
    def __init__(self, **kwargs):
        super().__init__(nlayers=1, **kwargs)


class Graphene(ConventionalCellGraphene):
    """Graphene structure class.

    .. versionchanged:: 0.3.11

       `Graphene` is now a sub-class of the `ConventionalCellGraphene`
       class to maintain backwards compatibility and also includes 2 new
       classmethods: :meth:`~Graphene.from_primitive_cell`
       and :meth:`~Graphene.from_conventional_cell`.

    Parameters
    ----------
    armchair_edge_length : float, optional
        Length of armchair edge in **Angstroms**

        .. versionadded:: 0.3.10

        .. versionchanged:: 0.4.0

           Changed units from nanometers to **Angstroms**

    zigzag_edge_length : float, optional
        Length of zigzag edge in **Angstroms**

        .. versionadded:: 0.3.10

        .. versionchanged:: 0.4.0

           Changed units from nanometers to **Angstroms**

    length : float, optional
        Length of armchair edge in **Angstroms**

        .. deprecated:: 0.3.10
           Use `armchair_edge_length` instead

    width : float, optional
        Width of graphene sheet in **Angstroms**

        .. deprecated:: 0.3.10
           Use `zigzag_edge_length` instead

    edge : {'AC', 'armchair', 'ZZ', 'zigzag'}, optional
        **A**\ rm\ **C**\ hair or **Z**\ ig\ **Z**\ ag edge along
        the `length` of the sheet.

        .. deprecated:: 0.3.10
           No longer used!

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
    nlayers : int, optional
        Number of graphene layers (default: 1)
    layer_spacing : float, optional
        Distance between layers in **Angstroms** (default: 3.4).
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers.
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
    verbose : bool, optional
        verbose output

    """
    @classmethod
    def from_primitive_cell(cls, **kwargs):
        """See the `PrimitiveCellGraphene` structure class documentation."""
        return PrimitiveCellGraphene(**kwargs)

    @classmethod
    def from_conventional_cell(cls, **kwargs):
        """See the `ConventionalCellGraphene` structure class documentation."""
        return ConventionalCellGraphene(**kwargs)
