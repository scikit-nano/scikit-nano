# -*- coding: utf-8 -*-
"""
==============================================================================
Graphene structure generators (:mod:`sknano.generators._graphene_generator`)
==============================================================================

.. currentmodule:: sknano.generators._graphene_generator

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import copy

import numpy as np

from sknano.core.math import Vector
from sknano.core.crystallography import SuperCell
from sknano.core.structures import PrimitiveCellGraphene, \
    ConventionalCellGraphene
from ._base import Atom, Atoms, NanoStructureGenerator

__all__ = ['GrapheneGenerator',
           'GrapheneGeneratorBase',
           'PrimitiveCellGrapheneGenerator',
           'HexagonalGrapheneGenerator',
           'HexagonalCellGrapheneGenerator',
           'ConventionalCellGrapheneGenerator',
           'RectangularGrapheneGenerator',
           'RectangularCellGrapheneGenerator']


class GrapheneGeneratorBase(NanoStructureGenerator):
    """`N`-layer graphene generator class.

    Parameters
    ----------
    basis : {:class:`python:list`}, optional
        List of :class:`python:str`\ s of element symbols or atomic number
        of the two atom basis (default: ['C', 'C'])
    bond : float, optional
        bond length between nearest-neighbor atoms in **Angstroms**.
    nlayers : int, optional
        Number of graphene layers.
    layer_spacing : float, optional
        Distance between layers in **Angstroms**.
    stacking_order : {'AB', 'AA'}, optional
        Stacking order of graphene layers
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
    autogen : bool, optional
        automatically generate unit cell and full structure
    verbose : bool, optional
        verbose output

    """
    def generate(self, finalize=True):
        """Generate the full structure coordinates."""
        self.structure.clear()
        layer0 = Atoms()
        supercell = SuperCell(self.unit_cell, [self.n1, self.n2, 1])

        for atom in supercell:
            layer0.append(Atom(**atom.todict()))
        layer0.center_centroid()

        lattice_shift = Vector(p0=supercell.basis.centroid,
                               p=layer0.centroid)
        lattice_shift.z = self.nlayers * lattice_shift.z
        self.lattice_shift = Vector(lattice_shift)
        # self.Natoms_per_layer = layer0.Natoms

        self.layers = []
        for nlayer in range(self.nlayers):
            layer = copy.deepcopy(layer0)
            layer.translate(Vector([0, 0, nlayer * self.layer_spacing]))
            if (nlayer % 2) != 0:
                layer.translate(self.layer_shift)
            [setattr(atom, 'mol', nlayer + 1) for atom in layer]
            layer.rotate(angle=self.layer_rotation_angles[nlayer], axis='z')
            self.atoms.extend(layer)
            self.layers.append(layer)

        if finalize:
            self.finalize()

    @classmethod
    def generate_fname(cls, nlayers=None, basis=None, **kwargs):
        """Generate a filename string."""
        nlayer = '{}layer'.format(nlayers)
        fname_wordlist = [nlayer, 'graphene']
        basis = '-'.join(basis)
        if not basis == 'C-C':
            fname_wordlist.append('_'.join((basis, 'basis')))
        fname = '_'.join(fname_wordlist)
        return fname

    def save(self, fname=None, outpath=None, structure_format=None,
             center_centroid=True, rotation_angle=np.pi / 2, rotation_axis='x',
             **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save` method
        for documentation.

        """
        if fname is None:
            fname = \
                self.generate_fname(nlayers=self.nlayers, basis=self.basis)

        if center_centroid and self.nlayers > 1:
            self.atoms.center_centroid()

        super().save(fname=fname, outpath=outpath,
                     structure_format=structure_format, center_centroid=False,
                     angle=rotation_angle, axis=rotation_axis, **kwargs)


class PrimitiveCellGrapheneGenerator(GrapheneGeneratorBase,
                                     PrimitiveCellGraphene):
    """`N`-layer graphene generator class using a primitive unit cell.

    Parameters
    ----------
    edge_length : float
        Length of graphene edges in **Angstroms**
    basis : {:class:`python:list`}, optional
        List of :class:`python:str`\ s of element symbols or atomic number
        of the two atom basis (default: ['C', 'C'])
    bond : float, optional
        bond length between nearest-neighbor atoms in **Angstroms**.
    nlayers : int, optional
        Number of graphene layers.
    layer_spacing : float, optional
        Distance between layers in **Angstroms**.
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers
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
    autogen : bool, optional
        automatically generate unit cell and full structure
    verbose : bool, optional
        verbose output

    Examples
    --------

    >>> from sknano.generators import PrimitiveCellGrapheneGenerator
    >>> graphene = PrimitiveCellGrapheneGenerator(edge_length=10)
    >>> graphene.save()

    .. image:: /images/10.0Å_1layer_graphene-1.png

    """
    @classmethod
    def generate_fname(cls, edge_length=None, **kwargs):
        """Generate a filename string."""
        dimensions = '{:.1f}Å'.format(edge_length)
        fname = '_'.join((dimensions, super().generate_fname(**kwargs)))
        return fname

    def save(self, fname=None, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save` method
        for documentation.

        """
        if fname is None:
            fname = \
                self.generate_fname(edge_length=self.edge_length,
                                    nlayers=self.nlayers, basis=self.basis)

        super().save(fname=fname, **kwargs)

HexagonalGrapheneGenerator = HexagonalCellGrapheneGenerator = \
    PrimitiveCellGrapheneGenerator


class ConventionalCellGrapheneGenerator(GrapheneGeneratorBase,
                                        ConventionalCellGraphene):
    """`N`-layer graphene generator class using a conventional unit cell.

    Parameters
    ----------
    armchair_edge_length : float, optional
        Length of armchair edge in **Angstroms**
    zigzag_edge_length : float, optional
        Length of zigzag edge in **Angstroms**
    basis : {:class:`python:list`}, optional
        List of :class:`python:str`\ s of element symbols or atomic number
        of the two atom basis (default: ['C', 'C'])
    bond : float, optional
        bond length between nearest-neighbor atoms in **Angstroms**.
    nlayers : int, optional
        Number of graphene layers.
    layer_spacing : float, optional
        Distance between layers in **Angstroms**.
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers
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
    autogen : bool, optional
        automatically generate unit cell and full structure
    verbose : bool, optional
        verbose output

    """
    @classmethod
    def generate_fname(cls, armchair_edge_length=None,
                       zigzag_edge_length=None, **kwargs):
        """Generate a filename string."""
        dimensions = '{:.1f}Åx{:.1f}Å'.format(armchair_edge_length,
                                              zigzag_edge_length)
        fname = '_'.join((dimensions, super().generate_fname(**kwargs)))
        return fname

    def save(self, fname=None, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save` method
        for documentation.

        """
        if fname is None:
            fname = \
                self.generate_fname(
                    armchair_edge_length=self.armchair_edge_length,
                    zigzag_edge_length=self.zigzag_edge_length,
                    nlayers=self.nlayers, basis=self.basis)

        super().save(fname=fname, **kwargs)

RectangularGrapheneGenerator = RectangularCellGrapheneGenerator = \
    ConventionalCellGrapheneGenerator


class GrapheneGenerator(ConventionalCellGrapheneGenerator):
    """`N`-layer graphene generator class.

    .. versionchanged:: 0.3.11

       `GrapheneGenerator` is now a sub-class of the
       `ConventionalCellGrapheneGenerator` class to maintain
       backwards compatibility and also includes 2 new
       classmethods: :meth:`~GrapheneGenerator.from_primitive_cell`
       and :meth:`~GrapheneGenerator.from_conventional_cell`.

    Parameters
    ----------
    armchair_edge_length : float, optional
        Length of armchair edge in **Angstroms**

        .. versionadded:: 0.3.10

    zigzag_edge_length : float, optional
        Length of zigzag edge in **Angstroms**

        .. versionadded:: 0.3.10

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
        bond length between nearest-neighbor atoms in **Angstroms**.
    nlayers : int, optional
        Number of graphene layers.
    layer_spacing : float, optional
        Distance between layers in **Angstroms**.
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers
    autogen : bool, optional
        automatically generate unit cell and full structure
    verbose : bool, optional
        verbose output

    Notes
    -----
    The `GrapheneGenerator` class and its subclasses generate graphene with
    either an armchair or zigzag edge using a 4-atom conventional unit cell.
    If you want to generate graphene as an *unrolled nanotube*, see the
    :class:`~sknano.generators.UnrolledSWNTGenerator` class.

    .. seealso:: :class:`~sknano.generators.UnrolledSWNTGenerator`

    Examples
    --------

    Start an interactive python or ipython session, then import the
    `GrapheneGenerator` class.

    >>> from sknano.generators import GrapheneGenerator

    Now generate a **100 Å AC x 10 Å ZZ** graphene nano-ribbon.

    >>> armchair_nanoribbon = GrapheneGenerator(armchair_edge_length=100,
    ...                                         zigzag_edge_length=10)

    Save structure data in default `xyz` format:

    >>> armchair_nanoribbon.save()

    The rendered structure look like:

    .. image:: /images/100.0Åx10.0Å_1layer_graphene.png

    Now let's generate a **10 Å ZZ x 100 Å AC** graphene nano-ribbon.

    >>> zigzag_nanoribbon = GrapheneGenerator(armchair_edge_length=10,
    ...                                       zigzag_edge_length=100)
    >>> zigzag_nanoribbon.save()

    The rendered structure looks like:

    .. image:: /images/10.0Åx100.0Å_1layer_graphene.png

    Now generate **100 Å AC x 25 Å ZZ**, 5 layer, `AB`-stacked graphene.

    >>> five_layer_graphene = GrapheneGenerator(armchair_edge_length=100,
    ...                                         zigzag_edge_length=25,
    ...                                         nlayers=5)
    >>> five_layer_graphene.save()

    The rendered structure looks like:

    .. image:: /images/100.0Åx25.0Å_5layer_graphene-1.png

    Now generate single layer, **10 Å x 10 Å** sheet of BN Graphene.

    >>> BN_graphene = GrapheneGenerator(armchair_edge_length=10,
    ...                                 zigzag_edge_length=10,
    ...                                 basis=['B', 'N'])
    >>> BN_graphene.save()

    The rendered structure looks like:

    .. image:: /images/10.0Åx10.0Å_1layer_graphene_B-N_basis-1.png

    """
    @classmethod
    def from_primitive_cell(cls, **kwargs):
        """`classmethod <https://docs.python.org/3/library/functions.html#classmethod>`_
        to call :class:`PrimitiveCellGrapheneGenerator`.

        """
        return PrimitiveCellGrapheneGenerator(**kwargs)

    @classmethod
    def from_conventional_cell(cls, **kwargs):
        """`classmethod <https://docs.python.org/3/library/functions.html#classmethod>`_
        to call :class:`ConventionalCellGrapheneGenerator`.

        """
        return ConventionalCellGrapheneGenerator(**kwargs)
