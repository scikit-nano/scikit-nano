# -*- coding: utf-8 -*-
"""
==============================================================================
Graphene structure generator (:mod:`sknano.generators._graphene_generator`)
==============================================================================

.. currentmodule:: sknano.generators._graphene_generator

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core.math import Vector
from sknano.structures import Graphene
from ._base import Atom, Atoms, GeneratorBase

__all__ = ['GrapheneGenerator']


class GrapheneGenerator(Graphene, GeneratorBase):
    """`N`-layer graphene generator class.

    Parameters
    ----------
    armchair_edge_length : float, optional
        Length of armchair edge in **nanometers**

        .. versionadded:: 0.3.10

    zigzag_edge_length : float, optional
        Length of zigzag edge in **nanometers**

        .. versionadded:: 0.3.10

    length : float, optional
        Length of armchair edge in **nanometers**

        .. deprecated:: 0.3.10
           Use `armchair_edge_length` instead

    width : float, optional
        Width of graphene sheet in **nanometers**

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

    Now generate a **20 nm AC x 1 nm ZZ** graphene nano-ribbon.

    >>> ACG = GrapheneGenerator(armchair_edge_length=20, zigzag_edge_length=1)

    Save structure data in `xyz` format:

    >>> ACG.save()

    The rendered structure look like:

    .. image:: /images/20nmx1nm_AC_edge.png

    Now let's generate a **20 nm ZZ x 1 nm AC** graphene nano-ribbon.

    >>> ZZG = GrapheneGenerator(armchair_edge_length=20, zigzag_edge_length=1)
    >>> ZZG.save()

    The rendered structure looks like:

    .. image:: /images/20nmx1nm_ZZ_edge.png

    Now generate **25 nm AC x 5 nm ZZ**, 5 layer, `AB`-stacked graphene.

    >>> ACG_5layers = GrapheneGenerator(armchair_edge_length=25,
    ...                                 zigzag_edge_length=5,
    ...                                 nlayers=5)
    >>> ACG_5layers.save()

    The rendered structure looks like:

    .. image:: /images/25nmx5nm_5layer_AC_graphene.png

    Now generate single layer, **10 nm x 10 nm** sheet of BN Graphene.

    >>> BN_graphene = GrapheneGenerator(armchair_edge_length=10,
    ...                                 zigzag_edge_length=10,
    ...                                 basis=['B', 'N'])
    >>> BN_graphene.save()

    The rendered structure looks like:

    .. image:: /images/10nmx10nm_single_layer_BN_graphene.png

    Now, just because we can, generate a **5 nm x 5 nm** sheet of
    Uranium-Einsteinium Graphene.

    >>> UEs_graphene = GrapheneGenerator(armchair_edge_length=5,
    ...                                  zigzag_edge_length=5,
    ...                                  basis=['U', 'Es'])
    >>> UEs_graphene.save()

    The rendered structure looks like:

    .. image:: /images/5nmx5nm_single_layer_UEs_graphene.png

    """

    def __init__(self, autogen=True, **kwargs):

        super().__init__(**kwargs)

        if autogen:
            self.generate_structure_data()

    def generate_structure_data(self):
        """Generate the full structure coordinates."""

        self.structure_data.clear()
        for nlayer in range(self.nlayers):
            layer = Atoms()
            for nx in range(self.Nx):
                for ny in range(self.Ny):
                    dr = Vector([nx * self.unit_cell.a,
                                 ny * self.unit_cell.b,
                                 0.0])
                    for atom in self.unit_cell:
                        layer_atom = Atom(atom.symbol)
                        layer_atom.r = atom.r + dr
                        layer.append(layer_atom)

            if self.Natoms_per_layer is None:
                self.Natoms_per_layer = layer.Natoms

            # translate layer to put its center of mass at the origin
            layer.center_CM()
            layer.translate(Vector([0, 0, nlayer * self.layer_spacing]))
            if (nlayer % 2) != 0:
                layer.translate(self.layer_shift)
            layer.rotate(angle=self.layer_rotation_angles[nlayer], axis='z')
            self.atoms.extend(layer)

    @classmethod
    def generate_fname(cls, armchair_edge_length, zigzag_edge_length,
                       nlayers, basis):
        dimensions = '{}nmx{}nm'.format(armchair_edge_length,
                                        zigzag_edge_length)
        nlayer = '{}layer'.format(nlayers)
        basis = ''.join(basis)
        fname_wordlist = (dimensions, nlayer, basis, 'graphene')
        fname = '_'.join(fname_wordlist)
        return fname

    def save(self, fname=None, outpath=None, structure_format=None,
             center_CM=True, rotation_angle=-np.pi/2, rotation_axis='x',
             **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save` method
        for documentation.

        """
        if fname is None:
            fname = self.generate_fname(self.armchair_edge_length,
                                        self.zigzag_edge_length,
                                        self.nlayers, self.basis)

        if center_CM and self.nlayers > 1:
            self.atoms.center_CM()

        super().save(fname=fname, outpath=outpath,
                     structure_format=structure_format, center_CM=False,
                     angle=rotation_angle, axis=rotation_axis, **kwargs)
