# -*- coding: utf-8 -*-
"""
==============================================================================
Graphene structure generator (:mod:`sknano.generators._graphene_generator`)
==============================================================================

.. currentmodule:: sknano.generators._graphene_generator

"""
from __future__ import absolute_import, division, print_function
from six.moves import range
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.structures import Graphene
from ._base import Atom, Atoms, GeneratorBase

__all__ = ['GrapheneGenerator']


class GrapheneGenerator(Graphene, GeneratorBase):
    """`N`-layer graphene generator class.

    Parameters
    ----------
    length : float
        Length of graphene sheet in **nanometers**
    width : float
        Width of graphene sheet in **nanometers**
    edge : {'AC', 'armchair', 'ZZ', 'zigzag'}, optional
        **A**\ rm\ **C**\ hair or **Z**\ ig\ **Z**\ ag edge along
        the `length` of the sheet.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2
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

    Now generate a **20 nm x 1 nm** armchair edge graphene nano-ribbon.

    >>> ACG = GrapheneGenerator(length=20, width=1, edge='AC')

    Save structure data in `xyz` format:

    >>> ACG.save_data()

    The rendered structure look like:

    .. image:: /images/20nmx1nm_AC_edge.png

    Now let's generate a **20 nm x 1 nm** zigzag edge graphene nano-ribbon.

    >>> ZZG = GrapheneGenerator(length=20, width=1, edge='ZZ')
    >>> ZZG.save_data()

    The rendered structure looks like:

    .. image:: /images/20nmx1nm_ZZ_edge.png

    Now generate **25 nm x 5 nm**, `armchair` edge,
    5 layer, `AB`-stacked graphene.

    >>> ACG_5layers = GrapheneGenerator(length=25, width=5, edge='AC',
    ...                                 nlayers=5)
    >>> ACG_5layers.save_data()

    The rendered structure looks like:

    .. image:: /images/25nmx5nm_5layer_AC_graphene.png

    Now generate single layer, **10 nm x 10 nm** sheet of BN Graphene.

    >>> BN_graphene = GrapheneGenerator(length=10, width=10, edge='AC',
    ...                                 element1='B', element2='N')
    >>> BN_graphene.save_data()

    The rendered structure looks like:

    .. image:: /images/10nmx10nm_single_layer_BN_graphene.png

    Now, just because we can, generate a **5 nm x 5 nm** sheet of
    Uranium-Einsteinium Graphene.

    >>> UEs_graphene = GrapheneGenerator(width=5, length=5, edge='zigzag',
    ...                                  element1='U', element2='Es')
    >>> UEs_graphene.save_data()

    The rendered structure looks like:

    .. image:: /images/5nmx5nm_single_layer_UEs_graphene.png

    """

    def __init__(self, autogen=True, **kwargs):

        super(GrapheneGenerator, self).__init__(**kwargs)

        if autogen:
            self.generate_unit_cell()
            self.generate_structure_data()

    def generate_unit_cell(self):
        """Generate the graphene unit cell.

        Called automatically if `autogen` is `True`.

        """
        bond = self.bond
        edge = self.edge
        e1 = self.element1
        e2 = self.element2

        self.unit_cell = Atoms()
        # Set up 4 atom basis
        # Leave atom 1 at the origin
        atom1 = Atom(element=e1)
        atom2 = Atom(element=e2)
        atom3 = Atom(element=e1)
        atom4 = Atom(element=e2)
        if edge == 'AC':
            # Move atom 2 to 2nd basis position
            atom2.x = -np.sqrt(3) / 2 * bond
            atom2.y = bond / 2

            # Move atom 3 from atom 1 along a1 primitive vector
            atom3.x = -np.sqrt(3) / 2 * bond
            atom3.y = 3 / 2 * bond

            # Move atom 4 from atom 2 along a2 primitive vector
            atom4.y = 2 * bond
        else:
            # Move atom 2 to 2nd basis position
            atom2.x = bond

            # Move atom 3 from atom 1 along a1 primitive vector
            atom3.x = 3 / 2 * bond
            atom3.y = np.sqrt(3) / 2 * bond

            # Move atom 4 from atom 2 along a2 primitive vector
            atom4.x = -bond / 2
            atom4.y = np.sqrt(3) / 2 * bond

        self.unit_cell.extend([atom1, atom2, atom3, atom4])
        #self.Natoms = self.unit_cell.Natoms

    def generate_structure_data(self):
        """Generate the full structure coordinates."""
        #self.atoms = Atoms()
        self.structure_data.clear()
        for nlayer in range(self.nlayers):
            layer = Atoms()
            for nx in range(self.Nx):
                for ny in range(self.Ny):
                    dr = np.array([nx * self.cell.x,
                                   ny * self.cell.y,
                                   nlayer * self.layer_spacing])
                    for atom in self.unit_cell:
                        layer_atom = Atom(element=atom.symbol)
                        layer_atom.r = atom.r + dr
                        layer.append(layer_atom)

            if self.Natoms_per_layer is None:
                self.Natoms_per_layer = layer.Natoms

            # translate layer to put its center of mass at the origin
            dr = layer.CM
            dr.z = 0
            layer.translate(dr)
            if (nlayer % 2) != 0:
                layer.translate(self.layer_shift)

            self.atoms.extend(layer)

    def save_data(self, fname=None, outpath=None, structure_format=None,
                  rotation_angle=-90, rot_axis='x', deg2rad=True,
                  center_CM=True, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save_data` method
        for documentation.

        """
        if fname is None:
            dimensions = '{}nmx{}nm'.format(self.length, self.width)
            nlayer = '{}layer'.format(self.nlayers)
            edge = 'AC' if self.edge in ('AC', 'armchair') else 'ZZ'
            atombond = '{}{}'.format(self.element1, self.element2)
            fname_wordlist = (dimensions, nlayer, edge, atombond, 'graphene')
            fname = '_'.join(fname_wordlist)

        if center_CM and self.nlayers > 1:
            self.atoms.center_CM()

        super(GrapheneGenerator, self).save_data(
            fname=fname, outpath=outpath, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            deg2rad=deg2rad, center_CM=False, **kwargs)
