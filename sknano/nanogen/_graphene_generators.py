# -*- coding: utf-8 -*-
"""
==========================================================================
Graphene structure generators (:mod:`sknano.nanogen._graphene_generators`)
==========================================================================

.. currentmodule:: sknano.nanogen._graphene_generators

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import copy

import numpy as np

from ..chemistry import Atom, Atoms
from ..tools import rotation_matrix
from ..tools.refdata import CCbond

from ._graphene import Graphene
from ._structure_generator import StructureGenerator

__all__ = ['GrapheneGenerator', 'BiLayerGrapheneGenerator']


class GrapheneGenerator(Graphene, StructureGenerator):
    """Class for generating `N`-layer graphene nanostructures.

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
        :class:`~sknano.chemistry.Atoms` 1 and 2
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
    :class:`~sknano.nanogen.UnrolledNanotubeGenerator` class.

    .. seealso:: :class:`~sknano.nanogen.UnrolledNanotubeGenerator`

    Examples
    --------

    Start an interactive python or ipython session, then import the
    `GrapheneGenerator` class.

    >>> from sknano.nanogen import GrapheneGenerator

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

    def __init__(self, length=None, width=None, edge=None,
                 element1='C', element2='C', bond=CCbond,
                 nlayers=1, layer_spacing=3.35, stacking_order='AB',
                 autogen=True, with_units=False, units=None,
                 verbose=False):

        super(GrapheneGenerator, self).__init__(
            length=length, width=width, edge=edge,
            element1=element1, element2=element2, bond=bond,
            nlayers=nlayers, layer_spacing=layer_spacing,
            stacking_order=stacking_order, with_units=with_units, units=units,
            verbose=verbose)

        if autogen:
            self.generate_unit_cell()
            self.generate_structure_data()

    def generate_unit_cell(self):
        """Generate the graphene unit cell.

        Called automatically if `autogen` is `True`.

        """
        bond = self._bond
        edge = self._edge
        e1 = self._element1
        e2 = self._element2

        self._unit_cell = Atoms()
        # Set up 4 atom basis
        # Leave atom 1 at the origin
        atom1, atom2, atom3, atom4 = Atom(e1), Atom(e2), Atom(e1), Atom(e2)
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

        for atom in (atom1, atom2, atom3, atom4):
            self._unit_cell.append(atom)
        self._Natoms = self._unit_cell.Natoms

    def generate_structure_data(self):
        """Generate the full structure coordinates."""
        self._structure_atoms = Atoms()
        for nlayer in xrange(self._nlayers):
            layer = Atoms()
            for nx in xrange(self._Nx):
                for ny in xrange(self._Ny):
                    dr = np.array([nx * self._cell.x,
                                   ny * self._cell.y,
                                   nlayer * self._layer_spacing])
                    for atom in self._unit_cell:
                        layer_atom = Atom(atom.symbol)
                        layer_atom.r = atom.r + dr
                        layer.append(layer_atom)

            if self._Natoms_per_layer is None:
                self._Natoms_per_layer = layer.Natoms

            # translate layer to put its center of mass at the origin
            layer.center_CM(r_indices=[0, 1])
            if (nlayer % 2) != 0:
                layer.translate(self._layer_shift.components)

            self._structure_atoms.extend(layer.atoms)

    def save_data(self, fname=None, outpath=None, structure_format=None,
                  rotation_angle=-90, rot_axis='x', deg2rad=True,
                  center_CM=True):
        """Save structure data.

        See :meth:`~sknano.nanogen.StructureGenerator.save_data` method
        for documentation.

        """
        if fname is None:
            dimensions = '{}nmx{}nm'.format(self._length, self._width)
            nlayer = '{}layer'.format(self._nlayers)
            edge = 'AC' if self._edge in ('AC', 'armchair') else 'ZZ'
            atombond = '{}{}'.format(self._element1, self._element2)
            fname_wordlist = (dimensions, nlayer, edge, atombond, 'graphene')
            fname = '_'.join(fname_wordlist)

        #structure_atoms = list(itertools.chain(*self._structure_atoms))
        #structure_atoms = Atoms(structure_atoms)
        if center_CM and self._nlayers > 1:
            self._structure_atoms.center_CM()

        super(GrapheneGenerator, self).save_data(
            fname=fname, outpath=outpath, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            deg2rad=deg2rad, center_CM=False)


class BiLayerGrapheneGenerator(GrapheneGenerator):
    """Class for generating bi-layer graphene structures.

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
        :class:`~sknano.chemistry.Atoms` 1 and 2
    bond : float, optional
        bond length between nearest-neighbor atoms in **Angstroms**.
    layer_spacing : float, optional
        Distance between layers in **Angstroms**.
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers
    rotation_angle : {None, float}, optional
        Rotation angle of second layer specified in degrees.
        If specified in radians, then you must set `deg2rad=False`
    deg2rad : bool, optional
        The rotation angle is specified in degrees and needs to be converted
        to radians.
    autogen : bool, optional
        if `True`, automatically generate unit cell and full structure
    verbose : bool, optional
        verbose output

    Examples
    --------

    Import the BiLayerGrapheneGenerator class

    >>> from sknano.nanogen import BiLayerGrapheneGenerator

    Generate **10 nm** wide by **1 nm** long `AB` stacked
    bilayer-graphene with a `ZZ` edge:

    >>> bi_graphene = BiLayerGrapheneGenerator(length=10, width=1, edge='ZZ')

    Save structure data in `xyz` format:

    >>> bi_graphene.save_data()

    The rendered structure looks like (after rotating 90 degrees so that
    it better fits the page):

    .. image:: /images/10nmx1nm_bilayer.png

    Now generate bilayer-graphene with top layer rotated by 45 degrees.

    >>> rotated_bilayer = BiLayerGrapheneGenerator(length=10, width=10,
    ...                                            edge='armchair',
    ...                                            rotation_angle=45)
    >>> rotated_bilayer.save_data(fname='rotated_bilayer.xyz')

    The rendered structure looks like:

    .. image:: /images/rotated_bilayer.png

    Now generate BN bilayer-graphene with top layer rotated 45 degrees.

    >>> rotated_BN_bilayer = BiLayerGrapheneGenerator(length=10, width=10,
    ...                                               edge='zigzag',
    ...                                               element1='B',
    ...                                               element2='N',
    ...                                               rotation_angle=45)
    >>> rotated_BN_bilayer.save_data(fname='BN_bilayer_rotated_45deg.xyz')

    The rendered structure looks like:

    .. image:: /images/BN_bilayer_rotated_45deg.png

    """

    def __init__(self, length=None, width=None, edge=None,
                 element1='C', element2='C', bond=CCbond,
                 layer_spacing=3.35, stacking_order='AB',
                 rotation_angle=None, deg2rad=True, autogen=True,
                 verbose=False):

        super(BiLayerGrapheneGenerator, self).__init__(
            length=length, width=width, edge=edge,
            element1=element1, element2=element2, bond=bond,
            nlayers=2, layer_spacing=layer_spacing, autogen=False,
            verbose=verbose)

        self._rotation_matrix = None
        if rotation_angle is not None:
            self._rotation_matrix = rotation_matrix(rotation_angle,
                                                    rot_axis='z',
                                                    deg2rad=deg2rad)

        if autogen:
            super(BiLayerGrapheneGenerator, self).generate_unit_cell()
            self.generate_structure_data()

    def generate_structure_data(self):
        """Generate the full structure coordinates."""
        super(BiLayerGrapheneGenerator, self).generate_structure_data()

        if self._rotation_matrix is not None:
            bilayer = copy.deepcopy(self._structure_atoms)
            self._structure_atoms = Atoms()

            z_coords = bilayer.get_coords(as_dict=True)['z']
            z_set = np.asarray(sorted(list(set(z_coords))))
            epsilon = 1e-10

            for n, z in enumerate(z_set):
                layer = Atoms(atoms=bilayer.get_atoms(asarray=True)[
                    np.where(np.abs(z_coords - z) < epsilon)].tolist(),
                    deepcopy=True)
                if (n % 2) != 0:
                    layer.rotate(self._rotation_matrix)

                self._structure_atoms.extend(layer.atoms)
