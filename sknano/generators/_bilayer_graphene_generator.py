# -*- coding: utf-8 -*-
"""
===============================================================================
Bilayer Graphene (:mod:`sknano.generators._bilayer_graphene_generator`)
===============================================================================

.. currentmodule:: sknano.generators._bilayer_graphene_generator

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import copy

import numpy as np

from sknano.structures import BilayerGraphene
from ._base import Atoms, GeneratorBase

__all__ = ['BilayerGrapheneGenerator']


class BilayerGrapheneGenerator(BilayerGraphene, GeneratorBase):
    """Bilayer graphene structure generator class.

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
    layer_spacing : float, optional
        Distance between layers in **Angstroms**.
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers
    layer_rotation_angle : {None, float}, optional
        Rotation angle of second layer specified in degrees.
        If specified in radians, then you must set `deg2rad=False`
    deg2rad : bool, optional
        The `layer_rotation_angle` is specified in degrees and needs to be
        converted to radians.
    autogen : bool, optional
        if `True`, automatically generate unit cell and full structure
    verbose : bool, optional
        verbose output

    Examples
    --------

    Import the BilayerGrapheneGenerator class

    >>> from sknano.generators import BilayerGrapheneGenerator

    Generate **10 nm** wide by **1 nm** long `AB` stacked
    bilayer-graphene with a `ZZ` edge:

    >>> bi_graphene = BilayerGrapheneGenerator(length=10, width=1, edge='ZZ')

    Save structure data in `xyz` format:

    >>> bi_graphene.save_data()

    The rendered structure looks like (after rotating 90 degrees so that
    it better fits the page):

    .. image:: /images/10nmx1nm_bilayer.png

    Now generate bilayer-graphene with top layer rotated by 45 degrees.

    >>> rotated_bilayer = BilayerGrapheneGenerator(length=10, width=10,
    ...                                            edge='armchair',
    ...                                            layer_rotation_angle=45)
    >>> rotated_bilayer.save_data(fname='rotated_bilayer.xyz')

    The rendered structure looks like:

    .. image:: /images/rotated_bilayer.png

    Now generate BN bilayer-graphene with top layer rotated 45 degrees.

    >>> rotated_BN_bilayer = BilayerGrapheneGenerator(length=10, width=10,
    ...                                               edge='zigzag',
    ...                                               element1='B',
    ...                                               element2='N',
    ...                                               layer_rotation_angle=45)
    >>> rotated_BN_bilayer.save_data(fname='BN_bilayer_rotated_45deg.xyz')

    The rendered structure looks like:

    .. image:: /images/BN_bilayer_rotated_45deg.png

    """

    def __init__(self, autogen=True, **kwargs):

        super(BilayerGrapheneGenerator, self).__init__(**kwargs)

        if autogen:
            super(BilayerGrapheneGenerator, self).generate_unit_cell()
            self.generate_structure_data()

    def generate_structure_data(self):
        """Generate the full structure coordinates."""
        super(BilayerGrapheneGenerator, self).generate_structure_data()

        if self.layer_rotation_angle is not None:
            bilayer = copy.deepcopy(self.atoms)
            self.atoms = Atoms()

            z_coords = bilayer.get_coords(asdict=True)['z']
            z_set = np.asarray(sorted(list(set(z_coords))))
            epsilon = 1e-10

            for n, z in enumerate(z_set):
                layer = Atoms(atoms=bilayer.get_atoms(asarray=True)[
                    np.where(np.abs(z_coords - z) < epsilon)].tolist(),
                    deepcopy=True)
                if (n % 2) != 0:
                    layer.rotate(angle=self.layer_rotation_angle, rot_axis='z')

                self.atoms.extend(layer.atoms)
