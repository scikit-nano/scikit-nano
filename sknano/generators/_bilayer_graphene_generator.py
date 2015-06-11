# -*- coding: utf-8 -*-
"""
===============================================================================
BLG structure generator (:mod:`sknano.generators._bilayer_graphene_generator`)
===============================================================================

.. currentmodule:: sknano.generators._bilayer_graphene_generator

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from sknano.structures import BilayerGraphene
from ._graphene_generator import GrapheneGenerator

__all__ = ['BilayerGrapheneGenerator']


class BilayerGrapheneGenerator(BilayerGraphene, GrapheneGenerator):
    """Bilayer graphene structure generator class.

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
    layer_spacing : float, optional
        Distance between layers in **Angstroms**.
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers
    layer_rotation_angle : {None, float}, optional
        Rotation angle of second layer specified in degrees.
        If specified in degrees, then you must set `degrees=True`
    degrees : bool, optional
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
    bilayer-graphene:

    >>> blg = BilayerGrapheneGenerator(armchair_edge_length=10,
    ...                                zigzag_edge_length=1)

    Save structure data in `xyz` format:

    >>> blg.save()

    The rendered structure looks like (after rotating 90 degrees so that
    it better fits the page):

    .. image:: /images/10nmx1nm_bilayer.png

    Now generate bilayer-graphene with top layer rotated by 45 degrees.

    >>> rotated_bilayer = BilayerGrapheneGenerator(armchair_edge_length=10,
    ...                                            zigzag_edge_length=10,
    ...                                            layer_rotation_angle=45,
    ...                                            degrees=True)
    >>> rotated_bilayer.save(fname='rotated_bilayer.xyz')

    The rendered structure looks like:

    .. image:: /images/rotated_bilayer.png

    Now generate BN bilayer-graphene with top layer rotated 45 degrees.

    >>> rotated_BN_bilayer = BilayerGrapheneGenerator(armchair_edge_length=10,
    ...                                               zigzag_edge_length=10,
    ...                                               basis=['B', 'N'],
    ...                                               layer_rotation_angle=45,
    ...                                               degrees=True)
    >>> rotated_BN_bilayer.save(fname='BN_bilayer_rotated_45deg.xyz')

    The rendered structure looks like:

    .. image:: /images/BN_bilayer_rotated_45deg.png

    """
    pass
