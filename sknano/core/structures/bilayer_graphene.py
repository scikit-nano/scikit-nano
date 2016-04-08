# -*- coding: utf-8 -*-
"""
===================================================================================
Bilayer Graphene structure class (:mod:`sknano.core.structures.bilayer_graphene`)
===================================================================================

.. currentmodule:: sknano.core.structures.bilayer_graphene

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from .graphene import Graphene

__all__ = ['BilayerGraphene']


class BilayerGraphene(Graphene):
    """Bilayer Graphene structure class.

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
    def __init__(self, **kwargs):
        super().__init__(nlayers=2, **kwargs)
