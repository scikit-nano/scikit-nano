# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class for structure analysis (:mod:`sknano.core.atoms._structure_atom`)
===============================================================================

.. currentmodule:: sknano.core.atoms._structure_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from ._charged_atom import ChargedAtom
from ._image_atom import ImageAtom
from ._type_atom import TypeAtom
from ._velocity_atom import VelocityAtom
from ._poav_atom import POAVAtom

__all__ = ['StructureAtom']


class StructureAtom(POAVAtom, ChargedAtom, VelocityAtom, ImageAtom, TypeAtom):
    """An `Atom` class for structure analysis.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    id : int, optional
        atom ID
    mol : int, optional
        molecule ID
    type : int, optional
        atom type
    x, y, z : float, optional
        :math:`x, y, z` components of `StructureAtom` position vector relative to
        origin.
    vx, vy, vz : float, optional
        :math:`v_x, v_y, v_z` components of `StructureAtom` velocity.

    """
    pass
