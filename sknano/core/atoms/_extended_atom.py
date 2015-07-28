# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class with extended feature set (:mod:`sknano.core.atoms._extended_atom`)
===============================================================================

An "eXtended" `Atom` class for structure analysis.

.. currentmodule:: sknano.core.atoms._extended_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from ._xyz_atom import XYZAtom
from ._velocity_atom import VelocityAtom
from ._force_atom import ForceAtom
from ._cn_atom import CNAtom
from ._energy_atom import EnergyAtom
from ._charged_atom import ChargedAtom
from ._id_atom import IDAtom
from ._image_atom import ImageAtom
from ._type_atom import TypeAtom

__all__ = ['XAtom']


class XAtom(IDAtom, TypeAtom, XYZAtom, ImageAtom, ChargedAtom, VelocityAtom,
            ForceAtom, EnergyAtom, CNAtom):
    """An `Atom` class with an eXtended set of attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    x, y, z : float, optional
        :math:`x, y, z` components of `XAtom` position vector relative to
        origin.
    id : int, optional
        atom ID
    mol : int, optional
        molecule ID
    type : int, optional
        atom type
    q : {int, float}, optional
        Net charge of `XAtom`.
    vx, vy, vz : float, optional
        :math:`v_x, v_y, v_z` components of `XAtom` velocity.

    """
    pass
