# -*- coding: utf-8 -*-
"""
===============================================================================
Atoms with several attributes (:mod:`sknano.core.atoms._extended_atoms`)
===============================================================================

An "eXtended" `Atoms` class.

.. currentmodule:: sknano.core.atoms._extended_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

from ._xyz_atoms import XYZAtom, XYZAtoms
from ._velocity_atoms import VelocityAtom, VelocityAtoms
from ._force_atoms import ForceAtom, ForceAtoms
from ._cn_atoms import CNAtom, CNAtoms
from ._energy_atoms import EnergyAtom, EnergyAtoms
from ._charged_atoms import ChargedAtom, ChargedAtoms
from ._id_atoms import IDAtom, IDAtoms
from ._image_atoms import ImageAtom, ImageAtoms
from ._type_atoms import TypeAtom, TypeAtoms

__all__ = ['XAtom', 'XAtoms']


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


class XAtoms(IDAtoms, TypeAtoms, XYZAtoms, ImageAtoms, ChargedAtoms,
             VelocityAtoms, ForceAtoms, EnergyAtoms, CNAtoms):
    """An eXtended `Atoms` class.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.XAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `XAtoms`}, optional
        if not `None`, then a list of `XAtom` instance objects or an
        existing `XAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return XAtom

    def sort(self, key=attrgetter('element', 'Z', 'mass', 'id', 'mol', 'type',
                                  'x', 'y', 'z', 'q', 'etotal',
                                  'CN'), reverse=False):
        super().sort(key=key, reverse=reverse)
