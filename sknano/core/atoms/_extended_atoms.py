# -*- coding: utf-8 -*-
"""
===============================================================================
Extended Atoms class feature set (:mod:`sknano.core.atoms._extended_atoms`)
===============================================================================

An "eXtended" `Atoms` class.

.. currentmodule:: sknano.core.atoms._extended_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

from ._extended_atom import XAtom
from ._xyz_atoms import XYZAtoms
from ._velocity_atoms import VelocityAtoms
from ._force_atoms import ForceAtoms
from ._cn_atoms import CNAtoms
from ._energy_atoms import EnergyAtoms
from ._charged_atoms import ChargedAtoms
from ._id_atoms import IDAtoms
from ._image_atoms import ImageAtoms
from ._type_atoms import TypeAtoms

__all__ = ['XAtoms']


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
