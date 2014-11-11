# -*- coding: utf-8 -*-
"""
===============================================================================
StructureAtom(s) class references (:mod:`sknano.core.atoms._structure_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._structure_atoms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._poav_atom import POAVAtom as StructureAtom
from ._poav_atoms import POAVAtoms as StructureAtoms

__all__ = ['StructureAtom', 'StructureAtoms']
