# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class for NN analysis (:mod:`sknano.core.atoms._neighbor_atom`)
===============================================================================

.. currentmodule:: sknano.core.atoms._neighbor_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import sknano.core.atoms
# from ._cn_atom import CNAtom
# from ._id_atom import IDAtom
# from ._xyz_atom import XYZAtom
# from ._bond import Bond
# from ._bonds import Bonds
from ._kdtree_atom import KDTAtom

__all__ = ['NeighborAtom']


NeighborAtom = KDTAtom
