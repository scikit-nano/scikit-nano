# -*- coding: utf-8 -*-
"""
====================================================================================
Atom classes for topological defects (:mod:`sknano.core.atoms.mixins._defect_atoms`)
====================================================================================

.. currentmodule:: sknano.core.atoms.mixins._defect_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import numbers

# import numpy as np

# from ._atoms import Atoms
# from ._bonds import Bonds, BondList

__all__ = ['DefectAtomMixin', 'DefectAtomsMixin']


class DefectAtomMixin:
    """An `Atom` class for defect analysis."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._holes = []

    @property
    def holes(self):
        return self._holes


class DefectAtomsMixin:
    """Atoms class for defect analysis."""
    pass
