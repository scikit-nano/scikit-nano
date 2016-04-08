# -*- coding: utf-8 -*-
"""
====================================================================================
Atom classes for topological defects (:mod:`sknano.core.atoms.mixins.defect_atoms`)
====================================================================================

.. currentmodule:: sknano.core.atoms.mixins.defect_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import Counter, OrderedDict
# import numbers

# import numpy as np

# from .atoms import Atoms
# from .bonds import Bonds, BondList

__all__ = ['DefectAtomMixin', 'DefectAtomsMixin']


class DefectAtomMixin:
    """An `Atom` class for defect analysis."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.defects = OrderedDict()
        self.defect_counter = Counter()


class DefectAtomsMixin:
    """Atoms class for defect analysis."""

    def __init__(self, *args, **kwargs):
        self.defects = OrderedDict()
        self.defect_counter = Counter()

    def analyze_defects(self, defect_type=None, defect_condition=None,
                        **kwargs):
        """Analyze defects."""
        if self.verbose:
            print('Analyzing defects...')
        atoms = self.atoms
        rings, ring_cntr = \
            atoms.analyze_network(cutoff=kwargs.get('cutoff', 1.5),
                                  max_ring_size=kwargs.get('max_ring_size', 20),
                                  retcodes=('rings', 'rings_counter'))
        # for n, list_of_atoms in rings.items():
        #     for ring_atoms in list_of_atoms:
