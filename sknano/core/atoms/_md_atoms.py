# -*- coding: utf-8 -*-
"""
===============================================================================
Atoms class for MD analysis (:mod:`sknano.core.atoms._md_atoms`)
===============================================================================

An `Atoms` class for molecular dynamics structure analysis.

.. currentmodule:: sknano.core.atoms._md_atoms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#import numbers
#import numpy as np

#from sknano.core.math import Vector, vector as vec
#from ._bond import Bond
#from ._bonds import Bonds
#from ._extended_atoms import XAtoms
from ._kdtree_atoms import KDTAtoms

__all__ = ['MDAtoms']


class MDAtoms(KDTAtoms):
    """An `Atoms` class for molecular dynamics structure analysis."""
    def __init__(self, **kwargs):
        super(MDAtoms, self).__init__(**kwargs)
