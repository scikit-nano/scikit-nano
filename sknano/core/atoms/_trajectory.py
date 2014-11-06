# -*- coding: utf-8 -*-
"""
===============================================================================
Trajectory class for MD Atoms analysis (:mod:`sknano.core.atoms._trajectory`)
===============================================================================

Class for analyzing the trajectories of molecular dynamics `Atoms`.

.. currentmodule:: sknano.core.atoms._trajectory

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#import numbers

#import numpy as np

#from sknano.core.math import Vector, vector as vec
#from ._bond import Bond
#from ._bonds import Bonds
#from ._extended_atoms import XAtoms
#from ._neighbor_atoms import NeighborAtoms

__all__ = ['Trajectory']


class Trajectory(object):
    """Base class for trajectory analysis."""

    def __init__(self, **kwargs):
        #super(Trajectory, self).__init__(**kwargs)
        pass
