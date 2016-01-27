# -*- coding: utf-8 -*-
"""
===============================================================================
MoS2 structure class (:mod:`sknano.structures._molybdenum_disulphide`)
===============================================================================

.. currentmodule:: sknano.structures._molybdenum_disulphide

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

#import itertools
import numbers

import numpy as np

#from sknano.core.atoms import Atom
from sknano.core.math import Vector
from ._base import NanoStructureBase

__all__ = ['MoS2']


class MoS2(NanoStructureBase):
    """Molybdenum disulphide structure class.

    Parameters
    ----------

    """

    def __init__(self, nlayers=1, layer_spacing=3.14,
                 sulphur_layer_spacing=3.01, **kwargs):

        super().__init__(**kwargs)

        self.layer_mass = None
        self.Natoms = 0
        self.Natoms_per_layer = 0

        self.nlayers = nlayers
        self.layer_spacing = layer_spacing

        self.layer_shift = Vector()

        self.cell = Vector()
