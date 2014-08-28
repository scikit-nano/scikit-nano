# -*- coding: utf-8 -*-
"""
===============================================================================
Graphene structure classes (:mod:`sknano.structures._bilayer_graphene`)
===============================================================================

.. currentmodule:: sknano.structures._bilayer_graphene

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#import itertools

import numpy as np

from ._graphene import Graphene

edge_types = {'armchair': 'AC', 'zigzag': 'ZZ'}

__all__ = ['BilayerGraphene']


class BilayerGraphene(Graphene):
    def __init__(self, layer_rotation_angle=None, deg2rad=True, **kwargs):

        if layer_rotation_angle is not None and deg2rad:
            layer_rotation_angle = np.radians(layer_rotation_angle)
        kwargs['layer_rotation_angles'] = layer_rotation_angle
        kwargs['nlayers'] = 2

        super(BilayerGraphene, self).__init__(**kwargs)
        self._layer_rotation_angle = self._layer_rotation_angles

    @property
    def layer_rotation_angle(self):
        return self._layer_rotation_angle

    @layer_rotation_angle.setter
    def layer_rotation_angle(self, value):
        self._layer_rotation_angle = value
