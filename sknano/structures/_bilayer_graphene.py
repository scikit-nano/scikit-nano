# -*- coding: utf-8 -*-
"""
===============================================================================
Graphene structure classes (:mod:`sknano.structures._bilayer_graphene`)
===============================================================================

.. currentmodule:: sknano.structures._bilayer_graphene

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._graphene import Graphene

__all__ = ['BilayerGraphene']


class BilayerGraphene(Graphene):
    """Bilayer Graphene structure class."""
    def __init__(self, layer_rotation_angle=None, **kwargs):

        super(BilayerGraphene, self).__init__(
            nlayers=2, layer_rotation_angles=layer_rotation_angle, **kwargs)

        self.layer_rotation_angle = self.layer_rotation_angles
