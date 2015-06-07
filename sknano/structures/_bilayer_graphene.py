# -*- coding: utf-8 -*-
"""
===============================================================================
Bilayer Graphene structure class (:mod:`sknano.structures._bilayer_graphene`)
===============================================================================

.. currentmodule:: sknano.structures._bilayer_graphene

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._graphene import Graphene

__all__ = ['BilayerGraphene']


class BilayerGraphene(Graphene):
    """Bilayer Graphene structure class."""
    def __init__(self, **kwargs):

        super().__init__(nlayers=2, **kwargs)
