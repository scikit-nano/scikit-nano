# -*- coding: utf-8 -*-
"""
==============================================================================
MWNT Bundle structure class (:mod:`sknano.structures._mwnt_bundle`)
==============================================================================

.. currentmodule:: sknano.structures._mwnt_bundle

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._nanotube_bundle import NanotubeBundle
from ._mwnt import MWNT

__all__ = ['MWNTBundle']


class MWNTBundle(NanotubeBundle, MWNT):
    def __init__(self, **kwargs):

        super(MWNTBundle, self).__init__(**kwargs)
