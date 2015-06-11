# -*- coding: utf-8 -*-
"""
==============================================================================
MWNT bundle class (:mod:`sknano.structures._mwnt_bundle`)
==============================================================================

.. currentmodule:: sknano.structures._mwnt_bundle

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._nanotube_bundle import NanotubeBundleBase
from ._mwnt import MWNT

__all__ = ['MWNTBundle']


class MWNTBundle(NanotubeBundleBase, MWNT):
    """MWNT bundle structure class."""
    pass
