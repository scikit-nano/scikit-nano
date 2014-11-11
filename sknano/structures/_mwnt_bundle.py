# -*- coding: utf-8 -*-
"""
==============================================================================
MWNT bundle structure class (:mod:`sknano.structures._mwnt_bundle`)
==============================================================================

.. currentmodule:: sknano.structures._mwnt_bundle

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._nanotube_bundle import NanotubeBundleBase
from ._mwnt import MWNT
#from ._swnt_bundle import SWNTBundle

__all__ = ['MWNTBundle']


class MWNTBundle(NanotubeBundleBase, MWNT):
    """MWNT bundle structure class."""
    def __init__(self, **kwargs):

        super(MWNTBundle, self).__init__(**kwargs)

        #self.shell_bundles = [SWNTBundle(**swnt.todict()
