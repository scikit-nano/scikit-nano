# -*- coding: utf-8 -*-
"""
==============================================================================
Unrolled SWNT structure class (:mod:`sknano.structures._unrolled_swnt`)
==============================================================================

.. currentmodule:: sknano.structures._unrolled_swnt

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from sknano.core.refdata import dVDW

from ._mixins import UnrolledSWNTMixin
from ._base import StructureBase

__all__ = ['UnrolledSWNT']


class UnrolledSWNT(UnrolledSWNTMixin, StructureBase):

    def __init__(self, nx=1, ny=1, Lx=None, Ly=None, nlayers=1,
                 layer_spacing=dVDW, stacking_order='AB', **kwargs):

        super(UnrolledSWNT, self).__init__(**kwargs)

        self.nx = nx
        self.ny = ny

    @property
    def Lx(self):
        return self.nx * self.Ch / 10

    @property
    def Ly(self):
        return self.ny * dVDW / 10
