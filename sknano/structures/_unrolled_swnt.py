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

from ._mixins import NanotubeMixin, UnrolledSWNTMixin
from ._base import StructureBase

__all__ = ['UnrolledSWNT']


class UnrolledSWNT(UnrolledSWNTMixin, NanotubeMixin, StructureBase):
    """Unrolled SWNT structure class."""

    def __init__(self, n=10, m=0, nx=1, nz=1, Nlayers=1, layer_spacing=dVDW,
                 stacking_order='AB', Lx=None, fix_Lx=False, Lz=None,
                 fix_Lz=False, **kwargs):

        super(UnrolledSWNT, self).__init__(**kwargs)

        self.n = n
        self.m = m

        self.fix_Lx = fix_Lx
        if Lx is not None:
            self.nx = 10 * float(Lx) / self.Ch
        elif nx is not None:
            self.nx = nx
        else:
            self.nx = 1

        self.fix_Lz = fix_Lz
        if Lz is not None:
            self.nz = 10 * float(Lz) / self.T
        elif nz is not None:
            self.nz = nz
        else:
            self.nz = 1

        self.Nlayers = Nlayers
        self.layer_spacing = layer_spacing
        self.stacking_order = stacking_order
