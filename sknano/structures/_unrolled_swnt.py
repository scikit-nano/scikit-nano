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

from ._mixins import UnrolledSWNTMixin, SWNT

__all__ = ['UnrolledSWNT']


class UnrolledSWNT(UnrolledSWNTMixin, SWNT):

    def __init__(self, nx=1, ny=1, Lx=None, Ly=None, nlayers=1,
                 layer_spacing=dVDW, stacking_order='AB', **kwargs):

        super(UnrolledSWNT, self).__init__(**kwargs)

        self._nx = int(nx)
        self._ny = int(ny)
        self._Lx = Lx
        self._Ly = Ly

        self.compute_tube_params()

    def compute_tube_params(self):

        super(UnrolledSWNT, self).compute_tube_params()

        self._Lx = \
            self.compute_Lx(n=self._n, m=self._m, nx=self._nx, bond=self._bond)

        self._Ly = \
            self.compute_Ly(n=self._n, m=self._m, ny=self._ny, bond=self._bond)
