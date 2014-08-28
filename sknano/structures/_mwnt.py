# -*- coding: utf-8 -*-
"""
==============================================================================
MWNT structure class (:mod:`sknano.structures._mwnt`)
==============================================================================

.. currentmodule:: sknano.structures._mwnt

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core.refdata import dVDW

from ._mixins import MWNTMixin
from ._swnt import SWNT

__all__ = ['MWNT']


class MWNT(MWNTMixin, SWNT):

    def __init__(self, add_inner_shells=True, add_outer_shells=False,
                 max_shells=None, max_shell_diameter=np.inf,
                 min_shells=None, min_shell_diameter=0.0,
                 new_shell_type=None, shell_spacing=dVDW, **kwargs):

        super(MWNT, self).__init__(**kwargs)

        self._add_inner_shells = add_inner_shells
        self._add_outer_shells = add_outer_shells
        self._starting_shell_position = 'outer'

        self._max_shells = max_shells
        if max_shells is None:
            self._max_shells = 10
        self._max_shell_diameter = max_shell_diameter

        self._min_shells = min_shells
        if min_shells is None:
            self._min_shells = 2
        self._min_shell_diameter = min_shell_diameter

        self._new_shell_type = new_shell_type
        self._shell_spacing = shell_spacing

        self._Nshells_per_tube = 1
        self._Natoms_per_tube = 0

        self.compute_tube_params()

    def compute_tube_params(self):
        super(MWNT, self).compute_tube_params()
