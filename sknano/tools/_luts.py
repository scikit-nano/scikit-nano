# -*- coding: utf-8 -*-
"""
==============================================================
Preset data structures (:mod:`sknano.tools._luts`)
==============================================================

.. currentmodule:: sknano.tools._luts

"""
from __future__ import division, absolute_import, print_function
__docformat__ = 'restructuredtext'

__all__ = ['chiral_type_name_mappings', 'xyz_axes']

chiral_type_name_mappings = \
    {'achiral': 'aCh', 'armchair': 'AC', 'zigzag': 'ZZ', 'chiral': 'Ch'}

xyz_axes = ('x', 'y', 'z')
