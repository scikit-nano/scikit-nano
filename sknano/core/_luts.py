# -*- coding: utf-8 -*-
"""
==============================================================
Preset data structures (:mod:`sknano.core._luts`)
==============================================================

.. currentmodule:: sknano.core._luts

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

__all__ = ['chiral_type_name_mappings', 'components', 'dimensions', 'xyz',
           'xyz_axes', 'filter_key_type_mappings']

chiral_type_name_mappings = \
    {'achiral': 'aCh', 'armchair': 'AC', 'zigzag': 'ZZ', 'chiral': 'Ch'}

components = dimensions = xyz = xyz_axes = ('x', 'y', 'z')

filter_key_type_mappings = {}
filter_key_type_mappings['Ch_type'] = str
for k in ('even', 'odd'):
    filter_key_type_mappings[k + '_only'] = bool

for k in ('min_index', 'max_index',
          'min_n', 'max_n',
          'min_m', 'max_m',
          'n', 'm'):
    filter_key_type_mappings[k] = int
