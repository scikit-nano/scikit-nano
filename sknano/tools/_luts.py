# -*- coding: utf-8 -*-
"""
==============================================================
Preset data structures (:mod:`sknano.tools._luts`)
==============================================================

.. currentmodule:: sknano.tools._luts

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

__all__ = ['chiral_type_name_mappings',
           'components', 'dimensions',
           'xyz', 'xyz_axes']

chiral_type_name_mappings = \
    {'achiral': 'aCh', 'armchair': 'AC', 'zigzag': 'ZZ', 'chiral': 'Ch'}

components = dimensions = xyz = xyz_axes = ('x', 'y', 'z')
