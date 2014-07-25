# -*- coding: utf-8 -*-
"""
=============================================================
Parameter LUTs (:mod:`sknano.nanogen._parameter_luts`)
=============================================================

.. currentmodule:: sknano.nanogen._parameter_luts

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

param_units = {}
param_units['dt'] = \
    param_units['rt'] = \
    param_units['Ch'] = \
    param_units['T'] = \
    param_units['bond'] = u' \u212B'
param_units['chiral_angle'] = u'\u00b0'

param_symbols = {}
param_symbols['t1'] = u't\u2081'
param_symbols['t2'] = u't\u2082'
param_symbols['chiral_angle'] = u'\u03b8c'

param_strfmt = {}
param_strfmt['Ch'] = \
    param_strfmt['T'] = \
    param_strfmt['dt'] = \
    param_strfmt['rt'] = \
    param_strfmt['chiral_angle'] = '{:.2f}'
param_strfmt['bond'] = '{:.3f}'

__all__ = ['param_units', 'param_symbols', 'param_strfmt']
