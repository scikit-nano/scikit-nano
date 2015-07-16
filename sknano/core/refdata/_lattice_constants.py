# -*- coding: utf-8 -*-
"""
==============================================================================
Lattice constants (:mod:`sknano.core.refdata._lattice_constants`)
==============================================================================

.. currentmodule:: sknano.core.refdata._lattice_constants

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

# import json

__all__ = ['lattice_constants', 'lattice_parameters']


lattice_constants = lattice_parameters = {}
lattice_parameters.update({'graphene': {'a': 2.46}})
lattice_parameters.update({'alpha_quartz': {'a': 4.916, 'c': 5.405}})
lattice_parameters.update(dict.fromkeys(['MoS2', 'molybdenum_disulphide'],
                                        {'a': 3.160, 'c': 12.294}))
lattice_parameters.update({'diamond': 3.567})
lattice_parameters.update({'caesium_chloride': 4.123})
lattice_parameters.update({'rock_salt': 5.406})
lattice_parameters.update({'zincblende': 5.406})

# lattice_parameters.update({'gold': 4.078})
# lattice_parameters.update({'copper': 3.615})
lattice_parameters.update({'cubic': {'FCC': {'Au': 4.078, 'Cu': 3.615}}})
