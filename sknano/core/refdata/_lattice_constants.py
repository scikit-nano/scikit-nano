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

__all__ = ['lattice_constants', 'lattice_parameters',
           'BCC_lattice_parameters', 'FCC_lattice_parameters',
           'trigonal_lattice_parameters']


lattice_constants = lattice_parameters = {}
lattice_parameters.update({'graphene': {'a': 2.46}})
lattice_parameters.update({'alpha_quartz': {'a': 4.916, 'c': 5.405}})
lattice_parameters.update(dict.fromkeys(['MoS2', 'molybdenum_disulphide'],
                                        {'a': 3.160, 'c': 12.294}))
lattice_parameters.update({'diamond': 3.567})
lattice_parameters.update({'caesium_chloride': 4.123})
lattice_parameters.update({'rock_salt': 5.406})
lattice_parameters.update({'zincblende': 5.406})

BCC_lattice_parameters = {}
BCC_lattice_parameters.update({'Ba': 5.028})
BCC_lattice_parameters.update({'Fe': 2.8665})

FCC_lattice_parameters = {}
FCC_lattice_parameters.update({'Ac': 5.67})
FCC_lattice_parameters.update({'Al': 4.0495})
FCC_lattice_parameters.update({'Ar': 5.256})
FCC_lattice_parameters.update({'Au': 4.078})
FCC_lattice_parameters.update({'Cu': 3.615})

HCP_lattice_parameters = {}
HCP_lattice_parameters.update({'Am': {'a': 3.4681, 'c': 11.241}})
HCP_lattice_parameters.update({'Bk': {'a': 3.416, 'c': 11.069}})
HCP_lattice_parameters.update({'Be': {'a': 2.2858, 'c': 3.5843}})

lattice_parameters.update({'cubic': {}})
lattice_parameters['cubic'].update({'BCC': BCC_lattice_parameters})
lattice_parameters['cubic'].update({'FCC': FCC_lattice_parameters})

trigonal_lattice_parameters = {}
trigonal_lattice_parameters.update({'Sb': {'a': 4.307, 'c': 11.273}})
trigonal_lattice_parameters.update({'As': {'a': 3.7598, 'c': 10.5475}})

lattice_parameters.update({'trigonal': trigonal_lattice_parameters})
