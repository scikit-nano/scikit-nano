# -*- coding: utf-8 -*-
"""
==============================================================================
Lattice constants (:mod:`sknano.core.refdata.lattice_constants`)
==============================================================================

.. currentmodule:: sknano.core.refdata.lattice_constants

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

# import json

import numpy as np

# from sknano.core import flatten
from . import aCC, element_data

__all__ = ['lattice_constants', 'lattice_parameters',
           'caesium_chloride_lattice_parameters', 'diamond_lattice_parameters',
           'rocksalt_lattice_parameters', 'zincblende_lattice_parameters',
           'wurtzite_lattice_parameters',
           'sc_lattice_parameters', 'bcc_lattice_parameters',
           'fcc_lattice_parameters', 'hexagonal_lattice_parameters',
           'trigonal_lattice_parameters']

r_CC_vdw = element_data['C']['VanDerWaalsRadius']


lattice_constants = lattice_parameters = {}

sc_lattice_parameters = {}
caesium_chloride_lattice_parameters = {}
caesium_chloride_lattice_parameters.update(
    dict.fromkeys(['caesium_chloride', 'CsCl'], 4.123))
sc_lattice_parameters.update(caesium_chloride_lattice_parameters)

bcc_lattice_parameters = {}
bcc_lattice_parameters.update({'Ba': 5.028})
bcc_lattice_parameters.update({'Fe': 2.8665})

fcc_lattice_parameters = {}
fcc_lattice_parameters.update({'Ac': 5.67})
fcc_lattice_parameters.update({'Al': 4.0495})
fcc_lattice_parameters.update({'Ar': 5.256})
fcc_lattice_parameters.update({'Au': 4.0782})
fcc_lattice_parameters.update({'Ag': 4.0853})
fcc_lattice_parameters.update({'Cu': 3.615})
diamond_lattice_parameters = {}
diamond_lattice_parameters.update(dict.fromkeys(['C', 'diamond'], 3.567))
diamond_lattice_parameters.update({'Si': 5.431})
diamond_lattice_parameters.update({'Ge': 5.658})
fcc_lattice_parameters.update(diamond_lattice_parameters)

rocksalt_lattice_parameters = {}
rocksalt_lattice_parameters.update(
    dict.fromkeys(['rocksalt', 'rock_salt', 'NaCl'], 5.640))
fcc_lattice_parameters.update(rocksalt_lattice_parameters)

zincblende_lattice_parameters = {}
zincblende_lattice_parameters.update(
    dict.fromkeys(['zincblende', 'zinc_blende', 'ZnS'], 5.420))
fcc_lattice_parameters.update(zincblende_lattice_parameters)

lattice_parameters.update({'cubic': {}})
lattice_parameters['cubic'].update(
    dict.fromkeys(['SC', 'sc', 'PC', 'pc', 'simple_cubic', 'primitive_cubic',
                   'cP', 'cp'], sc_lattice_parameters))
lattice_parameters['cubic'].update(
    dict.fromkeys(['BCC', 'bcc', 'body-centered_cubic', 'cI', 'ci'],
                  bcc_lattice_parameters))
lattice_parameters['cubic'].update(
    dict.fromkeys(['FCC', 'fcc', 'face-centered_cubic', 'cF', 'cf'],
                  fcc_lattice_parameters))

hexagonal_lattice_parameters = {}
hexagonal_lattice_parameters.update(
    dict.fromkeys(['MoS2', 'molybdenum_disulphide'], (3.160, 12.294)))
hexagonal_lattice_parameters.update(
    dict.fromkeys(['graphene', 'C'], (np.sqrt(3) * aCC, 2 * r_CC_vdw)))
hexagonal_lattice_parameters.update({'graphite': (2.461, 6.708)})
hexagonal_lattice_parameters.update({'alpha_quartz': (4.9134, 5.4052)})
hexagonal_lattice_parameters.update({'beta_quartz': (4.9965, 5.4546)})
hexagonal_lattice_parameters.update({'Am': (3.4681, 11.241)})
hexagonal_lattice_parameters.update({'Bk': (3.416, 11.069)})
hexagonal_lattice_parameters.update({'Be': (2.2858, 3.5843)})

wurtzite_lattice_parameters = {}
wurtzite_lattice_parameters.update(
    dict.fromkeys(['wurtzite', 'ZnS'], (3.82, 6.26)))
hexagonal_lattice_parameters.update(wurtzite_lattice_parameters)
lattice_parameters.update({'hexagonal': hexagonal_lattice_parameters})

trigonal_lattice_parameters = {}
trigonal_lattice_parameters.update({'Sb': (4.307, 11.273)})
trigonal_lattice_parameters.update({'As': (3.7598, 10.5475)})

lattice_parameters.update({'trigonal': trigonal_lattice_parameters})
