# -*- coding: utf-8 -*-
"""
=========================================================
Reference data (:mod:`sknano.core.refdata._bonds`)
=========================================================

.. currentmodule:: sknano.core.refdata._bonds

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

# import json
import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

# from sknano.core.atoms import Bond

aCC = C_C = CCbond = 1.42  # angstroms
C_H = CHbond = 1.09  # angstroms
dVDW = 3.35  # angstroms

__all__ = ['dump_bond_data', 'load_bond_data',
           'aCC', 'C_C', 'CCbond', 'C_H', 'CHbond', 'dVDW']

_bond_file = 'bonds.yaml'


def dump_bond_data(data):
    with open(_bond_file, 'w') as f:
        # json.dump(data, f, indent=1, separators=(',', ': '))
        yaml.dump(data, f, Dumper=Dumper)


def load_bond_data():
    data = None
    with open(_bond_file) as f:
        data = yaml.load(f, Loader=Loader)
    return data
