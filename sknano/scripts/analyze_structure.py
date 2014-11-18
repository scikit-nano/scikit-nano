# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
================================================================
Command line script (:mod:`sknano.scripts.analyze_structure`)
================================================================

.. currentmodule:: sknano.scripts.analyze_structure

"""
from __future__ import absolute_import, division, print_function

from collections import OrderedDict
import argparse
import os
import sys

import numpy as np

from sknano.io import StructureReader
from sknano.structures import get_Ch_indices

__all__ = ['analyze_structure']


def argparser():
    parser = argparse.ArgumentParser()

    parser.add_argument('--cfg', '--config', dest='analysis_config',
                        metavar='CONFIG_FILE', default=None,
                        help='config file with analysis '
                        'settings. (default: %(default)s)')
    parser.add_argument(
        '--structure-format', default=None, choices=('data', 'dump', 'xyz'),
        help='input structure data file format. If `None`, then guess '
        'format from file name. (default: %(default)s)')
    subparsers = parser.add_subparsers(title='sub-commands')

    poav_parser = subparsers.add_parser('POAV')
    poav_parser.add_argument('--atom-ids', metavar='ATOM_ID', type=int,
                             nargs='+',
                             help='One or more atom IDs to analyze. '
                             '(default: %(default)s)')
    poav_parser.add_argument('--include-NN', action='store_true',
                             help='analyze nearest neighbor atoms for '
                             'each atom in list of `atom-ids`. '
                             '(default: %(default)s)')
    poav_parser.set_defaults(analyze='POAV')

    parser.add_argument('structure_file',
                        help='input structure data files.')

    return parser


def analyze_structure(analyze=None, analysis_config=None,
                      structure_format=None, structure_file=None,
                      **kwargs):
    pass


def analyze_POAVs(atoms, atom_ids=None, include_NN=False):
    pass


def main():
    args = argparser().parse_args()
    analyze_structure(**vars(args))

if __name__ == '__main__':
    sys.main(main())
