#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
====================================================
Command line script (:mod:`sknano.scripts.sknano`)
====================================================

CLI to :mod:`sknano` tools.

.. currentmodule:: sknano.scripts.sknano

.. code-block:: python

   > sknano --help

.. autofunction:: sknano

Examples
--------


"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import argparse
import importlib
import sys

from sknano.core.refdata import CCbond, dVDW
from ._parser import add_default_arguments

__all__ = ['nanogen']


def argparser():
    parser = argparse.ArgumentParser()

    parser.add_argument('--basis', type=str, nargs=2,
                        metavar=('ELEMENT1', 'ELEMENT2'),
                        default=['C', 'C'],
                        help='2 element symbols or atomic numbers of the two '
                        'atom basis (default: %(default)s)')
    parser.add_argument('--bond', type=float, default=CCbond,
                        help='Bond length between nearest-neighbor atoms. '
                        'Must be in units of Angstroms. '
                        '(default: %(default)s)')
    parser.add_argument('--fname', help='structure file name')
    parser.add_argument('--structure-format', default='xyz',
                        choices=('data', 'xyz'),
                        help='structure file format (default: %(default)s)')
    parser = add_default_arguments(parser)

    subparsers = parser.add_subparsers(title='sub-commands')

    graphene_parent_parser = argparse.ArgumentParser(add_help=False)
    graphene_parent_parser.add_argument('--layer-spacing', type=float,
                                        default=dVDW, help='distance between '
                                        'graphene layers in **Angstroms**. '
                                        '(default: %(default)s)')
    graphene_parent_parser.add_argument('--layer-rotation-increment',
                                        type=float,
                                        default=None,
                                        help='Layer rotation angles of '
                                        'each layer in **degrees** '
                                        '(default: %(default)s)')
    graphene_parent_parser.add_argument('--layer-rotation-angles', type=list,
                                        default=None,
                                        help='Layer rotation angles of '
                                        'each layer in **degrees** '
                                        '(default: %(default)s)')
    graphene_parent_parser.add_argument('--stacking-order', default='AB',
                                        choices=('AA', 'AB'),
                                        help='Stacking order of graphene '
                                        'layers (default: %(default)s)')

    graphene_parser = \
        subparsers.add_parser('graphene', parents=[graphene_parent_parser])
    graphene_parser.add_argument('--nlayers', type=int, default=1,
                                 help='Number of graphene layers. '
                                 '(default: %(default)s)')
    graphene_subparsers = \
        graphene_parser.add_subparsers(title='graphene sub-commands')

    primitive_cell_graphene_parser = \
        graphene_subparsers.add_parser('from_primitive_cell')
    primitive_cell_graphene_parser.add_argument(
        '--edge-length', type=float, default=5,
        help='graphene edge length in **nanometers**')
    primitive_cell_graphene_parser.set_defaults(
        generator_class='PrimitiveCellGrapheneGenerator')

    conventional_cell_graphene_parser = \
        graphene_subparsers.add_parser('from_conventional_cell')
    conventional_cell_graphene_parser.add_argument(
        '--armchair-edge-length', type=float, default=5,
        help='length of graphene armchair edge in **nanometers**')
    conventional_cell_graphene_parser.add_argument(
        '--zigzag-edge-length', type=float, default=5,
        help='length of graphene zigzag edge in **nanometers**')
    conventional_cell_graphene_parser.set_defaults(
        generator_class='ConventionalCellGrapheneGenerator')

    # graphene_parser.set_defaults(generator_class='GrapheneGenerator')

    bilayer_graphene_parser = \
        subparsers.add_parser('bilayer_graphene',
                              parents=[graphene_parent_parser])
    bilayer_graphene_parser.set_defaults(
        generator_class='BilayerGrapheneGenerator')

    swnt_parent_parser = argparse.ArgumentParser(add_help=False)
    swnt_parent_parser.add_argument('--Ch', type=int, nargs=2,
                                    metavar=('n', 'm'),
                                    default=(10, 10),
                                    help='Chiral indices (`n`, `m`) '
                                    '(default: %(default)s)')
    swnt_parent_parser_group = \
        swnt_parent_parser.add_mutually_exclusive_group()
    swnt_parent_parser_group.add_argument('--nz', type=int, default=1,
                                          help='Number of repeat unit cells '
                                          'along `z` axis '
                                          '(default: %(default)s)')
    swnt_parent_parser_group.add_argument('--Lz', type=float, default=None,
                                          help='Length of nanotube along `z` '
                                          'axis in **nanometers**. '
                                          '(default: %(default)s)')
    swnt_parent_parser.add_argument(
        '--fix-Lz', action='store_true', help='Generate the nanotube with '
        'length as close to the specified `Lz` as possible. If `True`, then '
        'non integer `nz` cells are permitted. (default: %(default)s)')

    swnt_parser = \
        subparsers.add_parser('swnt', parents=[swnt_parent_parser])
    swnt_parser.set_defaults(generator_class='SWNTGenerator')

    unrolled_swnt_parser = \
        subparsers.add_parser('unrolled_swnt', parents=[swnt_parent_parser])

    unrolled_swnt_parser_group = \
        unrolled_swnt_parser.add_mutually_exclusive_group()
    unrolled_swnt_parser_group.add_argument(
        '--nx', type=int, default=1,
        help='Number of repeat unit cells along `x` axis '
        '(default: %(default)s)')
    unrolled_swnt_parser_group.add_argument(
        '--Lx', type=float, default=None,
        help='Length of unrolled nanotube along `x` axis in **nanometers**. '
        '(default: %(default)s)')
    unrolled_swnt_parser.add_argument(
        '--fix-Lx', action='store_true', help='Generate the nanotube with '
        'length as close to the specified `Lx` as possible. If `True`, then '
        'non integer `nx` cells are permitted. (default: %(default)s)')

    unrolled_swnt_parser.set_defaults(generator_class='UnrolledSWNTGenerator')

    mwnt_parent_parser = argparse.ArgumentParser(add_help=False)
    mwnt_parent_parser.add_argument('--Ch-list', type=list,
                                    default=None,
                                    help='list of (`n`, `m`) chiralities '
                                    'for each `SWNT` wall in `MWNT` '
                                    '(default: %(default)s)')
    mwnt_parent_parser.add_argument('--Nwalls', type=int,
                                    default=3,
                                    help='Number of `SWNT` walls in `MWNT` '
                                    '(default: %(default)s)')
    mwnt_parent_parser.add_argument('--Lz', type=float, default=None,
                                    help='Length of nanotube along `z` axis '
                                    'in **nanometers**. '
                                    '(default: %(default)s)')
    mwnt_parent_parser.add_argument('--min-wall-diameter', type=float,
                                    default=None,
                                    help='Minimum `MWNT` wall diameter, in '
                                    'units of **Angstroms** '
                                    '(default: %(default)s)')
    mwnt_parent_parser.add_argument('--max-wall-diameter', type=float,
                                    default=None,
                                    help='Maximum `MWNT` wall diameter, in '
                                    'units of **Angstroms** '
                                    '(default: %(default)s)')
    mwnt_parent_parser.add_argument('--max-walls', type=int, default=None,
                                    help='Maximum number of `MWNT` walls '
                                    '(default: %(default)s)')

    mwnt_parser = \
        subparsers.add_parser('mwnt', parents=[mwnt_parent_parser])
    mwnt_parser.set_defaults(generator_class='MWNTGenerator')

    bundle_parent_parser = argparse.ArgumentParser(add_help=False)
    bundle_parent_parser.add_argument('--nx', type=int, default=1,
                                      help='Number of repeat unit cells '
                                      'along `x` axis (default: %(default)s)')
    bundle_parent_parser.add_argument('--ny', type=int, default=1,
                                      help='Number of repeat unit cells '
                                      'along `y` axis (default: %(default)s)')
    bundle_parent_parser.add_argument('--vdw-spacing', type=float,
                                      default=dVDW, help='van der Waals '
                                      'distance between nanotubes '
                                      '(default: %(default)s)')

    swnt_bundle_parser = \
        subparsers.add_parser('swnt_bundle', parents=[swnt_parent_parser,
                                                      bundle_parent_parser])
    swnt_bundle_parser.set_defaults(generator_class='SWNTBundleGenerator')

    mwnt_bundle_parser = \
        subparsers.add_parser('mwnt_bundle', parents=[mwnt_parent_parser,
                                                      bundle_parent_parser])
    mwnt_bundle_parser.set_defaults(generator_class='MWNTBundleGenerator')

    return parser


def nanogen(generator_class=None, fname=None, structure_format='xyz',
            **kwargs):
    """Generate nano-structure data.

    Parameters
    ----------
    generator_class : str
        nano-structure generator class
    fname : str, optional
        structure file name
    structure_format : str, optional
        output file format

    """
    generator = getattr(importlib.import_module('sknano.generators'),
                        generator_class)
    generator(**kwargs).save(fname=fname,
                             structure_format=structure_format)


def main():
    args = argparser().parse_args()
    if hasattr(args, 'generator_class'):
        nanogen(**vars(args))
    else:
        argparser().print_help()


if __name__ == '__main__':
    sys.exit(main())
