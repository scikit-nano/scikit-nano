#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
====================================================
Command line script (:mod:`sknano.scripts.nanogen`)
====================================================

CLI to :mod:`sknano.generators` tools.

.. currentmodule:: sknano.scripts.nanogen

This module allows for easy structure generation from the command line.

.. seealso::

   :class:`~sknano.generators.FullereneGenerator`

   :class:`~sknano.generators.GrapheneGenerator`

   :class:`~sknano.generators.BilayerGrapheneGenerator`

   :class:`~sknano.generators.UnrolledSWNTGenerator`

   :class:`~sknano.generators.SWNTGenerator`

   :class:`~sknano.generators.MWNTGenerator`

   :class:`~sknano.generators.LayeredStructureGenerator`

.. code-block:: sh

   > nanogen --help

     usage: nanogen [-h] [--fname FNAME] [--structure-format {data,xyz}]
               [--debug] [--verbose] [--version]
               {fullerene,graphene,bilayer_graphene,unrolled_swnt,swnt,mwnt,layered_structure}
               ...

     optional arguments:
       -h, --help            show this help message and exit
       --fname FNAME         structure file name
       --structure-format {data,xyz}
                             structure file format (default: xyz)
       --debug               debug output
       --verbose             verbose output
       --version             show nanogen's version number and exit

     sub-commands:
       {fullerene,graphene,bilayer_graphene,unrolled_swnt,swnt,mwnt,layered_structure}

.. autofunction:: nanogen

Examples
--------

The following command generates a graphene sheet with a 10 Å long
armchair edge and 1 Å long zigzag edge and saves the data in the
LAMMPS data format.::

    > nanogen --structure-format data graphene 10 1

This command will generate a :math:`(20, 0)` SWNT, 5 unit cells
long and saves the data in `xyz` format.::

    > nanogen swnt --Ch 20 0 --nz 5

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import argparse
import importlib
import sys

from sknano.core.refdata import aCC, element_data
# from sknano.core.structures import get_chiral_indices_from_str
from ._parser import add_default_arguments

_r_CC_vdw = element_data['C']['VanDerWaalsRadius']

__all__ = ['nanogen', 'nanogen_parser']


def nanogen_parser():
    """:mod:`~sknano.scripts.nanogen` script \
        :class:`~python:argparse.ArgumentParser`."""
    parser = argparse.ArgumentParser()

    parser.add_argument('--fname', help='structure file name')
    parser.add_argument('--structure-format', default='xyz',
                        choices=('data', 'xyz'),
                        help='structure file format (default: %(default)s)')
    parser = add_default_arguments(parser)

    subparsers = parser.add_subparsers(title='sub-commands')

    fullerene_parser = subparsers.add_parser('fullerene')
    fullerene_parser.add_argument('N', type=int, help='The N in C_N')
    fullerene_parser.set_defaults(generator_class='FullereneGenerator')

    nanostructure_parent_parser = argparse.ArgumentParser(add_help=False)
    nanostructure_parent_parser.add_argument(
        '--basis', type=str, nargs=2, metavar=('ELEMENT1', 'ELEMENT2'),
        default=['C', 'C'], help='2 element symbols or atomic numbers of the '
        'two atom basis (default: %(default)s)')
    nanostructure_parent_parser.add_argument(
        '--bond', type=float, default=aCC, help='Bond length between '
        'nearest-neighbor atoms. Must be in units of Angstroms. '
        '(default: %(default)s)')

    graphene_parent_parser = argparse.ArgumentParser(add_help=False)
    graphene_parent_parser.add_argument('--layer-spacing', type=float,
                                        default=2*_r_CC_vdw,
                                        help='distance between '
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
        subparsers.add_parser('graphene', parents=[nanostructure_parent_parser,
                                                   graphene_parent_parser])
    graphene_parser.add_argument('--nlayers', type=int, default=1,
                                 help='Number of graphene layers. '
                                 '(default: %(default)s)')
    graphene_subparsers = \
        graphene_parser.add_subparsers(title='graphene sub-commands')

    primitive_cell_graphene_parser = \
        graphene_subparsers.add_parser('from_primitive_cell')
    primitive_cell_graphene_parser.add_argument(
        '--edge-length', type=float, default=10,
        help='graphene edge length in **Angstroms**')
    primitive_cell_graphene_parser.set_defaults(
        generator_class='PrimitiveCellGrapheneGenerator')

    conventional_cell_graphene_parser = \
        graphene_subparsers.add_parser('from_conventional_cell')
    conventional_cell_graphene_parser.add_argument(
        '--armchair-edge-length', type=float, default=10,
        help='length of graphene armchair edge in **Angstroms**')
    conventional_cell_graphene_parser.add_argument(
        '--zigzag-edge-length', type=float, default=10,
        help='length of graphene zigzag edge in **Angstroms**')
    conventional_cell_graphene_parser.set_defaults(
        generator_class='ConventionalCellGrapheneGenerator')

    # graphene_parser.set_defaults(generator_class='GrapheneGenerator')

    bilayer_graphene_parser = \
        subparsers.add_parser('bilayer_graphene',
                              parents=[nanostructure_parent_parser,
                                       graphene_parent_parser])
    bilayer_graphene_parser.set_defaults(
        generator_class='BilayerGrapheneGenerator')

    swnt_parent_parser = argparse.ArgumentParser(add_help=False)
    swnt_parent_parser.add_argument('--Ch', type=int, nargs=2,
                                    metavar=('n', 'm'),
                                    default=(5, 0),
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
                                          'axis in **Angstroms**. '
                                          '(default: %(default)s)')
    swnt_parent_parser.add_argument(
        '--fix-Lz', action='store_true', help='Generate the nanotube with '
        'length as close to the specified `Lz` as possible. If `True`, then '
        'non integer `nz` cells are permitted. (default: %(default)s)')

    unrolled_swnt_parser = \
        subparsers.add_parser('unrolled_swnt',
                              parents=[nanostructure_parent_parser,
                                       swnt_parent_parser])

    unrolled_swnt_parser_group = \
        unrolled_swnt_parser.add_mutually_exclusive_group()
    unrolled_swnt_parser_group.add_argument(
        '--nx', type=int, default=1,
        help='Number of repeat unit cells along `x` axis '
        '(default: %(default)s)')
    unrolled_swnt_parser_group.add_argument(
        '--Lx', type=float, default=None,
        help='Length of unrolled nanotube along `x` axis in **Angstroms**. '
        '(default: %(default)s)')
    unrolled_swnt_parser.add_argument(
        '--fix-Lx', action='store_true', help='Generate the nanotube with '
        'length as close to the specified `Lx` as possible. If `True`, then '
        'non integer `nx` cells are permitted. (default: %(default)s)')

    unrolled_swnt_parser.set_defaults(generator_class='UnrolledSWNTGenerator')

    mwnt_parent_parser = argparse.ArgumentParser(add_help=False)
    mwnt_parent_parser.add_argument('--Ch', default=None, metavar=('n', 'm'),
                                    type=int, nargs=2, action='append',
                                    help='Use this command multiple times to '
                                    'generate a list of (`n`, `m`) '
                                    'chiralities for each `SWNT` wall in '
                                    '`MWNT` (default: %(default)s)')
    mwnt_parent_parser.add_argument('--Nwalls', type=int,
                                    default=3,
                                    help='Number of `SWNT` walls in `MWNT` '
                                    '(default: %(default)s)')
    mwnt_parent_parser.add_argument('--Lz', type=float, default=None,
                                    help='Length of nanotube along `z` axis '
                                    'in **Angstroms**. '
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

    bundle_parent_parser = argparse.ArgumentParser(add_help=False)
    bundle_parent_parser.add_argument('--nx', type=int, default=1,
                                      help='Number of repeat unit cells '
                                      'along `x` axis (default: %(default)s)')
    bundle_parent_parser.add_argument('--ny', type=int, default=1,
                                      help='Number of repeat unit cells '
                                      'along `y` axis (default: %(default)s)')
    bundle_parent_parser.add_argument('--vdw-radius', type=float,
                                      default=_r_CC_vdw,
                                      help='van der Waals '
                                      'radius of nanotube atoms '
                                      '(default: %(default)s)')

    swnt_parser = \
        subparsers.add_parser('swnt', parents=[nanostructure_parent_parser,
                                               swnt_parent_parser,
                                               bundle_parent_parser])
    swnt_parser.set_defaults(generator_class='SWNTGenerator')

    mwnt_parser = \
        subparsers.add_parser('mwnt', parents=[nanostructure_parent_parser,
                                               mwnt_parent_parser,
                                               bundle_parent_parser])
    mwnt_parser.set_defaults(generator_class='MWNTGenerator')

    layered_structure_parser = subparsers.add_parser('layered_structure')
    layered_structure_parser.add_argument('cfgfile')
    layered_structure_parser.set_defaults(
        generator_class='LayeredStructureGenerator')
    return parser


def nanogen(generator_class=None, fname=None, structure_format='xyz',
            **kwargs):
    """Function for the nanogen script.

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
    args = nanogen_parser().parse_args()
    if hasattr(args, 'generator_class'):
        nanogen(**vars(args))
    else:
        nanogen_parser().print_help()


if __name__ == '__main__':
    sys.exit(main())
