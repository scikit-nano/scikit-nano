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

   :class:`~sknano.generators.GrapheneGenerator`

   :class:`~sknano.generators.BilayerGrapheneGenerator`

   :class:`~sknano.generators.UnrolledSWNTGenerator`

   :class:`~sknano.generators.SWNTGenerator`

   :class:`~sknano.generators.SWNTBundleGenerator`

   :class:`~sknano.generators.MWNTGenerator`

   :class:`~sknano.generators.MWNTBundleGenerator`

.. code-block:: python

   > nanogen --help

     usage: nanogen [-h] [--basis ELEMENT1 ELEMENT2] [--bond BOND] [--fname FNAME]
                    [--structure-format {data,xyz}] [--verbose] [--debug]
                    [--version]
                    {graphene,bilayer_graphene,swnt,unrolled_swnt,mwnt,swnt_bundle,mwnt_bundle}
                    ...

     optional arguments:
       -h, --help            show this help message and exit
       --basis ELEMENT1 ELEMENT2
                             2 element symbols or atomic numbers of the two atom
                             basis (default: ['C', 'C'])
       --bond BOND           Bond length between nearest-neighbor atoms. Must be in
                             units of Angstroms. (default: 1.42)
       --fname FNAME         structure file name
       --structure-format {data,xyz}
                             structure file format (default: xyz)
       --verbose             verbose output
       --debug               debug output
       --version             show nanogen's version number and exit

     sub-commands:
       {graphene,bilayer_graphene,swnt,unrolled_swnt,mwnt,swnt_bundle,mwnt_bundle}

.. autofunction:: nanogen

Examples
--------

The following command generates a graphene sheet with a 10 nm long
armchair edge and 1 nm long zigzag edge and saves the data in the
LAMMPS data format.::

    > nanogen --structure-format data graphene 10 1

This command will generate a :math:`(20, 0)` SWNT, 5 unit cells
long and saves the data in `xyz` format.::

    > nanogen swnt --Ch 20 0 --nz 5

Notes
-----
The :mod:`nanogen` script calls the :func:`~sknano.scripts.nanogen.nanogen`
function and neither the script nor the function allow for generating
multiple structures per call. It's easy to do this using the scripting
capabilities of the shell or from within an interactive Python session.

For example, if you have a list of chiralities that you want structure data
for, you can do something like this from within an interactive session::

    >>> from sknano.structures import generate_Ch_list
    >>> from sknano.generators import SWNTGenerator
    >>> # Generate your list of (n, m) chiralities
    >>> Ch_list = generate_Ch_list(ni=5, nf=25, mi=0, mf=25, handedness='right')
    >>> for Ch in Ch_list:
    ...     SWNTGenerator(Ch).save(structure_format='data')

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import argparse
import importlib
import sys

from sknano.core.refdata import aCC, element_data
from ._parser import add_default_arguments

_r_CC_vdw = element_data['C']['VanDerWaalsRadius']

__all__ = ['nanogen', 'nanogen_parser']


def nanogen_parser():
    """:mod:`~sknano.scripts.nanogen` script \
        :class:`~python:argparse.ArgumentParser`."""
    parser = argparse.ArgumentParser()

    parser.add_argument('--basis', type=str, nargs=2,
                        metavar=('ELEMENT1', 'ELEMENT2'),
                        default=['C', 'C'],
                        help='2 element symbols or atomic numbers of the two '
                        'atom basis (default: %(default)s)')
    parser.add_argument('--bond', type=float, default=aCC,
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
    bundle_parent_parser.add_argument('--vdw-radius', type=float,
                                      default=_r_CC_vdw,
                                      help='van der Waals '
                                      'radius of nanotube atoms '
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
    args = nanogen_parser().parse_args()
    if hasattr(args, 'generator_class'):
        nanogen(**vars(args))
    else:
        nanogen_parser().print_help()


if __name__ == '__main__':
    sys.exit(main())
