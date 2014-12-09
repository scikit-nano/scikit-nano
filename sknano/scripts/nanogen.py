# -*- coding: utf-8 -*-
#!/usr/bin/env python
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
     usage: nanogen [-h] [--element1 ELEMENT1] [--element2 ELEMENT2]
                    [--bond BOND] [--fname FNAME]
                    [--structure-format {data,xyz}] [--verbose]
                    {graphene,BLG,swnt,unrolled_swnt,mwnt,swnt_bundle,
                    mwnt_bundle}
                    ...

     optional arguments:
       -h, --help            show this help message and exit
       --element1 ELEMENT1   element symbol or atomic number of basis atom 1
                             (default: C)
       --element2 ELEMENT2   element symbol or atomic number of basis atom 2
                             (default: C)
       --bond BOND           Bond length between nearest neighbor atoms. Must
                             in units of Angstroms. (default: 1.42)
       --fname FNAME         structure file name
       --structure-format {data,xyz}
                             structure file format (default: xyz)
       --verbose             verbose output
       --debug               debug output

     sub-commands:
       {graphene,BLG,swnt,unrolled_swnt,mwnt,swnt_bundle,mwnt_bundle}

.. autofunction:: nanogen

Examples
--------

The following command generates a graphene sheet 10 nm long by 1 nm wide
with an armchair edge pattern along its length and saves the data in the
LAMMPS data format.::

    > nanogen --structure-format data graphene 10 1 --edge AC

This command will generate a :math:`(20, 0)` SWNT, 5 unit cells
long and saves the data in `xyz` format.::

    > nanogen swnt -n 20 -m 0 --nz 5

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
    >>> for n, m in Ch_list:
    ...     SWNTGenerator(n=n, m=m).save_data(structure_format='data')

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import argparse
import importlib
import sys

from sknano.core.refdata import CCbond, dVDW
from sknano.version import short_version as version

__all__ = ['nanogen']


def argparser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--element1', default='C', help='element symbol or atomic number of '
        'basis atom 1 (default: %(default)s)')
    parser.add_argument(
        '--element2', default='C', help='element symbol or atomic number of '
        'basis atom 2 (default: %(default)s)')
    parser.add_argument(
        '--bond', type=float, default=CCbond, help='Bond length between '
        'nearest neighbor atoms. Must in units of Angstroms. '
        '(default: %(default)s)')
    parser.add_argument('--fname', help='structure file name')
    parser.add_argument(
        '--structure-format', default='xyz', choices=('data', 'xyz'),
        help='structure file format (default: %(default)s)')
    parser.add_argument(
        '--verbose', action='store_true', help='verbose output')
    parser.add_argument(
        '--debug', action='store_true', help='debug output')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(version),
                        help="show %(prog)s's version number and exit")

    subparsers = parser.add_subparsers(title='sub-commands')

    graphene_parent_parser = argparse.ArgumentParser(add_help=False)
    graphene_parent_parser.add_argument(
        'length', type=float, help='length of graphene sheet in '
        '**nanometers**')
    graphene_parent_parser.add_argument(
        'width', type=float, help='width of graphene sheet in **nanometers** ')
    graphene_parent_parser.add_argument(
        '--edge', choices=('AC', 'ZZ'), help='chirality along the '
        '**length** of sheet. If `None`, then one will chosen at random. '
        '(default: %(default)s)')
    graphene_parent_parser.add_argument(
        '--layer-spacing', type=float, default=3.35, help='distance between '
        'layers in **Angstroms**. (default: %(default)s)')
    graphene_parent_parser.add_argument(
        '--stacking-order', choices=('AA', 'AB'), default='AB',
        help='Stacking order of graphene layers (default: %(default)s)')

    graphene_parser = \
        subparsers.add_parser('graphene', parents=[graphene_parent_parser])
    graphene_parser.add_argument('--nlayers', type=int, default=1,
                                 help='Number of graphene layers. '
                                 '(default: %(default)s)')
    graphene_parser.set_defaults(generator_class='GrapheneGenerator')

    bilayer_graphene_parser = \
        subparsers.add_parser('bilayer_graphene',
                              parents=[graphene_parent_parser])
    bilayer_graphene_parser.add_argument(
        '--layer-rotation-angle', type=float, help='Rotation angle of second '
        'layer in **degrees.** (default: %(default)s)')
    bilayer_graphene_parser.set_defaults(
        generator_class='BilayerGrapheneGenerator')

    nanotube_parent_parser = argparse.ArgumentParser(add_help=False)
    nanotube_parent_parser_group = \
        nanotube_parent_parser.add_mutually_exclusive_group()
    nanotube_parent_parser_group.add_argument(
        '--nz', type=int, default=1,
        help='Number of repeat unit cells along `z` axis '
        '(default: %(default)s)')
    nanotube_parent_parser_group.add_argument(
        '--Lz', type=float, default=None,
        help='Length of nanotube along `z` axis in **nanometers**. '
        '(default: %(default)s)')
    nanotube_parent_parser.add_argument(
        '--fix-Lz', action='store_true', help='Generate the nanotube with '
        'length as close to the specified `Lz` as possible. If `True`, then '
        'non integer `nz` cells are permitted. (default: ' '%(default)s)')

    swnt_parent_parser = argparse.ArgumentParser(add_help=False)
    swnt_parent_parser.add_argument('-n', type=int, default=10,
                                    help='Chiral index `n` '
                                    '(default: %(default)s)')
    swnt_parent_parser.add_argument('-m', type=int, default=10,
                                    help='Chiral index `m` '
                                    '(default: %(default)s)')

    mwnt_parent_parser = argparse.ArgumentParser(add_help=False)

    swnt_parser = \
        subparsers.add_parser('swnt', parents=[swnt_parent_parser,
                                               nanotube_parent_parser])
    swnt_parser.set_defaults(generator_class='SWNTGenerator')

    unrolled_swnt_parser = \
        subparsers.add_parser('unrolled_swnt',
                              parents=[swnt_parent_parser,
                                       nanotube_parent_parser])

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
        'non integer `nx` cells are permitted. (default: ' '%(default)s)')

    unrolled_swnt_parser.set_defaults(generator_class='UnrolledSWNTGenerator')

    mwnt_parser = \
        subparsers.add_parser('mwnt', parents=[nanotube_parent_parser])
    mwnt_parser.set_defaults(generator_class='MWNTGenerator')

    bundle_parent_parser = argparse.ArgumentParser(add_help=False)
    bundle_parent_parser.add_argument('--nx', type=int, default=1,
                                      help='Number of repeat unit cells '
                                      'along `x` axis (default: %(default)s)')
    bundle_parent_parser.add_argument('--ny', type=int, default=1,
                                      help='Number of repeat unit cells '
                                      'along `y` axis (default: %(default)s)')

    swnt_bundle_parser = \
        subparsers.add_parser('swnt_bundle', parents=[swnt_parent_parser,
                                                      bundle_parent_parser])
    swnt_bundle_parser.add_argument(
        '--vdw-spacing', type=float, default=dVDW, help='van der Waals '
        'distance between nearest neighbor tubes (default: %(default)s)')
    swnt_bundle_parser.set_defaults(generator_class='SWNTBundleGenerator')

    mwnt_bundle_parser = \
        subparsers.add_parser('mwnt_bundle', parents=[mwnt_parent_parser,
                                                      bundle_parent_parser])
    mwnt_bundle_parser.set_defaults(generator_class='MWNTBundleGenerator')

    return parser


def nanogen(generator_class=None, element1='C', element2='C', bond=CCbond,
            fname=None, structure_format='xyz', **kwargs):
    """Generate nano-structure data.

    Parameters
    ----------
    generator_class : str
        nano-structure generator class
    element1 : str, optional
        element symbol or atomic number of basis atom 1
    element2 : str, optional
        element symbol or atomic number of basis atom 2
    bond : float, optional
        element1-element2 bond length in ``units``
    fname : str, optional
        structure file name
    structure_format : str, optional
        output file format

    """
    generator = getattr(importlib.import_module('sknano.generators'),
                        generator_class)
    generator(**kwargs).save_data(fname=fname,
                                  structure_format=structure_format)


def main():
    args = argparser().parse_args()
    if hasattr(args, 'generator_class'):
        nanogen(**vars(args))
    else:
        argparser().print_help()


if __name__ == '__main__':
    sys.exit(main())
