# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
===================================================
Command line script (:mod:`sknano.scripts.nanogen`)
===================================================

CLI to :mod:`sknano.nanogen` tools.

.. currentmodule:: sknano.scripts.nanogen

This module allows for easy structure generation from the command line.

The :mod:`sknano.nanogen` package provides a python class
(:py:class:`~sknano.nanogen.TubeGen`) which wraps the
functionality of the TubeGen ([TG]_) binary executable.
TubeGen generates structure data for nanotubes and graphene.

You must have a compiled binary executable of TubeGen available in the
system :envvar:`PATH` and callable as :program:`tubegen`.

.. seealso::

   :py:class:`~sknano.nanogen.TubeGen`

.. argparse::
   :module: sknano.scripts.nanogen
   :func: _argparser
   :prog: nanogen

Examples
--------

::

    > nanogen --chirality 10 10 --cell-count 1 1 5

"""
from __future__ import absolute_import, division, print_function

import argparse
import sys

from pksci.chemistry import Atom
from pkshared.tools.refdata import ccbond

from ..nanogen import TubeGen, format_ext
from ..structure_io import XYZ2DATAConverter


def _argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--format', default='xyz', dest='fmt',
                        choices=tuple([fmt for fmt in format_ext.iterkeys()]),
                        help='output file format (default: %(default)s)')
    parser.add_argument('--units', default='angstrom',
                        choices=('angstrom', 'bohr'),
                        help='units (default: %(default)s)')
    parser.add_argument('--element', nargs=2, metavar=('1', '2'),
                        default=('C', 'C'),
                        help='atomic symbol or number of basis '
                        'atom elements (default: %(default)s)')
    parser.add_argument('--bond-length', type=float, default=ccbond,
                        help='element1-element2 bond length in Angstroms '
                        '(default: %(default)s)')
    parser.add_argument('--gutter', nargs=3, type=float,
                        metavar=('X', 'Y', 'Z'), default=(1.6735, 1.6735, 0),
                        help='cell gutter (default: %(default)s)')
    #parser.add_argument('--structure', default='nanotube',
    #        choices=('nanotube', 'graphene', 'hexagonal-bundle'),
    #        help='crystal structure (default: %(default)s)')
    parser.add_argument('--shape', default='hexagonal',
                        choices=('hexagonal', 'cubic', 'planar'),
                        help='crystal structure (default: %(default)s)')
    parser.add_argument('--chirality', nargs=2, type=int, metavar=('n', 'm'),
                        default=(10, 10), help='chiral indices (n,m) '
                        '(default: %(default)s)')
    parser.add_argument('--relax-tube', default='yes', choices=('yes', 'no'),
                        help='select tube relaxation (default: %(default)s)')

    parser_group = parser.add_mutually_exclusive_group()
    parser_group.add_argument('--cell-count', nargs=3, type=int,
                              metavar=('X', 'Y', 'Z'), default=(1, 1, 1),
                              help='number of unit cells '
                              '(default: %(default)s)')
    parser_group.add_argument('--tube-length', type=float,
                              metavar='L_tube', default=None,
                              help='tube length in nanometers '
                              '(default: %(default)s)')
    #parser.add_argument('--ion-element', default=None, help='ion element '
    #                    'chemical symbol (default: %(default)s)')
    parser.add_argument('--xyz2data', action='store_true',
                        help='convert xyz to LAMMPS data format')
    parser.add_argument('--add-atomtype', nargs=2, action='append',
                        metavar=('ELEMENT', 'ATOM-TYPE'),
                        dest='new_atomtypes',
                        help='add new atom to structure data.')

    return parser


def call_tubegen(args):
    """Call TubeGen."""

    tubegen = TubeGen(fmt=args.fmt, units=args.units, bond=args.bond_length,
                      element1=args.element[0], element2=args.element[1],
                      gutter=args.gutter, shape=args.shape,
                      chirality=args.chirality, cell_count=args.cell_count,
                      tube_length=args.tube_length, relax_tube=args.relax_tube)

    tubegen.generate()
    tubegen.cleanup()

    if args.tube_length is not None:
        tubegen.fix_length()

    if args.xyz2data:
        if args.fmt == 'xyz':
            xyzconverter = XYZ2DATAConverter(tubegen.output, pad_box=True)
            if args.new_atomtypes is not None:
                for element, atomtype in args.new_atomtypes:
                    atom = Atom(element, atomtype=int(atomtype))
                    xyzconverter.add_atomtype(atom=atom)
            xyzconverter.convert()
        else:
            print("can't convert {} format to data. ".format(args.fmt) +
                  "Must be xyz format.")


def main():

    args = _argparser().parse_args()
    call_tubegen(args)

if __name__ == '__main__':
    sys.exit(main())
