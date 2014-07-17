# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
====================================================
Command line script (:mod:`sknano.scripts.nanogen`)
====================================================

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

   :py:class:`~sknano.nanogen.GrapheneGenerator`

   :py:class:`~sknano.nanogen.NanotubeBundleGenerator`

   :py:class:`~sknano.nanogen.MWNTGenerator`

.. code-block:: python

   > nanogen --help
   usage: nanogen [-h]
                  [--format {pdb,pdb-pbc,cif,bgf,pov,wien,xyz,gaussian,
                  gaussian-pbc}]
                  [--units {angstrom,bohr}] [--element1 ELEMENT1]
                  [--element2 ELEMENT2] [--bond BOND] [--gutter X Y Z]
                  [--shape {hexagonal,cubic,planar}] [--chirality n m]
                  [--relax-tube {yes,no}]
                  [--cell-count X Y Z | --tube-length L_tube] [--xyz2data]
                  [--add-atomtype ELEMENT ATOM-TYPE]
                  [--generator {tubegen,nanogen}]

   optional arguments:
     -h, --help            show this help message and exit
     --format {pdb,pdb-pbc,cif,bgf,pov,wien,xyz,gaussian,gaussian-pbc}
     --units {angstrom,bohr}
     --element1 ELEMENT1
     --element2 ELEMENT2
     --bond BOND
     --gutter X Y Z
     --shape {hexagonal,cubic,planar}
     --chirality n m
     --relax-tube {yes,no}
     --cell-count X Y Z
     --tube-length L_tube
     --xyz2data
     --add-atomtype ELEMENT ATOM-TYPE
     --generator {tubegen,nanogen}

.. autofunction:: nanogen

Examples
--------

::

    > nanogen --chirality 10 10 --cell-count 1 1 5

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import argparse
import importlib
import sys

from ..structure_io.atoms import StructureAtom as Atom
from ..nanogen import TubeGen, tubegen_format_ext_map
from ..structure_io import XYZ2DATAConverter
from ..tools.refdata import CCbond, dVDW

__all__ = ['nanogen', 'tubegen']


def argparser():
    parser = argparse.ArgumentParser()

    parser.add_argument('--element1', default='C',
                        help='element symbol or atomic number of basis '
                        'atom 1 (default: %(default)s)')
    parser.add_argument('--element2', default='C',
                        help='element symbol or atomic number of basis '
                        'atom 2 (default: %(default)s)')
    parser.add_argument('--bond', type=float, default=CCbond,
                        help='Bond length between nearest neighbor atoms. '
                        'Must in units of Angstroms. (default: %(default)s)')

    subparsers = \
        parser.add_subparsers(title='sub-commands',
                              description='StructureGenerator classes',
                              dest='generator_class')

    nanogen_parent_parser = argparse.ArgumentParser(add_help=False)
    nanogen_parent_parser.add_argument('--fname', help='structure file name')
    nanogen_parent_parser.add_argument('--structure-format', default='xyz',
                                       choices=('data', 'xyz'),
                                       help='structure file format '
                                       '(default: %(default)s)')
    nanogen_parent_parser.add_argument('--verbose', action='store_true',
                                       help='verbose output')

    nanogen_graphene_parsers = argparse.ArgumentParser(add_help=False)
    nanogen_graphene_parsers.add_argument('length', type=float,
                                          help='length of graphene sheet in '
                                          '**nanometers** '
                                          '(default: %(default)s)')
    nanogen_graphene_parsers.add_argument('width', type=float,
                                          help='width of graphene '
                                          'sheet in **nanometers** '
                                          '(default: %(default)s)')
    nanogen_graphene_parsers.add_argument('--edge', choices=('AC', 'ZZ'),
                                          help='edge along the `length` '
                                          'of sheet (default: %(default)s)')
    nanogen_graphene_parsers.add_argument('--layer-spacing',
                                          type=float, default=3.35,
                                          help='distance between layers '
                                          'in **Angstroms**. '
                                          '(default: %(default)s)')
    nanogen_graphene_parsers.add_argument('--stacking-order',
                                          choices=('AA', 'AB'), default='AB',
                                          help='Stacking order of '
                                          'graphene layers '
                                          '(default: %(default)s)')

    nanogen_graphene_parser = \
        subparsers.add_parser('GrapheneGenerator',
                              parents=[nanogen_parent_parser,
                                       nanogen_graphene_parsers])
    nanogen_graphene_parser.add_argument('--nlayers', type=int, default=1,
                                         help='Number of graphene layers. '
                                         '(default: %(default)s)')
    nanogen_graphene_parser.set_defaults(func=nanogen)

    nanogen_bilayer_graphene_parser = \
        subparsers.add_parser('BiLayerGrapheneGenerator',
                              parents=[nanogen_parent_parser,
                                       nanogen_graphene_parsers])
    nanogen_bilayer_graphene_parser.add_argument('--layer-rotation-angle',
                                                 type=float,
                                                 help='Rotation angle of '
                                                 'second layer specified in '
                                                 'degrees. '
                                                 '(default: %(default)s)')
    nanogen_bilayer_graphene_parser.set_defaults(func=nanogen)

    nanogen_nanotube_parsers = argparse.ArgumentParser(add_help=False)
    nanogen_nanotube_parsers.add_argument('--n', type=int, default=10,
                                          help='Chiral index `n`'
                                          '(default: %(default)s)')
    nanogen_nanotube_parsers.add_argument('--m', type=int, default=10,
                                          help='Chiral index `m`'
                                          '(default: %(default)s)')
    nanogen_nanotube_parsers.add_argument('--nx', type=int, default=1,
                                          help='Number of repeat unit cells '
                                          'along `x` axis '
                                          '(default: %(default)s)')
    nanogen_nanotube_parsers.add_argument('--ny', type=int, default=1,
                                          help='Number of repeat unit cells '
                                          'along `y` axis '
                                          '(default: %(default)s)')
    nanogen_nanotube_parsers.add_argument('--nz', type=int, default=1,
                                          help='Number of repeat unit cells '
                                          'along `z` axis '
                                          '(default: %(default)s)')
    nanogen_nanotube_parsers.add_argument('--fix-Lz', action='store_true',
                                          help='Generate the nanotube with '
                                          'length as close to the specified '
                                          '`Lz` as possible. If `True`, then '
                                          'non integer `nz` cells are '
                                          'permitted. (default: '
                                          '%(default)s)')
    nanogen_nanotube_parser = \
        subparsers.add_parser('NanotubeGenerator',
                              parents=[nanogen_parent_parser,
                                       nanogen_nanotube_parsers])
    nanogen_nanotube_parser.set_defaults(func=nanogen)

    nanogen_unrolled_nanotube_parser = \
        subparsers.add_parser('UnrolledNanotubeGenerator',
                              parents=[nanogen_parent_parser,
                                       nanogen_nanotube_parsers])
    nanogen_unrolled_nanotube_parser.set_defaults(func=nanogen)

    nanogen_mwnt_parser = \
        subparsers.add_parser('MWNTGenerator',
                              parents=[nanogen_parent_parser,
                                       nanogen_nanotube_parsers])
    nanogen_mwnt_parser.set_defaults(func=nanogen)

    nanogen_nanotube_bundle_parser = \
        subparsers.add_parser('NanotubeBundleGenerator',
                              parents=[nanogen_parent_parser,
                                       nanogen_nanotube_parsers])
    nanogen_nanotube_bundle_parser.add_argument('--vdw-spacing', type=float,
                                                default=dVDW,
                                                help='van der Waals '
                                                'distance between nearest '
                                                'neighbor tubes '
                                                '(default: %(default)s)')
    nanogen_nanotube_bundle_parser.set_defaults(func=nanogen)

    nanogen_mwnt_bundle_parser = \
        subparsers.add_parser('MWNTBundleGenerator',
                              parents=[nanogen_parent_parser,
                                       nanogen_nanotube_parsers])
    nanogen_mwnt_bundle_parser.set_defaults(func=nanogen)

    tubegen_parser = subparsers.add_parser('TubeGen')
    tubegen_parser.add_argument('--chirality', nargs=2, type=int,
                                metavar=('n', 'm'), default=(10, 10),
                                help='(n, m) chiral indices '
                                '(default: %(default)s)')
    parser_group = tubegen_parser.add_mutually_exclusive_group()
    parser_group.add_argument('--cell-count', nargs=3, type=int,
                              metavar=('X', 'Y', 'Z'), default=(1, 1, 1),
                              help='number of unit cells '
                              '(default: %(default)s)')
    parser_group.add_argument('--Lz', type=float, default=None,
                              help='Length of nanotube in units of '
                              '**nanometers** (default: %(default)s)')

    tubegen_parser.add_argument('--format', default='xyz', dest='fmt',
                                choices=tubegen_format_ext_map.keys(),
                                help='output file format '
                                '(default: %(default)s)')
    tubegen_parser.add_argument('--units', default='angstrom',
                                choices=('angstrom', 'bohr'),
                                help='units (default: %(default)s)')
    tubegen_parser.add_argument('--gutter', nargs=3, type=float,
                                metavar=('X', 'Y', 'Z'),
                                default=(1.6735, 1.6735, 0),
                                help='cell gutter (default: %(default)s)')
    tubegen_parser.add_argument('--shape', default='hexagonal',
                                choices=('hexagonal', 'cubic', 'planar'),
                                help='crystal structure '
                                '(default: %(default)s)')
    tubegen_parser.add_argument('--relax-tube', default='yes',
                                choices=('yes', 'no'),
                                help='select tube relaxation '
                                '(default: %(default)s)')

    tubegen_parser.add_argument('--add-atomtype', nargs=2, action='append',
                                metavar=('ELEMENT', 'ATOM-TYPE'),
                                dest='new_atomtypes',
                                help='add new atom to structure data.')
    tubegen_parser.add_argument('--xyz2data', action='store_true',
                                help='convert xyz to LAMMPS data format')
    tubegen_parser.set_defaults(func=tubegen)

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
    generator = getattr(importlib.import_module('sknano.nanogen'),
                        generator_class)
    del kwargs['func']
    generator(**kwargs).save_data(fname=fname,
                                  structure_format=structure_format)


def tubegen(fmt='xyz', units='angstrom', element1='C', element2='C',
            bond=1.421, gutter=(1.6735, 1.6735, 0), shape='hexagonal',
            chirality=(10, 10), relax_tube='yes', cell_count=(1,1,1),
            Lz=None, xyz2data=False, new_atomtypes=None, **kwargs):
    """Generate nano-structure data using `TubeGen` class.

    Parameters
    ----------
    fmt : str, optional
        output file format
    units : str, optional
    element1 : str, optional
        element symbol or atomic number of basis atom 1
    element2 : str, optional
        element symbol or atomic number of basis atom 2
    bond : float, optional
        element1-element2 bond length in ``units``
    gutter : sequence, optional
        cell gutter
    shape : {'hexagonal', 'cubic', 'planar'}, optional
        crystal structure
    chirality : sequence, optional
        chiral indicies (n, m)
    relax_tube : {'yes', 'no'}, optional
    cell_count : sequence, optional
        number of unit cells
    Lz : {None, float}, optional
        tube length in **nanometers**
    xyz2data : bool, optional
        convert xyz to LAMMPS data format
    new_atomtypes : {None, sequence}, optional
        add new (element, atomtype) to LAMMPS data

    """
    if shape == 'planar':
        gutter = (0.0, gutter[1], gutter[2])

    tubegen = TubeGen(fmt=fmt, units=units, element1=element1,
                      element2=element2, bond=bond, gutter=gutter,
                      shape=shape, chirality=chirality,
                      cell_count=cell_count, Lz=Lz,
                      relax_tube=relax_tube)

    tubegen.generate()
    tubegen.cleanup()

    if Lz is not None:
        tubegen.fix_length()

    if xyz2data:
        if fmt == 'xyz':
            xyzconverter = XYZ2DATAConverter(tubegen.output, pad_box=True)
            if new_atomtypes is not None:
                for element, atomtype in new_atomtypes:
                    atom = Atom(element, atomtype=int(atomtype))
                    xyzconverter.add_atomtype(atom=atom)
            xyzconverter.convert()
        else:
            print("can't convert {} format to data. ".format(fmt) +
                  "Must be xyz format.")


def main():

    args = argparser().parse_args()
    args.func(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
