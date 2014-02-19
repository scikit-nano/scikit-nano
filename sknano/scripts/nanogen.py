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

import argparse
import sys

from pkshared.tools.refdata import CCbond

from ..chemistry import Atom
from ..nanogen import GrapheneGenerator, NanotubeBundleGenerator, \
    TubeGen, format_ext
from ..structure_io import XYZ2DATAConverter

__all__ = ['nanogen']


def _argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--format', default='xyz', dest='fmt',
                        choices=tuple([fmt for fmt in format_ext.iterkeys()]),
                        help='output file format (default: %(default)s)')
    parser.add_argument('--units', default='angstrom',
                        choices=('angstrom', 'bohr'),
                        help='units (default: %(default)s)')
    parser.add_argument('--element1', default='C',
                        help='element symbol or atomic number of basis '
                        'atom 1 (default: %(default)s)')
    parser.add_argument('--element2', default='C',
                        help='element symbol or atomic number of basis '
                        'atom 2 (default: %(default)s)')
    parser.add_argument('--bond', type=float, default=CCbond,
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

    parser.add_argument('--xyz2data', action='store_true',
                        help='convert xyz to LAMMPS data format')
    parser.add_argument('--add-atomtype', nargs=2, action='append',
                        metavar=('ELEMENT', 'ATOM-TYPE'),
                        dest='new_atomtypes',
                        help='add new atom to structure data.')

    parser.add_argument('--generator', choices=('tubegen', 'nanogen'),
                        default='tubegen', help='nano-structure generator '
                        '(default: %(default)s)')

    return parser


def nanogen(fmt='xyz', units='angstrom', element1='C', element2='C',
            bond=1.421, gutter=(1.6735, 1.6735, 0), shape='hexagonal',
            chirality=(10, 10), relax_tube='yes', cell_count=(1,1,1),
            tube_length=None, xyz2data=False, new_atomtypes=None,
            generator='tubegen'):
    """Generate nano-structure.

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
    tube_length : {None, float}, optional
        tube length in **nanometers**
    xyz2data : bool, optional
        convert xyz to LAMMPS data format
    new_atomtypes : {None, sequence}, optional
        add new (element, atomtype) to LAMMPS data
    generator : {'tubegen', 'nanogen'}, optional
        nano-structure generator

    """
    if generator == 'tubegen':
        tubegen = TubeGen(fmt=fmt, units=units, element1=element1,
                          element2=element2, bond=bond, gutter=gutter,
                          shape=shape, chirality=chirality,
                          cell_count=cell_count, tube_length=tube_length,
                          relax_tube=relax_tube)

        tubegen.generate()
        tubegen.cleanup()

        if tube_length is not None:
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
    else:
        n, m = chirality
        nx, ny, nz = cell_count

        if shape == 'planar':
            generator = GrapheneGenerator()
        else:
            generator = NanotubeBundleGenerator(n, m, nx=nx, ny=ny, nz=nz,
                                                Lz=tube_length,
                                                bundle_geometry=shape)

        if fmt not in ('xyz', 'data'):
            fmt = 'xyz'

        generator.save_data(structure_format=fmt)


def main():

    args = _argparser().parse_args()
    nanogen(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
