#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
========================================================================
Structure analysis script (:mod:`sknano.scripts.analyze_structure`)
========================================================================

.. currentmodule:: sknano.scripts.analyze_structure

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from collections import Counter  # , OrderedDict
import argparse
import os
import sys

# import numpy as np

from sknano.core import listdir_dirnames
from sknano.io import StructureReader
# from sknano.structures import get_chiral_indices_from_str
from ._parser import add_default_arguments

__all__ = ['analyze_structure']


def argparser():
    parser = argparse.ArgumentParser()

    # parser.add_argument('--cfg', '--config', dest='analysis_config',
    #                     metavar='CONFIG_FILE', default=None,
    #                     help='config file with analysis '
    #                     'settings. (default: %(default)s)')
    # parser.add_argument(
    #     '--structure-format', default=None, choices=('data', 'dump', 'xyz'),
    #     help='input structure data file format. If `None`, then guess '
    #     'format from file name. (default: %(default)s)')

    # subparsers = parser.add_subparsers(title='sub-commands')

    # poav_parser = subparsers.add_parser('POAV')
    # poav_parser.add_argument('--atom-ids', metavar='ATOM_ID', type=int,
    #                          nargs='+',
    #                          help='One or more atom IDs to analyze. '
    #                          '(default: %(default)s)')
    # poav_parser.add_argument('--include-NN', action='store_true',
    #                          help='analyze nearest neighbor atoms for '
    #                          'each atom in list of `atom-ids`. '
    #                          '(default: %(default)s)')
    # poav_parser.set_defaults(analyze='POAV')

    # parser.add_argument('structure_file',
    #                     help='input structure data files.')

    parser.add_argument('--rootdir', default=os.curdir,
                        help='root data dir (default: %(default)s)')
    parser.add_argument('--structure-file', default='system.dump',
                        help='structure file name (default: %(default)s)')
    parser = add_default_arguments(parser)

    return parser


def analyze_structure(rootdir=None, structure_file=None, **kwargs):
    if rootdir is None or not os.path.isdir(rootdir):
        rootdir = os.curdir
    rundirs = \
        listdir_dirnames(rootdir,
                         filterfunc=lambda name: name.startswith('run_'))
    ion_coordination_counter = Counter()
    graphene_coordination_counter = Counter()
    for rundir in rundirs:
        dumpdata = StructureReader.read(os.path.join(rundir, structure_file))
        dumpdata.remap_atomattr_names({'c_atom_ke': 'ke'})
        trj = dumpdata.trajectory
        final_snapshot = trj[-1]
        atoms = final_snapshot.atoms
        atoms.mapatomattr('type', 'element', {1: 'C', 2: 'N'})
        atoms.update_attrs()
        ion = atoms.filtered(atoms.elements == 'N')[0]
        graphene = atoms.filtered(atoms.elements == 'C')
        ion_coordination_counter.update([ion.CN])
        graphene_coordination_counter.update(graphene.coordination_numbers)
        statstr = 'run statistics:\n'
        statstr += 'rundir: {}\n'.format(rundir)
        statstr += 'ion CN: {}\n'.format(ion.CN)
        statstr += 'graphene CN counts: {}\n'.format(
            Counter(graphene.coordination_numbers))
        print(statstr)

    print('final stats')
    print('ion_coordination_counter: {}'.format(ion_coordination_counter))
    print('graphene_coordination_counter: {}'.format(
          graphene_coordination_counter))


def main():
    args = argparser().parse_args()
    analyze_structure(**vars(args))

if __name__ == '__main__':
    sys.main(main())
