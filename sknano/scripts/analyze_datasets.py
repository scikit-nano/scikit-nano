#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
================================================================
Command line script (:mod:`sknano.scripts.analyze_datasets`)
================================================================

.. currentmodule:: sknano.scripts.analyze_datasets

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
import argparse
import os
import sys

import h5py
#import tables as tbl
import numpy as np

from sknano.io import DATAReader
from sknano.structures import get_chiral_indices_from_str

__all__ = ['analyze_datasets']


def run_analysis(datafile, datalog=None, datagroup=None, Ch_subgroup=None):
    print('data file: {}'.format(datafile))
    structure_Ch = get_chiral_indices_from_str(datafile[:4])
    if structure_Ch is not None and datalog is not None:
        atoms = DATAReader(datafile).atoms
        Ch = str(structure_Ch)
        if datagroup is None:
            datagroup = '0'
        group = datalog.require_group(datagroup)

        gCh = group.require_group(Ch)
        #subgCh = gCh.require_group('generated')

        dset = 'mean_misalignment_angle'
        dset = gCh.require_dataset(dset, shape=(atoms.Natoms,),
                                   dtype=np.float64)
        mean_misalignment_angle = np.degrees(atoms.mean_misalignment_angle)
        print('mean_misalignment_angle: {}'.format(mean_misalignment_angle))


def analyze_datasets(datalist=None, log_data=False, datagroup=None,
                     hdf5_file=None):
    if not isinstance(datalist, list):
        raise TypeError('Expected a list of structure data files '
                        'and/or directories.')

    logkwargs = dict(datagroup=datagroup)

    datalog = None
    if log_data:
        if hdf5_file is None:
            hdf5_file = 'tmp.hdf5'
        datalog = h5py.File(hdf5_file)

    logkwargs['datalog'] = datalog

    for data in datalist:
        if os.path.isfile(data):
            run_analysis(data, **logkwargs)
        elif os.path.isdir(data):
            for dirpath, _, fnames in os.walk(data):
                fnames = [fname for fname in fnames if
                          fname.endswith(('data'))]

                for fname in fnames:
                    run_analysis(fname, **logkwargs)
    try:
        datalog.close()
    except AttributeError:
        pass


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('datalist', nargs='+', metavar='data',
                        help='List of structure data files and/or '
                        'directories.')
    parser.add_argument('--log-data', action='store_true',
                        help='log data to file.')
    parser.add_argument('--hdf5-file', default=None,
                        help='hdf5 file')
    parser.add_argument('--datagroup', default=None,
                        help='hdf5 data group to log data to')

    return parser


def main():
    args = argparser().parse_args()
    analyze_datasets(**vars(args))

if __name__ == '__main__':
    sys.main(main())
