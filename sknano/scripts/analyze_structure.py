# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
================================================================
Command line script (:mod:`sknano.scripts.analyze_structure`)
================================================================

.. currentmodule:: sknano.scripts.analyze_structure

"""
from __future__ import print_function, division, absolute_import
import argparse
import os
import sys

__all__ = ['analyze_structure']


def analyze_structure(data):
    if os.path.isfile(data):
        print('data: {}'.format(data))
    elif os.path.isdir(data):
        for dirpath, dirnames, fnames in os.walk(data):
            print('dirpath, dirnames, fnames: {}, {}, {}'.format(
                dirpath, dirnames, fnames))


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('data', help='Structure data file')

    return parser


def main():
    args = argparser().parse_args()
    analyze_structure(**vars(args))

if __name__ == '__main__':
    sys.main(main())
