# -*- coding: utf-8 -*-
"""
===========================================================================
Base argparser for comand line scripts (:mod:`sknano.scripts._parser`)
===========================================================================

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import argparse

from sknano.version import release, full_version, git_revision

__all__ = ['add_default_arguments', 'base_parser']


def add_default_arguments(parser):

    version = full_version
    if not release:
        version = '-'.join((full_version, git_revision[:7]))

    parser.add_argument('--debug', action='store_true', help='debug output')
    parser.add_argument('--verbose', action='store_true',
                        help='verbose output')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(version),
                        help="show %(prog)s's version number and exit")

    return parser


def base_parser():
    parser = add_default_arguments(argparse.ArgumentParser())
    return parser
