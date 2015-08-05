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

from sknano.version import version

__all__ = ['add_default_arguments', 'base_parser']


def add_default_arguments(parser):
    """Add a set of default arguments to an instance of \
        :class:`~python:argparse.ArgumentParser`.

    Parameters
    ----------
    parser : :class:`~python:argparse.ArgumentParser`

    Returns
    -------
    parser : class:`~python:argparse.ArgumentParser`

    """
    parser.add_argument('--debug', action='store_true', help='debug output')
    parser.add_argument('--verbose', action='store_true',
                        help='verbose output')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(version),
                        help="show %(prog)s's version number and exit")

    return parser


base_parser = add_default_arguments(argparse.ArgumentParser())
