#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
====================================================
Command line script (:mod:`sknano.scripts.sknano`)
====================================================

CLI to :mod:`sknano` tools.

.. currentmodule:: sknano.scripts.sknano

.. code-block:: python

   > sknano --help

.. autofunction:: sknano

Examples
--------


"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import argparse
# import importlib
import sys

from ._parser import add_default_arguments

__all__ = ['sknano', 'sknano_parser']


def sknano_parser():
    """:mod:`~sknano.scripts.sknano` script \
        :class:`~python:argparse.ArgumentParser`."""
    parser = argparse.ArgumentParser()
    parser = add_default_arguments(parser)


def sknano(**kwargs):
    """:mod:`~sknano.scripts.sknano` script function."""
    pass


def main():
    args = sknano_parser().parse_args()
    sknano(**vars(args))


if __name__ == '__main__':
    sys.exit(main())
