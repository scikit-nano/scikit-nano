#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
======================================================
NanoGen GUI CLI (:mod:`sknano.scripts.nanogenui`)
======================================================

.. currentmodule:: sknano.scripts.nanogenui

"""
from __future__ import absolute_import, print_function, division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import sys

try:
    from sknano.apps.nanogen_gui import NanoGenController, NanoGenModel
except ImportError as e:
    print(e)


__all__ = ['NanoGen']


class NanoGen:
    """Base class for instantiating the NanoGen MVC.

    .. versionadded:: 0.2.24

    .. seealso:: CLI module :py:mod:`sknano.scripts.nanogen`

    """
    def __init__(self, args):
        NanoGenController(args, model=NanoGenModel())


def main():
    NanoGen(sys.argv)

if __name__ == "__main__":
    sys.exit(main())
