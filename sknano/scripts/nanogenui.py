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


def main():
    try:
        from sknano.apps.nanogen_gui import NanoGenController, NanoGenModel
        NanoGenController(sys.argv, model=NanoGenModel())
    except ImportError:
        print('PyQt4 or PyQt5 required to run NanoGen GUI')


if __name__ == "__main__":
    sys.exit(main())
