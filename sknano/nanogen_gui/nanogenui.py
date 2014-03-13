# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
======================================================
NanoGen GUI CLI (:mod:`sknano.nanogen_gui.nanogenui`)
======================================================

.. currentmodule:: sknano.nanogen_gui.nanogenui

"""
from __future__ import absolute_import, print_function, division

import sys

from ._ng_controller import NGController
from ._ng_model import NGModel

__all__ = ['NanoGen']


class NanoGen(object):
    """Base class for instantiating the NanoGen MVC."""
    def __init__(self, args):
        model = NGModel()
        NGController(args, model)


def main():
    NanoGen(sys.argv)

if __name__ == "__main__":
    sys.exit(main())
