# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
======================================================
NanoGen GUI CLI (:mod:`sknano.scripts.nanogenui`)
======================================================

.. currentmodule:: sknano.scripts.nanogenui

"""
from __future__ import absolute_import, print_function, division
__docformat__ = 'restructuredtext en'

import sys

from ..nanogen_gui import NGController, NGModel

__all__ = ['NanoGen']


class NanoGen(object):
    """Base class for instantiating the NanoGen MVC.

    .. versionadded:: 0.2.24

    .. seealso:: CLI module :py:mod:`sknano.scripts.nanogen`

    """
    def __init__(self, args):
        model = NGModel()
        NGController(args, model)


def main():
    NanoGen(sys.argv)

if __name__ == "__main__":
    sys.exit(main())
