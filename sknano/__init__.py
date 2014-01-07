# -*- coding: utf-8 -*-
"""
=====================================
scikit-nano: Toolkit for nano-science
=====================================

Contents
========

Subpackages
-----------

::

 chemistry                --- abstract data structures for chemistry
 nanogen                  --- nano-structure generator
 scripts                  --- command-line scripts
 structure_io             --- helper classes for structure data I/O

Utilitiy tools
--------------

::

 __version__   --- sknano version string

"""
from __future__ import print_function, absolute_import

__all__ = ['chemistry', 'nanogen', 'scripts', 'structure_io']

from sknano.version import version as __version__


def get_path():
    try:
        from pkthemes import get_path
        return get_path()
    except ImportError:
        return None
