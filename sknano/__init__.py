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

 chemistry      --- abstract data structures for chemistry
 nanogen        --- modules for structure analysis and structure generation
 scripts        --- command-line scripts
 structure_io   --- classes for structure data I/O
 tools          --- tools for analysis, helper funcs, LUTs, etc.

Utilitiy tools
--------------

::

 __version__    --- sknano version string

"""
from __future__ import print_function, absolute_import

__all__ = ['chemistry', 'nanogen', 'scripts', 'structure_io', 'tools']

from sknano.version import version as __version__


def get_path():
    try:
        from pkthemes import get_path
        return get_path()
    except ImportError:
        return None
