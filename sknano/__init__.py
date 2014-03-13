# -*- coding: utf-8 -*-
"""
============================================================================
scikit-nano: Python toolkit for generating and analyzing nanostructure data
============================================================================

Contents
========

Subpackages
-----------

::

 chemistry      --- abstract data structures for chemistry
 nanogen        --- modules for structure analysis and structure generation
 nanogen_gui    --- NanoGen GUI front-end
 scripts        --- command-line scripts
 structure_io   --- classes for structure data I/O
 tools          --- tools for analysis, helper funcs, LUTs, etc.

Utilitiy tools
--------------

::

 __version__    --- sknano version string

"""
from __future__ import print_function, absolute_import

__all__ = ['chemistry',
           'nanogen',
           'nanogen_gui',
           'scripts',
           'structure_io',
           'tools']

from sknano.version import version as __version__
