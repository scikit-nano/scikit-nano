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

 analysis       --- functions and classes for structure analysis
 nanogen        --- classes for generating nanostructures
 nanogen_gui    --- NanoGen GUI front-end
 scripts        --- command-line scripts
 structure_io   --- classes for structure data I/O
 tools          --- helper funcs, LUTs, etc.

Utilitiy tools
--------------

::

 __version__    --- sknano version string

"""
from __future__ import absolute_import, print_function

__all__ = ['analysis',
           'nanogen',
           'nanogen_gui',
           'scripts',
           'structure_io',
           'tools']

from sknano.version import version as __version__
