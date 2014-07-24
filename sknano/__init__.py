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

 apps           --- GUI front-end app development
 core           --- core classes and functions
 generators     --- classes for generating abstract object representations
                    of nanostructures
 io             --- classes for reading/writing nanostructure data
 scripts        --- command-line scripts
 structures     ---
 utils          --- utility classes and functions

Utilitiy tools
--------------

::

 __version__    --- sknano version string

"""
from __future__ import absolute_import, division, print_function

__all__ = ['apps',
           'core',
           'generators',
           'io',
           'scripts',
           'structures',
           'utils']

from sknano.version import version as __version__
