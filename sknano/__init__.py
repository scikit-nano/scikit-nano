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

 apps           --- GUI front-end development
 core           --- core classes and functions
 generators     --- classes for generating nanostructure data
 io             --- classes for reading/writing nanostructure data
 scripts        --- command-line utilities
 structures     --- abstract representations of nanostructures
 testing        --- modules for testing
 utils          --- utility modules for analysis, testing, and general use

Utilitiy tools
--------------

::

 __version__    --- sknano version string

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__all__ = ['apps',
           'core',
           'generators',
           'io',
           'scripts',
           'structures',
           'testing',
           'utils']

from sknano.version import version as __version__
