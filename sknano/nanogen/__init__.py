# -*- coding: utf-8 -*-
"""
======================================================================
Package for generating nano-structures (:mod:`sknano.nanogen`)
======================================================================

.. currentmodule:: sknano.nanogen

sknano package for generating nano-structures including graphene and
nanotubes.

Contents
========

Abstract nanostructure objects
------------------------------

.. autosummary::
   :toctree: generated/

   Graphene
   Nanotube
   NanotubeBundle

Structure generator classes
---------------------------

.. autosummary::
   :toctree: generated/

   GrapheneGenerator
   BiLayerGrapheneGenerator

   NanotubeGenerator
   NanotubeBundleGenerator
   SWNTGenerator
   MWNTGenerator

   TubeGen

.. note::
   The :py:class:`~sknano.nanogen.NanotubeGenerator` class
   does not yet generate nanotube *bundles*. Only single tubes.

.. seealso:: CLI module :py:mod:`sknano.scripts.nanogen`

Classes for changing structures
-------------------------------

.. autosummary::
   :toctree: generated/

   VacancyGenerator

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from .graphene import *
from .nanotube import *
from .tubegen import *
from .vacancygenerator import *

__all__ = [s for s in dir() if not s.startswith('_')]
