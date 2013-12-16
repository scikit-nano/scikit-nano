# -*- coding: utf-8 -*-
"""
======================================================================
Package for generating nano-structures (:mod:`sknano.nanogen`)
======================================================================

.. currentmodule:: sknano.nanogen

sknano package for generating nano-structures including graphene,
bi-layer graphene, :math:`N`-layer graphene, single-walled carbon nanotubes
and single-walled carbon nanotube bundles.

.. versionchanged:: 0.2.4
   NanotubeBundleGenerator implemented

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
