# -*- coding: utf-8 -*-
"""
===============================================================
Package for generating nano-structures (:mod:`sknano.nanogen`)
===============================================================

.. currentmodule:: sknano.nanogen

sknano package for generating nano-structures including graphene,
bi-layer graphene, :math:`N`-layer graphene,
single-walled nanotubes (:abbr:`SWNTs (single-walled nanotubes)`),
multi-walled nanotubes (:abbr:`MWNTs (multi-walled nanotubes)`),
and bundles of :abbr:`SWNTs` and :abbr:`MWNTs`.

.. note::
   The default basis atoms are both carbon and therefore the default
   nano-structures are carbon nano-structures. However,
   you can use any atom in the periodic table of elements as a basis atom.

.. versionadded:: 0.2.8
   `MWNTGenerator` implemented

.. versionadded:: 0.2.6
   `GrapheneVacancyGenerator` and `NanotubeVacancyGenerator` implemented

.. versionadded:: 0.2.4
   `NanotubeBundleGenerator` implemented

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

   StructureGenerator

   GrapheneGenerator
   BiLayerGrapheneGenerator

   NanotubeGenerator
   NanotubeBundleGenerator
   MWNTGenerator

   TubeGen

.. seealso:: CLI module :py:mod:`sknano.scripts.nanogen`

Tools for manipulating nano-structures
--------------------------------------

.. autosummary::
   :toctree: generated/

   VacancyGenerator
   GrapheneVacancyGenerator
   NanotubeVacancyGenerator

Custom exception classes for handling errors
--------------------------------------------

.. autosummary::
   :toctree: generated/

   StructureGeneratorError
   GrapheneGeneratorError
   NanotubeGeneratorError
   VacancyGeneratorError

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from ._defect_generator import *
from ._fullerene import *
from ._fullerene_generator import *
from ._graphene_defect_generator import *
from ._graphene_generator import *
from ._graphene import *
from ._nanotube_defect_generator import *
from ._nanotube_generator import *
from ._nanotube_junction_generator import *
from ._nanotube import *
from ._structure_generator import *
from ._tubegen import *
from ._twisted_bundle_generator import *

__all__ = [s for s in dir() if not s.startswith('_')]
