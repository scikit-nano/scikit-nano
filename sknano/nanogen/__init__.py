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

.. versionadded:: 0.2.23
   `UnrolledNanotubeGenerator` implemented.

.. versionadded:: 0.2.20
   `MWNTBundleGenerator` implemented.

.. versionchanged:: 0.2.20
   `MWNTGenerator` no longer generates MWNT *bundles*, only *single* MWNTs.
   To generate bundled MWNT structure data, use the `MWNTBundleGenerator`
   class.

.. versionadded:: 0.2.8
   `MWNTGenerator` implemented.

.. versionadded:: 0.2.6
   `GrapheneVacancyGenerator` and `NanotubeVacancyGenerator` implemented.

.. versionadded:: 0.2.4
   `NanotubeBundleGenerator` implemented.

.. seealso:: CLI module :py:mod:`sknano.scripts.nanogen`

Contents
========

Classes for creating abstract physical representations of nano-structures
--------------------------------------------------------------------------
.. autosummary::
   :toctree: generated/

   Fullerene
   Graphene
   GraphenePrimitiveCell
   Nanotube
   NanotubeBundle

Classes for generating fullerene structure data
------------------------------------------------
.. autosummary::
   :toctree: generated/

   FullereneGenerator
   FullereneBravaisLatticeGenerator

Classes for generating graphene structure data
-----------------------------------------------
.. autosummary::
   :toctree: generated/

   GrapheneGenerator
   BiLayerGrapheneGenerator
   UnrolledNanotubeGenerator

Classes for generating nanotube structure data
-----------------------------------------------
.. autosummary::
   :toctree: generated/

   NanotubeGenerator
   NanotubeBundleGenerator
   UnrolledNanotubeGenerator
   MWNTGenerator
   MWNTBundleGenerator
   TubeGen

Classes for generating nano-structure data with vacancies
----------------------------------------------------------
.. autosummary::
   :toctree: generated/

   GrapheneVacancyGenerator
   NanotubeVacancyGenerator

Base classes for creating new structure data generators
--------------------------------------------------------
.. autosummary::
   :toctree: generated/

   StructureGenerator
   VacancyGenerator

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from ._defect_generators import *
from ._fullerene_bravais_lattice_generators import *
from ._fullerene_generators import *
from ._fullerenes import *
from ._graphene_defect_generators import *
from ._graphene_generators import *
from ._graphene import *
from ._nanotube_defect_generators import *
from ._nanotube_bundle_generators import *
from ._nanotube_generators import *
from ._nanotube_junction_generators import *
from ._nanotubes import *
from ._structure_generator import *
from ._tubegen import *
from ._twisted_bundle_generators import *

__all__ = [s for s in dir() if not s.startswith('_')]
