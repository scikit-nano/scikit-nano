# -*- coding: utf-8 -*-
"""
===============================================================
Package for generating nano-structures (:mod:`sknano.nanogen`)
===============================================================

.. currentmodule:: sknano.nanogen

Collection of modules for generating nano-structures including graphene,
bi-layer graphene, :math:`N`-layer graphene,
single-walled nanotubes (:abbr:`SWNTs (single-walled nanotubes)`),
multi-walled nanotubes (:abbr:`MWNTs (multi-walled nanotubes)`),
and bundles of :abbr:`SWNTs` and :abbr:`MWNTs`.

See the individual class doc pages for detailed documentation and
examples on how to use them.

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

.. seealso::

   CLI script: :mod:`sknano.scripts.nanogen`
      The :mod:`~sknano.scripts.nanogen` script provides a
      :abbr:`CLI (command-line interface)` to the `nanogen` structure
      generator classes.

   GUI front-end: :mod:`sknano.scripts.nanogenui`
      The :mod:`~sknano.scripts.nanogenui` script launches a
      :abbr:`GUI (graphical user interface)` front-end to the
      `nanogen` structure generator classes.

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

Helper functions
----------------
.. autosummary::
   :toctree: generated/

   generate_Ch_list
   generate_Ch_property_grid

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

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
from ._nanogen_funcs import *

__all__ = [s for s in dir() if not s.startswith('_')]
