# -*- coding: utf-8 -*-
"""
==========================================================
Tools for structure data I/O (:mod:`sknano.structure_io`)
==========================================================

.. currentmodule:: sknano.structure_io

Contents
========

Classes for creating abstract structure data objects for reading/writing data
------------------------------------------------------------------------------

.. autosummary::
   :toctree: generated/

   StructureData
   LAMMPSDATA
   XYZDATA

Classes defining format specs/properties
----------------------------------------

.. autosummary::
   :toctree: generated/

   StructureSpecs
   LAMMPSDATASpecs

Classes for reading structure data
----------------------------------

.. autosummary::
   :toctree: generated/

   DATAReader
   XYZReader

Classes for saving structure data
---------------------------------

.. autosummary::
   :toctree: generated/

   DATAWriter
   XYZWriter

Classes for converting structure data formats
---------------------------------------------

.. autosummary::
   :toctree: generated/

   DATA2XYZConverter
   XYZ2DATAConverter

Custom exception classes for handling errors
--------------------------------------------

.. autosummary::
   :toctree: generated/

   StructureReaderError
   StructureInputError

   StructureWriterError
   StructureOutputError

   StructureConverterError
   StructureFormatError

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from ._structure_converters import *
from ._structure_data import *
from ._structure_specs import *
from ._lammps_data_structure_data import *
from ._xyz_structure_data import *

__all__ = [s for s in dir() if not s.startswith('_')]
