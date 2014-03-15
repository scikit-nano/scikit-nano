# -*- coding: utf-8 -*-
"""
==========================================================
Tools for structure data I/O (:mod:`sknano.structure_io`)
==========================================================

.. currentmodule:: sknano.structure_io

Contents
========

Classes for handling `LAMMPS data` structure data format
--------------------------------------------------------
.. autosummary::
   :toctree: generated/

   DATAData
   DATAReader
   DATAWriter
   DATA2XYZConverter
   DATAFormat
   DATAError

Classes for handling `xyz` structure data format
-------------------------------------------------
.. autosummary::
   :toctree: generated/

   XYZData
   XYZReader
   XYZWriter
   XYZ2DATAConverter
   XYZFormat

Base classes for creating new structure data readers/writers/converters
------------------------------------------------------------------------
.. autosummary::
   :toctree: generated/

   StructureData
   StructureReader
   StructureWriter
   StructureConverter
   StructureFormat

Base custom exception classes for handling errors
--------------------------------------------------
.. autosummary::
   :toctree: generated/

   StructureDataError
   StructureReaderError
   StructureWriterError
   StructureConverterError
   StructureFormatError

Sub-packages
=============
.. autosummary::
   :toctree: generated/

   atoms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._lammps_data_format import *
from ._xyz_format import *
from ._structure_data import *

__all__ = [s for s in dir() if not s.startswith('_')]
