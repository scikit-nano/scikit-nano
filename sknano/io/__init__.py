# -*- coding: utf-8 -*-
"""
==========================================================
Structure data I/O (:mod:`sknano.io`)
==========================================================

.. currentmodule:: sknano.io

Contents
========

Base I/O classes to inherit from when creating new I/O classes
----------------------------------------------------------------
.. autosummary::
   :toctree: generated/

   StructureIO
   StructureReader
   StructureWriter
   StructureIOMixin
   StructureReaderMixin
   StructureWriterMixin
   StructureConverter
   StructureFormatSpec
   StructureIOError

I/O classes for the `LAMMPS data` structure data format
--------------------------------------------------------
.. autosummary::
   :toctree: generated/

   DATAIO
   DATAData
   DATAReader
   DATAWriter
   DATAFormatSpec
   DATAIOError

I/O classes for the `LAMMPS dump` structure data format
--------------------------------------------------------
.. autosummary::
   :toctree: generated/

   DUMPIO
   DUMPData
   DUMPReader
   DUMPWriter
   DUMPFormatSpec
   DUMPIOError

I/O classes for the `pdb` structure data format
-------------------------------------------------
.. autosummary::
   :toctree: generated/

   PDBIO
   PDBData
   PDBReader
   PDBWriter
   PDBFormatSpec
   PDBIOError

I/O classes for the `xyz` structure data format
-------------------------------------------------
.. autosummary::
   :toctree: generated/

   XYZIO
   XYZData
   XYZReader
   XYZWriter
   XYZFormatSpec
   XYZIOError

Sub-packages
-------------

* tokenizers (:mod:`sknano.io.tokenizers`)

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._base import *
from ._lammps_data import *
from ._lammps_dump import *
from ._pdb import *
from ._xyz import *

__all__ = [s for s in dir() if not s.startswith('_')]
