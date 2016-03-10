# -*- coding: utf-8 -*-
"""
====================================================
PDB format (:mod:`sknano.io._pdb`)
====================================================

.. currentmodule:: sknano.io._pdb

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import OrderedDict

from monty.io import zopen
from pyparsing import ParseException

from sknano.core import get_fpath
from sknano.core.atoms import Atom
from .tokenizers import PDBTokenizer
from ._base import StructureIO, StructureIOError, StructureFormatSpec, \
    default_comment_line

__all__ = ['PDBData', 'PDBReader', 'PDBWriter', 'PDBFormatSpec']


class PDBReader(StructureIO):
    """Class for reading pdb chemical file format.

    Parameters
    ----------
    fpath : str
        pdb structure file

    """
    def __init__(self, fpath):
        super().__init__(fpath=fpath)

        self.records = OrderedDict()
        self.tokenizer = PDBTokenizer()

        if self.fpath is not None:
            self.read()

    def read(self):
        """Read PDB file."""
        tokenizer = self.tokenizer
        self.structure.clear()
        with zopen(self.fpath) as f:
            while True:
                line = f.readline().strip()
                record = line[:6].strip()
                if record == 'END':
                    break

                parse_method = '_parse_{}'.format(record.lower())
                if not hasattr(tokenizer, parse_method):
                    print('No parser for {} records'.format(record))
                    continue

                try:
                    record = getattr(tokenizer, parse_method)(line)
                except ParseException as e:
                    print('Error in PDB file format')
                    raise StructureIOError(e)
                else:
                    if isinstance(record, Atom):
                        self.structure.append(record)


class PDBWriter:
    """Class for writing pdb chemical file format."""

    @classmethod
    def write(cls, fname=None, outpath=None, fpath=None, structure=None,
              atoms=None, **kwargs):
        """Write structure data to file.

        Parameters
        ----------
        fname : str, optional
            Output file name.
        outpath : str, optional
            Output file path.
        fpath : str, optional
            Full path (directory path + file name) to output data file.
        atoms : :py:class:`~sknano.io.atoms.Atoms`
            An :py:class:`~sknano.io.atoms.Atoms` instance.

        """
        if structure is None and atoms is None:
            raise ValueError('Expected either `structure` or `atoms` object.')

        if structure is not None and atoms is None:
            atoms = structure.atoms

        if fpath is None:
            fpath = get_fpath(fname=fname, ext='pdb', outpath=outpath,
                              overwrite=True, add_fnum=False)

        pdb = PDBData()
        pdb.write(pdbfile=fpath, atoms=atoms, **kwargs)

        with zopen(fpath, 'wt') as f:
            for atom in atoms:
                f.write(PDBFormatter.format(atom))


class PDBData(PDBReader):
    """Class for reading and writing structure data in PDB data format.

    Parameters
    ----------
    fpath : str, optional

    """
    def __init__(self, fpath=None):
        try:
            # super(PDBData, self).__init__(fpath=fpath)
            super().__init__(fpath=fpath)
        except StructureIOError:
            pass

    def write(self, pdbfile=None, atoms=None, comment_line=None, **kwargs):
        try:
            kwargs.update(self.kwargs)

            if not pdbfile:
                if self.fpath is None:
                    error_msg = 'Invalid `pdbfile` {}'.format(pdbfile)
                    raise ValueError(error_msg)
                else:
                    pdbfile = self.fpath

            if comment_line is None:
                comment_line = default_comment_line

            if atoms is not None:
                self._atoms = atoms

            super()._update_atoms(**kwargs)

            try:
                with zopen(pdbfile, 'wt') as fh:
                    self._write_header(fh, comment_line)
                    self._write_atoms(fh)
            except OSError as e:
                print(e)

            self._atoms = self._atoms_copy

        except (TypeError, ValueError) as e:
            print(e)


class PDBFormatSpec(StructureFormatSpec):
    """Class defining the structure file format for PDB data.

    Parameters
    ----------

    """
    def __init__(self):
        super().__init__()

        # self._properties['RECORDS'] = records


class PDBFormatter:
    """Formatter class to convert a `PDBAtom` to a formatted string."""

    @classmethod
    def format(cls):
        pass


mandatory_records = ['HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS',
                     'EXPDTA', 'AUTHOR', 'REVDAT', 'REMARK 2',
                     'REMARK 3', 'SEQRES', 'CRYST1', 'ORIGX1',
                     'ORIGX2', 'ORIGX3', 'SCALE1', 'SCALE2', 'SCALE3',
                     'MASTER', 'END']
optional_records = ['OBSLTE', 'SPLIT', 'CAVEAT', 'NUMMDL', 'MDLTYP',
                    'SPRSDE', 'JRNL', 'REMARK 0', 'REMARK 1',
                    'REMARK N', 'DBREF', 'DBREF1', 'DBREF2',
                    'SEQADV', 'MODRES', 'HET', 'HETNAM', 'HETSYN',
                    'FORMUL', 'HELIX', 'SHEET', 'SSBOND', 'LINK',
                    'CISPEP', 'SITE', 'MTRIX1', 'MTRIX2', 'MTRIX3',
                    'MODEL', 'ATOM', 'ANISOU', 'TER', 'HETATM',
                    'ENDMDL', 'CONECT']
