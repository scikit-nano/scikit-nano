# -*- coding: utf-8 -*-
"""
====================================================
CIF format (:mod:`sknano.io.cif`)
====================================================

.. currentmodule:: sknano.io.cif

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import OrderedDict

from monty.io import zopen
from pyparsing import ParseException

from sknano.core import get_fpath
from sknano.core.atoms import Atom
from .tokenizers import CIFTokenizer
from .base import StructureData, StructureDataError, StructureDataFormatter, \
    default_comment_line

__all__ = ['CIF', 'CIFData', 'CIFReader', 'CIFWriter', 'CIFFormatter',
           'CIFError', 'CIFIO', 'CIFIOReader', 'CIFIOWriter', 'CIFIOFormatter',
           'CIFIOError', 'CIFFormatSpec']


class CIFReader(StructureData):
    """Class for reading cif chemical file format.

    Parameters
    ----------
    fpath : str
        cif structure file

    """
    def __init__(self, fpath, formatter=None, **kwargs):
        if formatter is None or not isinstance(formatter, CIFFormatter):
            formatter = CIFFormatter()
        super().__init__(fpath=fpath, formatter=formatter)

        self.records = OrderedDict()
        self.tokenizer = CIFTokenizer()

        if self.fpath is not None:
            self.read()

    def read(self):
        """Read CIF file."""
        tokenizer = self.tokenizer
        # formatter = self.formatter

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
                    print('Error in CIF file format')
                    raise StructureDataError(e)
                else:
                    if isinstance(record, Atom):
                        self.structure.append(record)

CIFIOReader = CIFReader


class CIFWriter:
    """Class for writing cif chemical file format."""

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
            fpath = get_fpath(fname=fname, ext='cif', outpath=outpath,
                              overwrite=True, add_fnum=False)

        cif = CIFData()
        cif.write(ciffile=fpath, atoms=atoms, **kwargs)

CIFIOWriter = CIFWriter


class CIFData(CIFReader):
    """Class for reading and writing structure data in CIF data format.

    Parameters
    ----------
    fpath : str, optional

    """
    def __init__(self, fpath=None):
        try:
            # super(CIFData, self).__init__(fpath=fpath)
            super().__init__(fpath=fpath)
        except StructureDataError:
            pass

    def write(self, ciffile=None, atoms=None, mode='w',
              comment_line=None, **kwargs):
        """Write cif file.

        Parameters
        ----------
        ciffile : {None, str}, optional
        atoms : :class:`~sknano.core.atoms.Atoms`, optional
        mode : {'w', 'a'}, optional

        """

        try:
            kwargs.update(self.kwargs)

            if not ciffile:
                if self.fpath is None:
                    error_msg = 'Invalid `ciffile` {}'.format(ciffile)
                    raise ValueError(error_msg)
                else:
                    ciffile = self.fpath

            if comment_line is None:
                comment_line = default_comment_line

            if atoms is not None:
                self._atoms = atoms

            super()._update_structure(**kwargs)

            mode += 't'

            try:
                with zopen(ciffile, mode) as stream:
                    self._write_header(stream, comment_line)
                    self._write_atoms(stream)
            except OSError as e:
                print(e)

            super()._reset_structure()

        except (TypeError, ValueError) as e:
            print(e)

    def _write_header(self, stream, comment_line):
        pass

    def _write_atoms(self, stream):
        formatter = self.formatter
        for atom in self.atoms:
            stream.write(formatter.format(atom))
            stream.write('\n')

CIF = CIFIO = CIFData


class CIFFormatter(StructureDataFormatter):
    """Class defining the structure file format for CIF data.

    Parameters
    ----------

    """
    def __init__(self, records=None):
        super().__init__()

        self.records = OrderedDict()

        # self._properties['RECORDS'] = records
        self.record_fmtstr = ('{0:6s}{1:5d} {2:4s}{3:1s}'
                              '{4:4s}{5:1s}{6:4d}{7:1s}   '
                              '{8:8.3f}{9:8.3f}{10:8.3f}'
                              '{11:6.2f}{12:6.2f}      '
                              '{13:4s}{14:2s}\n')

        self.fmtstr = "records={records!r}"

    def format(self, atom):
        """Return record :class:`~python:str`"""
        fmtstr = self.record_fmtstr
        line = fmtstr.format(**self.record_dict(atom))
        return line

    def record_dict(self, atom):
        """Return :class:`~python:dict` of record parameters."""
        # if hetfield != " ":
        #     record_type = "HETATM"
        # else:
        #     record_type = "ATOM  "
        # name = atom.get_fullname()
        # altloc = atom.get_altloc()
        # x, y, z = atom.get_coord()
        # bfactor = atom.get_bfactor()
        # occupancy = atom.get_occupancy()
        # return dict(record=record, )
        #     record_type, atom_number % 100000, name, altloc, resname, chain_id,
        #     resseq % 10000, icode, x, y, z, occupancy, bfactor, segid, element, charge)

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attr_dict = super().todict()
        attr_dict.update(dict(records=self.records))
        return attr_dict


CIFFormatSpec = CIFFormatter


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

CIFIOFormatter = CIFFormatter


class CIFError(StructureDataError):
    """Exception class for :class:`CIFData` I/O errors."""
    pass

CIFIOError = CIFError
