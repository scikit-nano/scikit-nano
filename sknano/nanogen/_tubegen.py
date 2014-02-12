# -*- coding: utf-8 -*-
"""
=========================================================
TubeGen wrapper class (:mod:`sknano.nanogen._tubegen`)
=========================================================

.. currentmodule:: sknano.nanogen._tubegen

"""
from __future__ import division, print_function, absolute_import

import os
import subprocess
import sys

#from collections import OrderedDict

from math import ceil

#import numpy as np

from pkshared.tools.strfuncs import plural_word_check
from pkshared.tools.refdata import CCbond
from ._nanotube import Nanotube
from ..structure_io import XYZReader, XYZWriter

format_ext = {'gaussian': '.com',
              'gaussian-pbc': '.com',
              'wien': '.struct',
              'xyz': '.xyz',
              'pdb': '.pdb',
              'pdb-pbc': '.pdb',
              'pov': '.pov',
              'bgf': '.bgf',
              'cif': '.cif'}

shape_structure_dict = {'nanotube': 'hexagonal',
                        'graphene': 'planar',
                        'hexagonal-bundle': 'hexagonal',
                        'cubic-bundle': 'cubic',
                        'custom': None}

__all__ = ['TubeGen', 'TubeGenArgs', 'format_ext', 'shape_structure_dict']


class TubeGenArgs(object):
    """A class for populating tubegen arg attributes"""
    pass


class TubeGen(object):

    """Wrapper class around the TubeGen binary.

    TubeGen ([TG]_) is developed by Jeffrey T. Frey.

    Parameters
    ----------
    fmt : str, optional
        Structure data format. Must be one of:

            - xyz
            - gaussian
            - gaussian-pbc
            - wien
            - pdb
            - pdb-pbc
            - pov
            - bgf
            - cif

    units : {'angstrom', 'bohr'}, optional
        The unit of length to generate the spatial coordinates in.
    bond : float, optional
        Bond length between ``element1`` and ``element2``.
    element1, element2 : {int, str}, optional
        Atomic {number, symbol} of two atom basis.
    gutter : sequence, optional
        A 3-tuple of floats specifying the unit cell gutter in **Angstroms**.
    shape : {'hexagonal', 'cubic', 'planar'}, optional
        Selects the structure shape.
    chirality : sequence, optional
        A 2-tuple of ints specifying the chiral indices defining the
        chirality :math:`\\mathbf{C}_{h} = (n, m)`
    cell_count : sequence, optional
        A 3-tuple of ints specifying the unit cell count along the
        :math:`(x, y, z)` axes.
    relax_tube : {'yes', 'no'}, optional
        Selects whether or not to apply energy minimization of the structure.
    tube_length : float, optional
        length of the nanotube in **nanometers**.

    Examples
    --------

    >>> from sknano.nanogen import TubeGen

    Generate a single unit cell of a :math:`\\mathbf{C}_{h} = (10, 5)`
    :abbr:`SWCNT (single-wall carbon nanotube)` in *xyz* format:

    .. code-block:: python

       >>> swcnt1005 = TubeGen(chirality=(10, 5))
       >>> swcnt1005.generate()
       TubeGen - Carbon Nanotube Stucture Generator
       Version 3.4 [03:16:28 Oct 29 2013]
       Copyright (c) 2000-2013, J. T. Frey

       > set element1 C
       > set element2 C
       > set gutter 1.6735,1.6735,0
       > set format xyz
       > set units angstrom
       > set cell_count 1,1,1
       > set chirality 10,5
       --- Construction of Graphitic Atomic-Basis Vectors ---------------------
       Nearest neighbor bonding distance as:  1.421
           a1 = < 2.1315 , 1.23062 >
           a2 = < 2.1315 , -1.23062 >
       ------------------------------------------------------------------------

       --- Construction of Chiral/Tubule Translation Vectors ------------------
       n = ( 10) and m = (  5):
       n'= (  4) and m'= (  5):
       Chiral vector Ch constructed as 10(a1) + 5(a2):
           Ch = < 31.9725 , 6.15311 >, |Ch| = 32.5592
       Chiral angle is 19.1066 degrees
       Tubule translation vector T constructed as 4(a1) - 5(a2):
           T = < -2.1315 , 11.0756 >, |T| = 11.2788
       Tubule radius: 5.18196
       Tubule height: 11.2788
       Reciprocal-space vectors:
           T_R = < -0.0167555 , 0.087064 >
           Ch_R = < 0.0301598 , 0.00580427 >
       ------------------------------------------------------------------------

       --- Relaxing tubule to appropriate bond lengths ------------------------
                                                   delta-radius           1e-15
                                                   scaling-factors        1e-15
                                                   error-function         1e-15
       ========================================================================
                                                                     Iterations
       Iter     delta-h      delta-r      Gammas                     a1  a2  a3
       ========================================================================
       0                               1.000000  1.000000  1.000000  45  39  30
       1     3.5044e-03   1.0648e-02   1.002936  1.000576  1.000036  26  21   8
       2    -1.5682e-05  -4.3187e-05   1.002925  1.000573  1.000036   6   0   0
       3     9.0760e-09   1.1291e-07   1.002925  1.000573  1.000036   0   0   0
       4     0.0000e+00   0.0000e+00   1.002925  1.000573  1.000036
       ========================================================================
       Convergence reached in 4 cycles
       New graphitic basis:
           a1 = < 2.13606 , 1.23133 >
           a2 = < 2.13568 , -1.23067 >
           cc-bond = 1.42516
       New chiral/tubule translation vectors:
           Ch = < 32.039 , 6.15994 >, |Ch| = 32.6258
           T = < -2.13415 , 11.0786 >, |T| = 11.2823
       Tubule radius: 5.19256     [total delta-r of 0.0106047]
       Tubule height: 11.2823     [total delta-h of 0.00348869]
       Angle between Ch and T:  97.6112 degrees
       ------------------------------------------------------------------------

       > set shape hexagonal
       > set relax_tube yes
       > set bond 1.421
       --- Construction of Graphitic Atomic-Basis Vectors ---------------------
       Nearest neighbor bonding distance as:  1.421
           a1 = < 2.1315 , 1.23062 >
           a2 = < 2.1315 , -1.23062 >
       ------------------------------------------------------------------------

       --- Construction of Chiral/Tubule Translation Vectors ------------------
       n = ( 10) and m = (  5):
       n'= (  4) and m'= (  5):
       Chiral vector Ch constructed as 10(a1) + 5(a2):
           Ch = < 31.9725 , 6.15311 >, |Ch| = 32.5592
       Chiral angle is 19.1066 degrees
       Tubule translation vector T constructed as 4(a1) - 5(a2):
           T = < -2.1315 , 11.0756 >, |T| = 11.2788
       Tubule radius: 5.18196
       Tubule height: 11.2788
       Reciprocal-space vectors:
           T_R = < -0.0167555 , 0.087064 >
           Ch_R = < 0.0301598 , 0.00580427 >
       ------------------------------------------------------------------------

       --- Relaxing tubule to appropriate bond lengths ------------------------
                                                   delta-radius           1e-15
                                                   scaling-factors        1e-15
                                                   error-function         1e-15
       ========================================================================
                                                                     Iterations
       Iter    delta-h      delta-r      Gammas                      a1  a2  a3
       ========================================================================
       0                               1.000000  1.000000  1.000000  45  39  30
       1     3.5044e-03   1.0648e-02   1.002936  1.000576  1.000036  26  21   8
       2    -1.5682e-05  -4.3187e-05   1.002925  1.000573  1.000036   6   0   0
       3     9.0760e-09   1.1291e-07   1.002925  1.000573  1.000036   0   0   0
       4     0.0000e+00   0.0000e+00   1.002925  1.000573  1.000036
       ========================================================================
       Convergence reached in 4 cycles
       New graphitic basis:
           a1 = < 2.13606 , 1.23133 >
           a2 = < 2.13568 , -1.23067 >
           cc-bond = 1.42516
       New chiral/tubule translation vectors:
           Ch = < 32.039 , 6.15994 >, |Ch| = 32.6258
           T = < -2.13415 , 11.0786 >, |T| = 11.2823
       Tubule radius: 5.19256     [total delta-r of 0.0106047]
       Tubule height: 11.2823     [total delta-h of 0.00348869]
       Angle between Ch and T:  97.6112 degrees
       ------------------------------------------------------------------------

       > generate
       Producing rolled, hexagonal nanotube lattice.

       Lattice consists of 70 hexagonal sub-cells.
       Cell generation complete.  140 basis points defined.
       > save 1005r_1cell.xyz
       > tubegen finished with no errors


    References
    ----------
    .. [TG] TubeGen 3.4 (web-interface,
            http://turin.nss.udel.edu/research/tubegenonline.html),
            J. T. Frey and D. J. Doren, University of Delaware,
            Newark DE, 2011

    """

    def __init__(self, fmt='xyz', units='angstrom', bond=CCbond,
                 element1='C', element2='C', gutter=(1.6735, 1.6735, 0),
                 shape='hexagonal', chirality=(10, 10), cell_count=(1, 1, 1),
                 relax_tube='yes', tube_length=None):

        self.n = int(chirality[0])
        self.m = int(chirality[1])

        self.T = Nanotube.compute_T(n=self.n, m=self.m)

        nx, ny, nz = cell_count

        self._tube_length = None
        if tube_length is not None and tube_length != 'None' \
           and shape == 'hexagonal' and \
           isinstance(tube_length, (int, float, str)):
            self._tube_length = 10 * float(tube_length)
            nz = int(ceil(self._tube_length / self.T))
        elif nz is not None and nz != 'None' and \
                isinstance(nz, (int, float, str)):
            nz = int(nz)
            self._tube_length = nz * self.T

        self._fmt = fmt
        self._units = units
        self._bond = bond
        self._e1 = element1
        self._e2 = element2
        self._gutter = ','.join([str(x) for x in gutter])
        self._shape = shape
        self._chirality = ','.join([str(x) for x in chirality])
        self._relax_tube = relax_tube
        self._cell_count = ','.join([str(x) for x in
                                     (nx, ny, nz)])
        self._genfile_name = ''
        self._kwargs = {'format': self._fmt,
                        'units': self._units,
                        'bond': self._bond,
                        'element1': self._e1,
                        'element2': self._e2,
                        'gutter': self._gutter,
                        'shape': self._shape,
                        'chirality': self._chirality,
                        'cell_count': self._cell_count,
                        'relax_tube': self._relax_tube}
        self._output = ''

        self._dt = Nanotube.compute_dt(n=self.n, m=self.m)

    @property
    def fmt(self):
        """Structure data format.

        Returns
        -------
        fmt : str
            Structure data format.

        """
        return self._fmt

    @property
    def cell_count(self):
        """Cell count.

        Returns
        -------
        (int, int, int)
            tuple of `nx`, `ny`, `nz`

        """
        return tuple([int(c) for c in self._cell_count.split(',')])

    @property
    def chirality(self):
        """Chirality.

        Returns
        -------
        (int, int)
            chiral indices (`n`, `m`)

        """
        return tuple([int(c) for c in self._chirality.split(',')])

    @property
    def dt(self):
        """Tube diameter.

        Returns
        dt : float
            nanotube diameter in Angstroms

        """
        return self._dt

    @property
    def gutter(self):
        """Gutter.

        Returns
        -------
        (float, float, float)
            tuple of gutter values

        """
        return tuple([float(c) for c in self._gutter.split(',')])

    @property
    def output(self):
        """Output file.

        Returns
        -------
        output : str
            output file name

        """
        return self._output

    @property
    def shape(self):
        """Structure shape.

        Returns
        -------
        str : {'hexagonal', 'cubic', 'planar'}
            structure shape

        """
        return self._shape

    def fix_length(self):
        """Crop nanotube ends for correct length."""
        if self._tube_length is not None:
            if self._fmt == 'xyz':
                xyzreader = XYZReader(self._output)
                atoms = xyzreader.atoms
                atoms.clip_bounds(abs_limit=self._tube_length / 2 + 0.1,
                                  r_indices=[2])
                XYZWriter.write(fname=self._output, atoms=atoms,
                                comment_line=xyzreader.comment_line)

    def cleanup(self):
        """Cleanup tubegen input file."""

        if os.path.isfile(self._genfile_name):
            try:
                os.remove(self._genfile_name)
            except OSError as e:
                print(e)
            else:
                self._genfile_name = ''

    def generate(self):
        """Generate structure."""
        self._genfile_name = self._name_genfile()
        try:
            self._write_genfile()
        except OSError as e:
            print(e)
            sys.exit(1)

        self._tubegen()
        self._output = \
            self._genfile_name.split('.')[0] + format_ext[self._fmt]

    def _name_genfile(self):

        genfile = None
        cell_count = self._cell_count
        chirality = self._chirality
        shape = self._shape
        xcells, ycells, zcells = tuple([int(c) for c in cell_count.split(',')])
        if cell_count == '1,1,1':
            genfile = ''.join([x.zfill(2) for x in chirality.split(',')]) \
                + ('f' if shape == 'planar' else 'r') + '_1cell' + '.gen'
        elif shape != 'planar' \
                and xcells == 1 and ycells == 1 and zcells != 1:
            genfile = ''.join([x.zfill(2) for x in chirality.split(',')]) \
                + 'r_' + str(zcells) + 'cells' + '.gen'
        elif shape == 'planar' \
                and (xcells != 1 or ycells != 1) and zcells == 1:
            genfile = ''.join([x.zfill(2) for x in chirality.split(',')]) \
                + 'f_' + 'x'.join([str(n) + plural_word_check('cell', n)
                                   for n in (xcells, ycells)]) + '.gen'
        else:
            genfile = ''.join([x.zfill(2) for x in chirality.split(',')]) \
                + ('f' if shape == 'planar' else 'r') + '_' + \
                'x'.join([n + plural_word_check('cell', n)
                          for n in cell_count.split(',')]) + '.gen'

        return genfile

    def _write_genfile(self):
        f = open(self._genfile_name, 'w')
        for k, v in self._kwargs.iteritems():
            f.write('set {!s} {!s}\n'.format(k, v))
        f.write('generate\n')
        f.write('save {!s}{!s}\n'.format(self._genfile_name.split('.')[0],
                                         format_ext[self._fmt]))
        f.close()

        return None

    def _tubegen(self):
        retcode = subprocess.call(["tubegen", self._genfile_name])
        if retcode != 0:
            print('tubegen failed')
        else:
            print('tubegen finished with no errors')

        return None
