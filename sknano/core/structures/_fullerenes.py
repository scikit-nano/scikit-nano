# -*- coding: utf-8 -*-
"""
===============================================================================
Fullerene structure classes (:mod:`sknano.core.structures._fullerenes`)
===============================================================================

.. currentmodule:: sknano.core.structures._fullerenes

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import Counter, OrderedDict
from operator import setitem
from pkg_resources import resource_filename

import numbers
import os

import numpy as np

from sknano.core import listdir_dirnames, listdir_fnames
from sknano.core.atoms import BasisAtom, BasisAtoms
from sknano.core.crystallography import Crystal3DLattice, UnitCell
from ._base import CrystalStructureBase

__all__ = ['Fullerene', 'Fullerenes', 'load_fullerene_data']


def load_fullerene_data():
    """Helper function to populate :class:`~python:dict` of \
        fullerene data files."""
    datadir = resource_filename('sknano', 'data/fullerenes')
    fullerenes = \
        listdir_dirnames(datadir, filterfunc=lambda name: name.startswith('C'))
    fullerenes = sorted(fullerenes, key=lambda s: int(s[1:]))

    datadict = OrderedDict()
    for fullerene in fullerenes:
        datapath = os.path.join('data', 'fullerenes', fullerene)
        datadir = resource_filename('sknano', datapath)
        dataset = datadict[fullerene] = {}
        files = dataset['files'] = \
            sorted(listdir_fnames(datadir, include_path=True))
        pg_counter = Counter()
        point_groups = dataset['point_groups'] = []
        # isomer_numbers = dataset['isomer_numbers'] = []

        for f in files:
            s = os.path.splitext(f)[0].split('-')
            if len(s) == 3:
                _, PG, Niso = s
                Niso = int(Niso)
            else:
                try:
                    Niso = int(s[-1])
                    PG = None
                except ValueError:
                    PG = s[-1]
                    Niso = None

            pg_counter[PG] += 1

            if PG not in point_groups:
                point_groups.append(PG)
                pg_dataset = dataset[PG] = {}
                pg_files = pg_dataset['files'] = []
                pg_isomer_numbers = pg_dataset['isomer_numbers'] = []

            pg_files.append(f)
            pg_isomer_numbers.append(Niso)
            pg_dataset[Niso] = f

        dataset['point_groups'] = \
            [t[0] for t in sorted(pg_counter.most_common(),
                                  key=lambda t: (-t[1], t[0]))]

        [setitem(dataset, PG_, OrderedDict([('files', dataset[PG_]['files']),
                 ('isomer_numbers', dataset[PG_]['isomer_numbers'])]))
         for PG_ in dataset['point_groups']]

        [[setitem(dataset[PG_], Niso_, dataset[PG_]['files']
                  [dataset[PG_]['isomer_numbers'].index(Niso_)])
          for Niso_ in dataset[PG_]['isomer_numbers']]
         for PG_ in dataset['point_groups']]

        if len(point_groups) == 1 and None in point_groups:
            dataset['point_groups'] = None
        else:
            dataset['point_group_counter'] = pg_counter

        # if len(isomer_numbers) == 1 and None in isomer_numbers:
        #     dataset['isomer_numbers'] = None
    return datadict


class Fullerene(CrystalStructureBase):
    """Fullerene structure class.

    The `fullerene data
    <http://www.nanotube.msu.edu/fullerene/fullerene-isomers.html>`_
    were downloaded from the
    `The Nanotube Site <http://www.nanotube.msu.edu>`_.

    Parameters
    ----------
    N : :class:`~python:int`
        The :math:`N` in :math:`C_N` where :math:`N` is the number
        of carbon atoms. Always an even integer.
    PG : :class:`~python:str`, optional
        Point Group symmetry of `Fullerene` isomer given in
        `Schoenflies notation
        <https://en.wikipedia.org/wiki/Schoenflies_notation>`_.
    Niso : :class:`~python:int`, optional
        Number specifying specific isomer with given `PG` symmetry.
        If not `None`, then `PG` must also be specified.

    Examples
    --------

    """
    _data = load_fullerene_data()

    def __init__(self, *args, N=None, PG=None, Niso=None, **kwargs):

        args = list(args)
        if N is None and len(args) == 1 and isinstance(args[0], int):
            N = args.pop()

        super().__init__(*args, **kwargs)

        if N is not None:
            self.N = N
            if PG is not None:
                self.PG = PG
                if Niso is not None:
                    self.Niso = Niso
        else:
            self._N = N
            self._PG = PG
            self._Niso = Niso

        self.generate_unit_cell()
        self.fmtstr = "N={N!r}, PG={PG!r}, Niso={Niso!r}"

    def __str__(self):
        strrep = "C{N!r}".format(N=self.N)
        return strrep

    @property
    def data(self):
        """Return :class:`~python:dict` of :class:`Fullerene` data."""
        try:
            return self._data[self.name]
        except KeyError:
            return None

    @property
    def datafile(self):
        """Return :class:`Fullerene` `xyz` structure data file."""
        try:
            return self.data[self.PG][self.Niso]
        except TypeError:
            return None

    @property
    def N(self):
        """Number of atoms."""
        return self._N

    @N.setter
    def N(self, value):
        if not (isinstance(value, numbers.Integral) or
                value > 0 or value % 2 == 0):
            raise TypeError('Expected an even, positive integer.')
        self._N = int(value)
        try:
            self._PG = self.point_groups[0]
        except TypeError:
            self._PG = None

        try:
            self._Niso = self.data[self.PG]['isomer_numbers'][0]
        except TypeError:
            self._Niso = None

    @property
    def Natoms(self):
        """Alias for :attr:`~Fullerene.N`."""
        return self.N

    @property
    def PG(self):
        """Point Group symmetry in Schoenflies notation."""
        return self._PG

    @PG.setter
    def PG(self, value):
        if ((self.point_groups is not None and
                value not in self.point_groups) or
                (self.point_groups is None and value is not None)):
            errmsg = 'No PG={} for {}'.format(value, self.name)
            raise ValueError(errmsg)
        self._PG = value

        try:
            self._Niso = self.data[self.PG]['isomer_numbers'][0]
        except TypeError:
            self._Niso = None

    @property
    def Niso(self):
        """Isomer number."""
        return self._Niso

    @Niso.setter
    def Niso(self, value):
        if value is not None and \
                value not in self.data[self.PG]['isomer_numbers']:
            errmsg = \
                'No isomer Niso={} for {}'.format(value, self.name)
            raise ValueError(errmsg)
        self._Niso = int(value) if isinstance(value, numbers.Integral) \
            else None

    @property
    def name(self):
        """Return string representation if :attr:`~Fullerene.N` is not `None` \
            else `None`."""
        if self.N is not None:
            return str(self)
        else:
            return None

    @property
    def Nisomers(self):
        """Return number of :attr:`~Fullerene.N` isomers."""
        try:
            return len(self.data['files'])
        except TypeError:
            return 0

    @property
    def point_groups(self):
        """List of point groups."""
        try:
            return self.data['point_groups']
        except TypeError:
            return None

    def generate_unit_cell(self):
        """Generate :class:`Fullerene` :class:`UnitCell`."""
        from sknano.io import XYZReader
        atoms = XYZReader(self.datafile).atoms
        atoms.center_centroid()
        bounding_box = atoms.bounding_box
        a = np.sqrt(2) * (max(bounding_box.abc) + 3.4)
        lattice = Crystal3DLattice.cubic(a)
        basis = BasisAtoms()
        for atom in atoms:
            xs, ys, zs = lattice.cartesian_to_fractional(atom.r)
            basis.append(BasisAtom(atom.element, lattice=lattice,
                                   xs=xs, ys=ys, zs=zs))

        self.unit_cell = UnitCell(lattice=lattice, basis=basis)

    def todict(self):
        """Return :class:`~python:dict` of :class:`Fullerene` parameters."""
        return dict(N=self.N, PG=self.PG, Niso=self.Niso)


class Fullerenes:
    data = load_fullerene_data()
    fullerenes = tuple(data.keys())
    N = tuple(int(fullerene[1:]) for fullerene in fullerenes)
