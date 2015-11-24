# -*- coding: utf-8 -*-
"""
===================================================================
Fullerene structure classes (:mod:`sknano.structures._fullerenes`)
===================================================================

.. currentmodule:: sknano.structures._fullerenes

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import Counter, OrderedDict
from operator import setitem
from pkg_resources import resource_filename

import numbers
import os

from sknano.core import BaseClass, listdir_dirnames, listdir_fnames
from sknano.core.crystallography import BaseStructure

__all__ = ['Fullerene', 'load_fullerene_data']


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
        files = dataset['files'] = sorted(listdir_fnames(datadir))
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

        [setitem(dataset, PG, OrderedDict([('files', dataset[PG]['files']),
                 ('isomer_numbers', dataset[PG]['isomer_numbers'])]))
         for PG in dataset['point_groups']]

        [[setitem(dataset[PG], Niso, dataset[PG]['files']
                  [dataset[PG]['isomer_numbers'].index(Niso)])
          for Niso in dataset[PG]['isomer_numbers']]
         for PG in dataset['point_groups']]

        if len(point_groups) == 1 and None in point_groups:
            dataset['point_groups'] = None
        else:
            dataset['point_group_counter'] = pg_counter

        # if len(isomer_numbers) == 1 and None in isomer_numbers:
        #     dataset['isomer_numbers'] = None

    return datadict


class Fullerene(BaseStructure, BaseClass):
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
    fullerene_data = load_fullerene_data()
    fullerenes = list(fullerene_data.keys())

    def __init__(self, *args, N=None, PG=None, Niso=None, **kwargs):

        super().__init__(*args, **kwargs)

        self._N = N
        self._PG = PG
        self._Niso = Niso

        self.fmtstr = "N={N!r}, PG={PG!r}, Niso={Niso!r}"

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

    def Natoms(self):
        """Alias for :attr:`~Fullerene.N`."""
        return self.N

    @property
    def PG(self):
        """Point Group symmetry in Schoenflies notation."""
        return self._PG

    @PG.setter
    def PG(self, value):
        if not isinstance(value, str):
            raise TypeError('Expected a string.')
        self._PG = value

    @property
    def Niso(self):
        """Isomer number."""
        return self._Niso

    @Niso.setter
    def Niso(self, value):
        if not isinstance(value, numbers.Integral):
            raise TypeError('Expected an integer.')
        self._Niso = int(value)

    def todict(self):
        """Return :class:`~python:dict` of `Fullerene` attributes."""
        return dict(N=self.N, PG=self.PG, Niso=self.Niso)
