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

from pkg_resources import resource_filename

import numbers
import os

from sknano.core import BaseClass, listdir_dirnames, listdir_fnames
from sknano.core.crystallography import BaseStructure

__all__ = ['Fullerene', 'Fullerenes', 'load_fullerene_data']


def load_fullerene_data():
    """Helper function to populate dict of fullerene data files."""
    datadir = resource_filename('sknano', 'data/fullerenes')
    fullerenes = \
        listdir_dirnames(datadir, filterfunc=lambda name: name.startswith('C'))
    fullerene_data = {}
    for fullerene in fullerenes:
        datapath = os.path.join('data', 'fullerenes', fullerene)
        datadir = resource_filename('sknano', datapath)
        fullerene_data[fullerene] = listdir_fnames(datadir)
    return fullerene_data


class Fullerene(BaseStructure, BaseClass):
    """Fullerene structure class.

    The `fullerene data
    <http://www.nanotube.msu.edu/fullerene/fullerene-isomers.html>`_
    were downloaded from the
    `The Nanotube Site <http://www.nanotube.msu.edu>`_.

    Parameters
    ----------
    N : int
        The :math:`N` in :math:`C_N` where :math:`N` is the number
        fullerene atoms. Always an even integer.
    PG : str, optional
        Symmetry

    Examples
    --------

    """
    def __init__(self, N=60, PG=None, Ni=None, **kwargs):

        super().__init__(**kwargs)

        # Check that N is a valid `N` and that `N` is
        self.N = N
        self.PG = PG
        self.Ni = Ni

    @property
    def N(self):
        return self._N

    @N.setter
    def N(self, value):
        if not (isinstance(value, numbers.Integral) or
                value > 0 or value % 2 == 0):
            raise TypeError('Expected an even, positive integer.')
        self._N = int(value)

    @N.deleter
    def N(self):
        del self._N

    def todict(self):
        return dict(N=self.N)


class Fullerenes(BaseStructure, BaseClass):
    def __init__(self):
        super().__init__()
        self.fullerene_data = load_fullerene_data()

        self.fullerenes = list(self.fullerene_data.keys())
        self.fullerene_files = list(self.fullerene_data.values())

    def todict(self):
        return dict()
