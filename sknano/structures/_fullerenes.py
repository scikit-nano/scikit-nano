# -*- coding: utf-8 -*-
"""
===================================================================
Fullerene structure classes (:mod:`sknano.structures._fullerenes`)
===================================================================

.. currentmodule:: sknano.structures._fullerenes

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numbers

__all__ = ['Fullerene']


class Fullerene(object):
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

        super(Fullerene, self).__init__(**kwargs)

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
