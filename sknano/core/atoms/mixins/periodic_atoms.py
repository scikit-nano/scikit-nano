# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin class for PBCs (:mod:`sknano.core.atoms.mixins.periodic_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins.periodic_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numpy as np

__all__ = ['PBCAtomsMixin']


class PBCAtomsMixin:
    """Mixin :class:`~sknano.core.atoms.mixins.Atoms` class for PBC.

    Parameters
    ----------
    xperiodic, yperiodic, zperiodic : :class:`~python:bool`, optional
        `True` for PBC along the given dimension. Default is `False`.

    """
    def __init__(self, *args, xperiodic=False, yperiodic=False,
                 zperiodic=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.xperiodic = xperiodic
        self.yperiodic = yperiodic
        self.zperiodic = zperiodic

    @property
    def xperiodic(self):
        """Return `True` if periodic along the `x` axis."""
        return self._xperiodic

    @xperiodic.setter
    def xperiodic(self, value):
        self._xperiodic = self.kwargs['xperiodic'] = bool(value)

    @property
    def yperiodic(self):
        """Return `True` if periodic along the `y` axis."""
        return self._yperiodic

    @yperiodic.setter
    def yperiodic(self, value):
        self._yperiodic = self.kwargs['yperiodic'] = bool(value)

    @property
    def zperiodic(self):
        """Return `True` if periodic along the `z` axis."""
        return self._zperiodic

    @zperiodic.setter
    def zperiodic(self, value):
        self._zperiodic = self.kwargs['zperiodic'] = bool(value)

    @property
    def pbc(self):
        """Return boolean array of periodic boundaries set along `x,y,z` axes.

        Returns
        -------
        :class:`~numpy:numpy.ndarray`
            Returns a :class:`~python:bool` :class:`~numpy:numpy.ndarray`
            of the :attr:`~PBCAtomsMixin.xperiodic`,
            :attr:`~PBCAtomsMixin.yperiodic`, :attr:`~PBCAtomsMixin.zperiodic`
            :class:`~python:bool` attributes.

        """
        return np.asarray([self.xperiodic, self.yperiodic, self.zperiodic],
                          dtype=bool)

    @pbc.setter
    def pbc(self, value):
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self.xperiodic, self.yperiodic, self.zperiodic = value

    def set_pbc(self, dims):
        """Set periodic boundaries along `dims`.

        Parameters
        ----------
        dims : :class:`~python:str`

        """
        for dim in 'xyz':
            if dim in dims:
                setattr(self, dim + 'periodic', True)

    def unset_pbc(self):
        """Turn of PBCs along all dimensions."""
        self.xperiodic = False
        self.yperiodic = False
        self.zperiodic = False
