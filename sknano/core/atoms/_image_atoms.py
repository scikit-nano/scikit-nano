# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes with image id attributes (:mod:`sknano.core.atoms._image_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._image_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from functools import total_ordering
from operator import attrgetter

import numbers
import numpy as np

from sknano.core.math import Point
from ._atoms import Atom, Atoms

__all__ = ['ImageAtom', 'ImageAtoms']


@total_ordering
class ImageAtom(Atom):
    """An `Atom` sub-class with image count attributes.

    Parameters
    ----------
    ix, iy, iz : int, optional
        :math:`x, y, z` `ImageAtom` image count

    """
    def __init__(self, *args, ix=None, iy=None, iz=None, **kwargs):

        super().__init__(*args, **kwargs)
        self._i = Point([ix, iy, iz], dtype=int)
        self.fmtstr = super().fmtstr + ", ix={ix:d}, iy={iy:d}, iz={iz:d}"

    def __eq__(self, other):
        return self.i == other.i and super().__eq__(other)

    def __lt__(self, other):
        return (self.i < other.i and super().__le__(other)) or \
            (self.i <= other.i and super().__lt__(other))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['ix', 'iy', 'iz'])
        return attrs

    @property
    def ix(self):
        """:math:`i_x` image flag."""
        return self.i.x

    @ix.setter
    def ix(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.i.x = int(value)

    @property
    def iy(self):
        """:math:`i_y` image flag."""
        return self.i.y

    @iy.setter
    def iy(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.i.y = int(value)

    @property
    def iz(self):
        """:math:`i_z` image flag."""
        return self.i.z

    @iz.setter
    def iz(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.i.z = int(value)

    @property
    def i(self):
        """:math:`i_x, i_y, i_z` image flags

        Returns
        -------
        `Point`

        """
        return self._i

    @i.setter
    def i(self, value):
        """Set :math:`i_x, i_y, i_z` image flags.

        Parameters
        ----------
        value : array_like

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._i[:] = Point(value, nd=3, dtype=int)

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(ix=self.ix, iy=self.iy, iz=self.iz))
        return super_dict


class ImageAtoms(Atoms):
    """An `Atoms` sub-class for `ImageAtom`\ s.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.ImageAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `ImageAtoms`}, optional
        if not `None`, then a list of `ImageAtom` instance objects or an
        existing `ImageAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return ImageAtom

    def sort(self, key=attrgetter('i'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def images(self):
        """:class:`~numpy:numpy.ndarray` of `ImageAtom` images."""
        return np.asarray([atom.i for atom in self])

    @property
    def i(self):
        """Alias for :attr:`~ImageAtoms.images`."""
        return self.images

    @property
    def ix(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`i_x` values."""
        return self.i[:, 0]

    @property
    def iy(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`i_y` values."""
        return self.i[:, 1]

    @property
    def iz(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`i_z` values."""
        return self.i[:, 2]
