# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class for crystal structure basis (:mod:`sknano.core.atoms._basis_atom`)
===============================================================================

An `Atom` sub-class for a crystal structure basis atom.

.. currentmodule:: sknano.core.atoms._basis_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from functools import total_ordering
import numbers
import numpy as np

from sknano.core.math import Vector
from ._atom import Atom

__all__ = ['BasisAtom']


@total_ordering
class BasisAtom(Atom):
    """An abstract object representation of a crystal structure basis atom.

    Parameters
    ----------
    basis_point : `Point`

    """

    def __init__(self, frac_coords=None, **kwargs):

        super().__init__(**kwargs)

        self.frac_coords = Vector(frac_coords)

        self.fmtstr = "Atom(frac_coords={frac_coords!s}, " + \
            "element={element!s}, mass={mass!s})"

    def __str__(self):
        """Return a nice string representation of `BasisAtom`."""
        try:
            return self.fmtstr.format(**self.todict())
        except KeyError:
            return repr(self)

    def __repr__(self):
        """Return canonical string representation of `BasisAtom`."""
        return "Atom(coords={!r}, element={!r}, mass={!r})".format(
            self.frac_coords, self.element, self.mass)

    def __eq__(self, other):
        return super().__eq__(other)

    def __lt__(self, other):
        """Test if `self` is *less than* `other`."""
        if not isinstance(type(self)):
            return NotImplemented
        if self.frac_coords.length < other.frac_coords.length:
            return True
        else:
            return super().__lt__(other)

    def __dir__(self):
        attrs = super().__dir__()
        attrs.append('frac_coords')
        return attrs

    @property
    def x(self):
        """:math:`x`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`x`-coordinate in units of **Angstroms**.

        """
        return self._r.x

    @x.setter
    def x(self, value):
        """Set `Atom` :math:`x`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`x`-coordinate in units of **Angstroms**.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._r.x = value

    @property
    def y(self):
        """:math:`y`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`y`-coordinate in units of **Angstroms**.

        """
        return self._r.y

    @y.setter
    def y(self, value):
        """Set `Atom` :math:`y`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`y`-coordinate in units of **Angstroms**.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._r.y = value

    @property
    def z(self):
        """:math:`z`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`z`-coordinate in units of **Angstroms**.

        """
        return self._r.z

    @z.setter
    def z(self, value):
        """Set `Atom` :math:`z`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`z`-coordinate in units of **Angstroms**.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._r.z = value

    @property
    def r(self):
        """:math:`x, y, z` components of `Atom` position vector.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x, y, z`] coordinates of `Atom`.

        """
        return self._r

    @r.setter
    def r(self, value):
        """Set :math:`x, y, z` components of `Atom` position vector.

        Parameters
        ----------
        value : array_like
            :math:`x, y, z` coordinates of `Atom` position vector relative to
            the origin.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._r = Vector(value)

    @property
    def frac_coords(self):
        return self.r

    @frac_coords.setter
    def frac_coords(self, value):
        self.r = value

    def rezero(self, epsilon=1.0e-10):
        """Re-zero position vector components.

        Set position vector components with absolute value less than
        `epsilon` to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        self._r.rezero(epsilon=epsilon)

    def rotate(self, angle=None, axis=None, anchor_point=None,
               rot_point=None, from_vector=None, to_vector=None, degrees=False,
               transform_matrix=None, verbose=False, **kwargs):
        """Rotate `Atom` position vector.

        Parameters
        ----------
        angle : float
        axis : :class:`~sknano.core.math.Vector`, optional
        anchor_point : :class:`~sknano.core.math.Point`, optional
        rot_point : :class:`~sknano.core.math.Point`, optional
        from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
        degrees : bool, optional
        transform_matrix : :class:`~numpy:numpy.ndarray`

        """
        self.r.rotate(angle=angle, axis=axis,
                      anchor_point=anchor_point,
                      rot_point=rot_point, from_vector=from_vector,
                      to_vector=to_vector, transform_matrix=transform_matrix,
                      degrees=degrees, verbose=verbose, **kwargs)

    def translate(self, t, fix_anchor_point=True):
        """Translate `Atom` position vector by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_point : bool, optional

        """
        self.r.translate(t, fix_anchor_point=fix_anchor_point)

    def todict(self):
        return dict(frac_coords=self.frac_coords, element=self.element,
                    mass=self.mass)
