# -*- coding: utf-8 -*-
"""
==============================================================================
Base molecule classes (:mod:`sknano.core.molecules._molecules`)
==============================================================================

.. currentmodule:: sknano.core.molecules._molecules

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

from sknano.core import BaseClass, UserList
from sknano.core.math import convert_condition_str, rotation_matrix

__all__ = ['Molecule', 'Molecules']


class Molecule(BaseClass):
    """Base class for abstract representation of a molecule.

    Per Wikipedia, "a molecule is an electrically neutral group of two
    or more atoms held together by chemical bonds."

    Parameters
    ----------
    atoms : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instances or an
        `Atoms` instance.

    """
    def __init__(self, atoms=None, **kwargs):
        super().__init__(**kwargs)
        self.atoms = atoms
        self.fmtstr = "atoms={atoms!r}"

    # def __eq__(self, other):
    #     """Test equality of two `Molecule` object instances."""
    #     return self is other or self.atoms == other.atoms

    # def __lt__(self, other):
    #     """Test if `self` is *less than* `other`."""
    #     return self.atoms < other.atoms

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.atoms == other.atoms and super().__eq__(other)

    def __le__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if not super().__le__(other) or self.atoms > other.atoms:
            return False
        return True

    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return (self.atoms < other.atoms and self.__le__(other))

    def __ge__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if not super().__ge__(other) or self.atoms < other.atoms:
            return False
        return True

    def __gt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.atoms > other.atoms and self.__ge__(other)

    def rezero(self, *args, **kwargs):
        assert not hasattr(super(), 'rezero')

    def rotate(self, **kwargs):
        assert not hasattr(super(), 'rotate')

    def translate(self, *args, **kwargs):
        assert not hasattr(super(), 'translate')

    def todict(self):
        """Return :class:`~python:dict` of `Molecule` constructor \
            parameters."""
        return dict(atoms=self.atoms)


class Molecules(UserList):
    """Base class for collection of `Molecule` objects.

    Parameters
    ----------
    molecules : {None, sequence, `Molecules`}, optional
        if not `None`, then a list of `Molecule` instance objects or an
        existing `Molecules` instance object.

    """
    def __init__(self, molecules=None, casttype=True, **kwargs):
        super().__init__(initlist=molecules)
        self.fmtstr = "{molecules!r}"

    @property
    def __item_class__(self):
        return Molecule

    def sort(self, key=attrgetter('id'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def Nmolecules(self):
        """Number of molecules in `Molecules`."""
        return len(self)

    @property
    def masses(self):
        """Return list of `Molecule` masses."""
        return np.asarray([molecule.mass for molecule in self])

    def filter(self, condition, invert=False):
        """Filter `Molecules` by `condition`.

        Parameters
        ----------
        condition : array_like, bool
        invert : bool, optional

        Returns
        -------
        filtered_molecules : `Molecules`

        """
        if isinstance(condition, str):
            condition = convert_condition_str(self, condition)
        if invert:
            condition = ~condition
        try:
            self.data = np.asarray(self)[condition].tolist()
        except AttributeError:
            self.data = np.asarray(self)[condition]

    def filtered(self, condition, invert=False):
        if isinstance(condition, str):
            condition = convert_condition_str(self, condition)
        if invert:
            condition = ~condition
        try:
            return self.__class__(atoms=np.asarray(self)[condition].tolist(),
                                  **self.kwargs)
        except AttributeError:
            return self.__class__(atoms=np.asarray(self)[condition],
                                  **self.kwargs)

    def get_molecules(self, asarray=False):
        """Return list of `Molecules`.

        Parameters
        ----------
        asarray : bool, optional

        Returns
        -------
        sequence or ndarray

        """
        if asarray:
            return np.asarray(self.data)
        else:
            return self.data

    def rezero(self, epsilon=1.0e-10):
        """Set really really small coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        [molecule.rezero(epsilon=epsilon) for molecule in self]

    def rotate(self, **kwargs):
        """Rotate `Molecule` position vectors.

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
        if kwargs.get('transform_matrix', None) is None:
            kwargs['transform_matrix'] = rotation_matrix(**kwargs)
        [molecule.rotate(**kwargs) for molecule in self]

    def translate(self, t, fix_anchor_points=True):
        """Translate `Molecule` position vectors by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_points : bool, optional

        """
        [molecule.translate(t, fix_anchor_point=fix_anchor_points)
         for molecule in self]

    def todict(self):
        return dict(molecules=self.data)
