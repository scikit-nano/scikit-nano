# -*- coding: utf-8 -*-
"""
==============================================================================
Base Atom classes (:mod:`sknano.core.atoms._atoms`)
==============================================================================

.. currentmodule:: sknano.core.atoms._atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from functools import total_ordering
from operator import attrgetter
import numbers
import numpy as np

from sknano.core import BaseClass, UserList
from sknano.core.math import convert_condition_str, rotation_matrix
from sknano.core.refdata import atomic_masses, atomic_mass_symbol_map, \
    atomic_numbers, atomic_number_symbol_map, element_symbols, element_names

__all__ = ['Atom', 'Atoms']


@total_ordering
class Atom(BaseClass):
    """Base class for abstract representation of structure atom.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number :math:`\\boldsymbol{Z}`.

    """
    # _fields = ['element']

    def __init__(self, *args, element=None, mass=None, Z=None, **kwargs):
        args = list(args)

        if 'm' in kwargs and mass is None:
            mass = kwargs['m']
            del kwargs['m']

        if element is None:

            if mass is not None or Z is not None:
                if Z is not None:
                    args.append(Z)
                else:
                    args.append(mass)

            if len(args) > 0:
                element = args.pop()

        super().__init__(*args, **kwargs)
        self.mass = mass
        self.element = element
        self.fmtstr = "{element!r}, Z={Z!r}, mass={mass!r}"

    def __eq__(self, other):
        """Test equality of two `Atom` object instances."""
        return self is other or (self.element == other.element and
                                 self.Z == other.Z and
                                 np.allclose(self.mass, other.mass))

    def __lt__(self, other):
        """Test if `self` is *less than* `other`."""
        if self.element == other.element == 'X' and \
                self.Z == other.Z == 0:
            return self.mass < other.mass
        else:
            return self.Z < other.Z

    def __dir__(self):
        return ['element', 'Z', 'mass']

    @property
    def Z(self):
        """Atomic number :math:`Z`.

        Returns
        -------
        int
            Atomic number :math:`Z`.
        """
        return self._Z

    @Z.setter
    def Z(self, value):
        """Set atomic number :math:`Z`.

        Parameters
        ----------
        value : int
            Atomic number :math:`Z`.

        """
        if not (isinstance(value, numbers.Real) and int(value) > 0):
            raise ValueError('Expected a real, positive integer.')
        try:
            Z = int(value)
            idx = Z - 1
            symbol = element_symbols[idx]
            mass = atomic_masses[symbol]
        except KeyError:
            print('unrecognized element number: {}'.format(value))
        else:
            self._Z = atomic_numbers[symbol]
            self._mass = mass
            self._symbol = symbol

    @property
    def element(self):
        """Element symbol.

        Returns
        -------
        str
            Symbol of chemical element.
        """
        return self._symbol

    @element.setter
    def element(self, value):
        """Set element symbol."""
        symbol = None
        if isinstance(value, numbers.Integral):
            try:
                Z = int(value)
                idx = Z - 1
                symbol = element_symbols[idx]
            except IndexError:
                print('unrecognized element number: {}'.format(value))

        if symbol is None:
            if isinstance(value, str):
                if value in element_symbols:
                    symbol = value
                elif value.capitalize() in element_names:
                    symbol = element_symbols[element_names.index(value)]
            elif isinstance(value, numbers.Number):
                if value in atomic_mass_symbol_map:
                    symbol = atomic_mass_symbol_map[value]
                elif int(value / 2) in atomic_number_symbol_map:
                    symbol = atomic_number_symbol_map[int(value / 2)]

        if symbol is None:
            symbol = 'X'

        try:
            self._Z = atomic_numbers[symbol]
            self._mass = atomic_masses[symbol]
        except KeyError:
            self._Z = 0
            if self.mass is None:
                self._mass = 0

        self._symbol = symbol

    @property
    def symbol(self):
        """Element symbol.

        Returns
        -------
        str
            Element symbol.
        """
        return self._symbol

    @property
    def mass(self):
        """Atomic mass :math:`m_a` in atomic mass units.

        Returns
        -------
        float
            Atomic mass :math:`m_a` in atomic mass units.
        """
        return self._mass

    @mass.setter
    def mass(self, value):
        self._mass = value

    @property
    def m(self):
        return self.mass

    @m.setter
    def m(self, value):
        self.mass = value

    def rezero(self, *args, **kwargs):
        assert not hasattr(super(), 'rezero')

    def rotate(self, **kwargs):
        assert not hasattr(super(), 'rotate')

    def translate(self, *args, **kwargs):
        assert not hasattr(super(), 'translate')

    def todict(self):
        """Return `dict` of `Atom` constructor parameters."""
        return dict(element=self.element, mass=self.mass, Z=self.Z)


class Atoms(UserList):
    """Base class for collection of `Atom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instance objects or an
        existing `Atoms` instance object.

    """
    def __init__(self, atoms=None, casttype=True, **kwargs):
        verbose = kwargs.get('verbose', False)
        if atoms is not None and (isinstance(atoms, str) or
                                  isinstance(atoms, Atom)):
            atoms = [atoms]
        if casttype and not isinstance(atoms, type(self)) and \
                isinstance(atoms, list):
            atoms = atoms[:]
            for i, atom in enumerate(atoms):
                try:
                    # atoms[i] = self.__atom_class__(**atom.todict())
                    atomdict = atom.todict()
                    if verbose and i in list(range(len(atoms), 100)):
                        print(type(atom))
                        print(atomdict)
                        filtered_atomdict = \
                            {k: atomdict[k] for k in set(dir(atom)) &
                             set(dir(self.__atom_class__()))}
                        print('filtered_atomdict: {}'.format(
                              filtered_atomdict))

                    atoms[i] = self.__atom_class__(
                        **{k: atomdict[k] for k in set(dir(atom)) &
                           set(dir(self.__atom_class__()))})
                except AttributeError:
                    atoms[i] = self.__atom_class__(atom)
        super().__init__(initlist=atoms, **kwargs)
        self.fmtstr = "{atoms!r}"

    @property
    def __atom_class__(self):
        return Atom

    def __repr__(self):
        """Return canonical string representation of `Atoms`."""
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    def sort(self, key=attrgetter('element', 'Z', 'mass'), reverse=False):
        super().sort(key=key, reverse=reverse)

    def __getitem__(self, index):
        data = super().__getitem__(index)
        if isinstance(data, list):
            return self.__class__(data, **self.kwargs)
        return data

    def __setitem__(self, index, item):
        if not isinstance(item, (self.__class__, self.__atom_class__)):
            if isinstance(index, slice):
                item = self.__class__(item)
            else:
                try:
                    atomdict = super().__getitem__(index).todict()
                    if isinstance(item, str):
                        atomdict['element'] = item
                    elif isinstance(item, int):
                        atomdict['Z'] = item
                    elif isinstance(item, float):
                        atomdict['mass'] = item
                    else:
                        raise ValueError
                    item = self.__atom_class__(**atomdict)
                except (IndexError, AttributeError, ValueError):
                    item = self.__atom_class__(item)
        super().__setitem__(index, item)

    def append(self, atom):
        if not isinstance(atom, self.__atom_class__):
            if isinstance(atom, dict):
                atom = self.__atom_class__(
                    **{k: atom[k] for k in set(atom) &
                       set(dir(self.__atom_class__()))})
            elif isinstance(atom, (str, int, float)):
                atom = self.__atom_class__(atom)
        super().append(atom)

    def __atoms_not_in_self(self, other):
        return [atom for atom in other if atom not in self]

    def __add__(self, other):
        try:
            addlist = self.__atoms_not_in_self(other)
        except TypeError:
            addlist = self.__atoms_not_in_self(self.__class__(other))
        return self.__class__(self.data + addlist, **self.kwargs)

    def __radd__(self, other):
        try:
            addlist = self.__atoms_not_in_self(other)
        except TypeError:
            addlist = self.__atoms_not_in_self(self.__class__(other))
        return self.__class__(addlist + self.data, **self.kwargs)

    def __iadd__(self, other):
        try:
            self.data += self.__atoms_not_in_self(other)
        except TypeError:
            self.data += self.__atoms_not_in_self(self.__class__(other))
        return self

    # def __atoms_not_in_other(self, other):
    #     return [atom for atom in self if atom not in other]

    # def __sub__(self, other):
    #     if isinstance(other, self.__class__):
    #         return self.__class__(self.__atoms_not_in_other(other),
    #                               **self.kwargs)
    #     return NotImplemented

    # def __rsub__(self, other):
    #     return NotImplemented

    # def __isub__(self, other):
    #     if isinstance(other, self.__class__):
    #         self.data -= self.__atoms_not_in_other(other)
    #     return self

    @property
    def fmtstr(self):
        return self._fmtstr

    @fmtstr.setter
    def fmtstr(self, value):
        self._fmtstr = value

    @property
    def Natoms(self):
        """Number of atoms in `Atoms`."""
        return len(self)

    @property
    def M(self):
        """Total mass of `Atoms`."""
        # return math.fsum(self.masses)
        return self.masses.sum()

    @property
    def elements(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`Atom.element`\ s \
            in `Atoms`."""
        return np.asarray([atom.element for atom in self])

    @property
    def masses(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`Atom.mass`\ s \
            in `Atoms`."""
        return np.asarray([atom.mass for atom in self])

    @property
    def symbols(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`Atom.symbol`\ s \
            in `Atoms`."""
        return np.asarray([atom.symbol for atom in self])

    def filter(self, condition, invert=False):
        """Filter `Atoms` by `condition`.

        .. versionchanged:: 0.3.11

           Filters the list of `Atoms` **in-place**. Use
           :meth:`~Atoms.filtered` to generate a new filtered list
           of `Atoms`.

        Parameters
        ----------
        condition : :class:`~python:str` or boolean array
            Boolean index array having same shape as the initial dimensions
            of the list of `Atoms` being indexed.
        invert : bool, optional
            If `True`, the boolean array `condition` is inverted element-wise.

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
        """Return new list of `Atoms` filtered by `condition`.

        .. versionadded:: 0.3.11

        Parameters
        ----------
        condition : array_like, bool
            Boolean index array having same shape as the initial dimensions
            of the list of `Atoms` being indexed.
        invert : bool, optional
            If `True`, the boolean array `condition` is inverted element-wise.

        Returns
        -------
        filtered_atoms : `Atoms`
            If `invert` is `False`, return the elements where `condition`
            is `True`.

            If `invert` is `True`, return the elements where `~condition`
            (i.e., numpy.invert(condition)) is `True`.

        Examples
        --------
        An example using the structure data of a 10 nm `(10, 0)`
        `SWCNT`:

        >>> from sknano.generators import SWNTGenerator
        >>> swnt = SWNTGenerator(10, 0, Lz=10, fix_Lz=True).atoms
        >>> # select 'left', 'middle', 'right' atoms
        >>> latoms = swnt.filtered(swnt.z <= 25)
        >>> matoms = swnt.filtered((swnt.z < 75) & (swnt.z > 25))
        >>> ratoms = swnt.filtered(swnt.z >= 75)
        >>> from pprint import pprint
        >>> pprint([getattr(atoms, 'bounds') for atoms in
        ...         (latoms, matoms, ratoms)])
        [Cuboid(pmin=Point([-3.914435, -3.914435, 0.0]), \
                pmax=Point([3.914435, 3.914435, 24.85])),
         Cuboid(pmin=Point([-3.914435, -3.914435, 25.56]), \
                pmax=Point([3.914435, 3.914435, 74.55])),
         Cuboid(pmin=Point([-3.914435, -3.914435, 75.97]), \
                pmax=Point([3.914435, 3.914435, 100.11]))]
        >>> latoms.Natoms + matoms.Natoms + ratoms.Natoms == swnt.Natoms
        True

        """
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

    def get_atoms(self, asarray=False):
        """Return list of `Atoms`.

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

    def getatomattr(self, attr):
        """Get :class:`~numpy:numpy.ndarray` of atom attributes `attr`.

        Parameters
        ----------
        attr : str
            Name of attribute to pass to `getattr` from each `XAtom` in
            `XAtoms`.

        Returns
        -------
        :class:`~numpy:numpy.ndarray`

        """
        try:
            return np.asarray([getattr(atom, attr) for atom in self])
        except AttributeError:
            return None

    def mapatomattr(self, from_attr=None, to_attr=None, attrmap=None):
        """Set/update atom attribute from another atom attribute with dict.

        .. versionchanged:: 0.3.11

           Made all arguments required keyword arguments and reversed the
           order of the former positional arguments `from_attr` and
           `to_attr` to be more natural and consistent with the key, value
           pairs in the `attrmap` dictionary.

        Parameters
        ----------
        from_attr, to_attr : :class:`python:str`
        attrmap : :class:`python:dict`

        Examples
        --------
        Suppose you have an `XAtoms` instance named ``atoms`` that has
        `XAtom` instances of two atom types `1` and `2` and we want to set
        all `XAtom`\ s with `type=1` to Nitrogen and all `XAtom`\ s with
        `type=2` to Argon. In other words, we want to map the
        `XAtom.type` attribute to the `XAtom.element` attribute.

        We'd call this method like so::

        >>> atoms.mapatomattr('type', 'element', {1: 'N', 2: 'Ar'})

        """
        try:
            [setattr(atom, to_attr, attrmap[getattr(atom, from_attr)])
             for atom in self if getattr(atom, from_attr) is not None]
        except (KeyError, TypeError) as e:
            print(e)

    def rezero(self, epsilon=1.0e-10):
        """Set values with absolute value less than `epsilon` to zero.

        Calls the `rezero` method on each `atom` in `self`.

        Parameters
        ----------
        epsilon : :class:`~python:float`
            values with absolute value less than `epsilon` are set to zero.

        """
        [atom.rezero(epsilon=epsilon) for atom in self]

    def rotate(self, **kwargs):
        """Rotate `Atom` vectors.

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
        [atom.rotate(**kwargs) for atom in self]

    def translate(self, t, fix_anchor_points=True):
        """Translate `Atom` vectors by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_points : bool, optional

        """
        [atom.translate(t, fix_anchor_point=fix_anchor_points)
         for atom in self]

    def todict(self):
        return dict(atoms=self.data)
