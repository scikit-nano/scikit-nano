# -*- coding: utf-8 -*-
"""
==============================================================================
Base Atom classes (:mod:`sknano.core.atoms.atoms`)
==============================================================================

.. currentmodule:: sknano.core.atoms.atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections.abc import Iterable
from operator import attrgetter
import copy
import numbers
import re
# import warnings

import numpy as np

from sknano.core import BaseClass, UserList, TabulateMixin, dedupe
from sknano.core.math import convert_condition_str
from sknano.core.refdata import atomic_masses, atomic_mass_symbol_map, \
    atomic_numbers, atomic_number_symbol_map, element_symbols, element_names
from .selections import AtomsSelectionMixin

__all__ = ['Atom', 'Atoms', 'update_atoms']


class Atom(BaseClass):
    """Base class for abstract representation of structure atom.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number :math:`\\boldsymbol{Z}`.

    """
    _fields = ['element', 'Z', 'mass']

    def __init__(self, *args, element=None, mass=None, Z=None, parent=None,
                 **kwargs):
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
        self.parent = parent
        self.fmtstr = "{element!r}, Z={Z!r}, mass={mass!r}"

    def _is_valid_operand(self, other):
        return isinstance(other, self.__class__)

    def __eq__(self, other):
        """Test equality of two `Atom` object instances."""
        if not self._is_valid_operand(other):
            return NotImplemented
        return (self is other or
                (np.allclose([self.Z, self.mass], [other.Z, other.mass]) and
                 self.element == other.element))

    def __le__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        for attr in Atom._fields:
            if getattr(self, attr) > getattr(other, attr):
                return False
        return True

    def __lt__(self, other):
        """Test if `self` is *less than* `other`."""
        if not self._is_valid_operand(other):
            return NotImplemented
        # if self.element == other.element == 'X' and \
        #         self.Z == other.Z == 0:
        #     test = self.mass < other.mass
        # else:
        #     test = self.Z < other.Z
        # return test and self.__le__(other)
        for attr in Atom._fields:
            if getattr(self, attr) >= getattr(other, attr):
                return False
        return True

    def __ge__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        for attr in Atom._fields:
            if getattr(self, attr) < getattr(other, attr):
                return False
        return True

    def __gt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        # if self.element == other.element == 'X' and \
        #         self.Z == other.Z == 0:
        #     test = self.mass > other.mass
        # else:
        #     test = self.Z > other.Z
        # return test and self.__ge__(other)
        for attr in Atom._fields:
            if getattr(self, attr) <= getattr(other, attr):
                return False
        return True

    def __dir__(self):
        return ['element', 'Z', 'mass', 'parent']

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
        element_error_msg = 'unrecognized element value: {}'.format(value)
        symbol = None

        if isinstance(value, str) and re.match(r"[\W]+", value):
            raise ValueError(element_error_msg)

        if isinstance(value, numbers.Integral):
            try:
                Z = int(value)
                idx = Z - 1
                symbol = element_symbols[idx]
            except IndexError:
                print(element_error_msg)

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
        """An alias for :attr:`~Atom.mass`."""
        return self.mass

    @m.setter
    def m(self, value):
        self.mass = value

    def getattr(self, attr, default=None, recursive=False):
        """Get atom attribute named `attr`.

        Parameters
        ----------
        attr : str
            Name of attribute
        default : :class:`~python:object`, optional
        recursive : :class:`~python:bool`, optional

        Returns
        -------
        val : :class:`~python:object`

        """
        if recursive:
            attr_list = attr.split('.')
            obj = self
            for attr in attr_list:
                obj = getattr(obj, attr, default)
            return obj
        else:
            return getattr(self, attr, default)

    def rezero(self, *args, **kwargs):
        assert not hasattr(super(), 'rezero')

    def reset_attrs(self, **kwargs):
        """Reset atom attributes."""
        assert not hasattr(super(), 'reset_attrs')

    def update_attrs(self, **kwargs):
        """Update atom attributes."""
        assert not hasattr(super(), 'update_attrs')

    def todict(self):
        """Return :class:`~python:dict` of `Atom` constructor parameters."""
        return dict(element=self.element, mass=self.mass, Z=self.Z,
                    parent=self.parent)


class Atoms(AtomsSelectionMixin, TabulateMixin, UserList):
    """Base class for collection of `Atom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instance objects or an
        existing `Atoms` instance object.

    """
    def __init__(self, atoms=None, update_item_class=True, **kwargs):
        verbose = kwargs.get('verbose', False)
        if atoms is not None and \
                (isinstance(atoms, str) or isinstance(atoms, Atom)):
            atoms = [atoms]
        # if update_item_class and not isinstance(atoms, type(self)) and \
        #         isinstance(atoms, list):
        if update_item_class and isinstance(atoms, Iterable) and \
                len(atoms) > 0 and not \
                isinstance(atoms[0], self.__atom_class__):
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

    @property
    def __atom_class__(self):
        return Atom

    @property
    def __item_class__(self):
        return self.__atom_class__

    def __str__(self):
        strrep = self._table_title_str()
        objstr = self._obj_mro_str()
        if self.data:
            items = ('Natoms', 'centroid', 'center_of_mass')
            values = [self.Natoms, self.centroid, self.center_of_mass]
            table = self._tabulate(list(zip(items, values)))
            strrep = '\n'.join((strrep, objstr, table))
        return strrep

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return (self is other or self.data == other.data)

    def __lt__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.data < other.data

    def __le__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.data <= other.data

    def __ne__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.data != other.data

    def __gt__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.data > other.data

    def __ge__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.data >= other.data

    def sort(self, key=attrgetter('element', 'Z', 'mass'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @classmethod
    def _from_iterable(cls, it, **kwargs):
        return cls(atoms=it, update_item_class=False, **kwargs)

    def _is_valid_operand(self, other):
        return isinstance(other, (self.__class__, self.__atom_class__))

    def __cast(self, other):
        if not isinstance(other, self.__class__):
            if not isinstance(other, Iterable):
                other = [self.__cast_item(other)]
            other = self._from_iterable(other, **self.kwargs)
        return other

    def __cast_item(self, item):
        if not isinstance(item, self.__item_class__):
            try:
                itemdict = item.todict()
                item = self.__item_class__(
                    **{k: itemdict[k] for k in set(dir(item)) &
                       set(dir(self.__item_class__()))})
            except AttributeError:
                item = self.__item_class__(item)
        return item

    def __setitem__(self, index, item):
        if not self._is_valid_operand(item):
            if isinstance(index, slice):
                item = self.__cast(item)
            else:
                # try:
                #     # atomdict = item.todict()
                #     # atomdict = super().__getitem__(index).todict()
                #     if isinstance(item, str):
                #         atomdict['element'] = item
                #     elif isinstance(item, int):
                #         atomdict['Z'] = item
                #     elif isinstance(item, float):
                #         atomdict['mass'] = item
                #     else:
                #         raise ValueError
                #     item = self.__atom_class__(**atomdict)
                # except (IndexError, AttributeError, ValueError):
                #     item = self.__atom_class__(item)
                item = self.__cast_item(item)
        super().__setitem__(index, item)

    def append(self, atom):
        if not self._is_valid_operand(atom):
            atom = self.__cast_item(atom)
        super().append(atom)

    def insert(self, i, atom):
        if not self._is_valid_operand(atom):
            atom = self.__cast_item(atom)
        super().insert(i, atom)

    def __add__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        atoms = self.data + \
            list(atom for atom in self.__cast(other) if atom not in self)
        return self._from_iterable(atoms, **self.kwargs)

    def __radd__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        other = self.__cast(other)
        atoms = other.data + list(atom for atom in self if atom not in other)
        return self._from_iterable(atoms, **self.kwargs)

    def __iadd__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        self.data += list(atom for atom in self.__cast(other)
                          if atom not in self)
        return self

    def __sub__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self._from_iterable((atom for atom in self
                                    if atom not in self.__cast(other)),
                                   **self.kwargs)

    def __rsub__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self._from_iterable((atom for atom in self.__cast(other)
                                    if atom not in self),
                                   **self.kwargs)

    def __isub__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        self.data -= list(atom for atom in self.__cast(other) if atom in self)
        return self

    def __and__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self._from_iterable((atom for atom in self.__cast(other)
                                    if atom in self), **self.kwargs)

    __rand__ = __and__

    def __iand__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        self.data -= list(atom for atom in self.__cast(other)
                          if atom not in self)
        return self

    def __or__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        # atoms = self.data + list(atom for atom in other if atom not in self)
        atoms = dedupe((atom for atoms in (self, self.__cast(other))
                        for atom in atoms), key=attrgetter('id'))
        return self._from_iterable(atoms, **self.kwargs)

    __ror__ = __or__

    def __ior__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        self.data += \
            list(atom for atom in self.__cast(other) if atom not in self)
        return self

    def __xor__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        other = self.__cast(other)
        return (self - other) | (other - self)

    __rxor__ = __xor__

    def __ixor__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if other is self:
            self.clear()
        else:
            other = self.__cast(other)
            for atom in other:
                if atom in self:
                    self.remove(atom)
                else:
                    self.append(atom)
        return self

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
        """:class:`~numpy:numpy.ndarray` of :attr:`Atom.element`\ s."""
        return np.asarray([atom.element for atom in self])

    @property
    def masses(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`Atom.mass`\ s."""
        return np.asarray([atom.mass for atom in self])

    @property
    def symbols(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`Atom.symbol`\ s."""
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

    def get_atoms(self, asarray=False, aslist=True):
        """Return `Atoms` either as list (default) or numpy array or self.

        Parameters
        ----------
        asarray, aslist : :class:`~python:bool`, optional
            Default values: `asarray=False`, `aslist=True`

        Returns
        -------
        :class:`~numpy:numpy.ndarray` if `asarray` is `True`
        :class:`~python:list` if `asarray` is `False` and `aslist` is `True`
        :class:`Atoms` object if `asarray` and `aslist` are `False`

        """
        if asarray:
            return np.asarray(self.data)
        elif aslist:
            return self.data
        else:
            return self

    def getattr(self, attr, default=None, recursive=False):
        """Get :class:`~numpy:numpy.ndarray` of atom attributes `attr`.

        Parameters
        ----------
        attr : str
            Name of attribute to pass to `getattr` from each `Atom` in
            `Atoms`.
        default : :class:`~python:object`, optional
        recursive : :class:`~python:bool`, optional

        Returns
        -------
        :class:`~numpy:numpy.ndarray`

        """
        if recursive:
            attr_values = []
            attr_list = attr.split('.')
            for atom in self:
                obj = atom
                for attr in attr_list:
                    obj = getattr(obj, attr, default)
                attr_values.append(obj)
            return np.asarray(attr_values)
        else:
            return np.asarray([getattr(atom, attr, default) for atom in self])

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
        Suppose you have an `Atoms` instance named ``atoms`` that contains
        `Atom` instances of two atom types `1` and `2` and we want to set
        all `Atom`\ s with `type=1` to Nitrogen and all `Atom`\ s with
        `type=2` to Argon. In other words, we want to use the
        :attr:`~sknano.core.atoms.TypeAtom.type` attribute to set the
        :attr:`~Atom.element` attribute, which we do by passing a
        :class:`~python:dict` mapping each `type` to the respective element
        symbol. For example::

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

    def reset_attrs(self, **kwargs):
        """Call corresponding `reset_attrs` method on each atom"""
        # [atom.reset_attrs(**kwargs) for atom in self]
        assert not hasattr(super(), 'reset_attrs')

    def update_attrs(self, **kwargs):
        """Call `update_attrs` method on each atom."""
        # [atom.update_attrs() for atom in self]
        # if len(kwargs) != 0:
        #     warnings.warn('`Atoms.update_attrs` received unused kwargs: \n'
        #                   '{}'.format(kwargs))
        assert not hasattr(super(), 'update_attrs')


def update_atoms(atoms, kwargs, deepcopy=True, update_kwargs=False):
    if not update_kwargs:
        kwargs = kwargs.copy()

    atoms = atoms[:]
    if deepcopy:
        atoms = copy.deepcopy(atoms)

    if any([kw in kwargs for kw
            in ('center_CM', 'center_center_of_mass')]):
        center_com = \
            kwargs.pop('center_CM', kwargs.pop('center_center_of_mass'))

    region_bounds = kwargs.pop('region_bounds', None)
    center_centroid = kwargs.pop('center_centroid', True)
    center_com = kwargs.pop('center_com', False)
    filter_condition = kwargs.pop('filter_condition', None)
    rotation_parameters = kwargs.pop('rotation_parameters', None)

    if region_bounds is not None:
        atoms.clip_bounds(region_bounds)

    if center_centroid:
        atoms.center_centroid()
    elif center_com:
        atoms.center_com()

    if filter_condition is not None:
        atoms.filter(filter_condition)
        # atoms = atoms.filtered(filter_condition)

    rotation_kwargs = ['rotation_angle', 'angle', 'rot_axis', 'axis',
                       'anchor_point', 'deg2rad', 'degrees', 'rot_point',
                       'from_vector', 'to_vector', 'transform_matrix']

    if rotation_parameters is None and \
            any([kw in kwargs for kw in rotation_kwargs]):
        rotation_parameters = {kw: kwargs.pop(kw) for kw in rotation_kwargs
                               if kw in kwargs}
        if 'rotation_angle' in rotation_parameters:
            rotation_parameters['angle'] = \
                rotation_parameters.pop('rotation_angle')
        if 'rot_axis' in rotation_parameters:
            rotation_parameters['axis'] = rotation_parameters.pop('rot_axis')
        if 'deg2rad' in rotation_parameters:
            rotation_parameters['degrees'] = rotation_parameters.pop('deg2rad')

    if rotation_parameters is not None and \
            isinstance(rotation_parameters, dict):
        atoms.rotate(**rotation_parameters)

    atoms.rezero()
    return atoms
