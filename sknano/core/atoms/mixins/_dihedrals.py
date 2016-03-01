# -*- coding: utf-8 -*-
"""
===============================================================================
Class representations of dihedral angles (:mod:`sknano.core.atoms.mixins._dihedrals`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins._dihedrals

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core.math import vector as vec
from ._bonds import Bond
from ._topology_base import AtomTopology, AtomsTopology, \
    check_operands as check_operands_


__all__ = ['compute_dihedral', 'Dihedral', 'Dihedrals']


def compute_dihedral(*atoms, check_operands=True, degrees=False):
    """Compute dihedral angle.

    Parameters
    ----------
    *atoms : {:class:`~python:list`, :class:`~sknano.core.atoms.Atoms`}
        :class:`~python:list` of :class:`~sknano.core.atoms.Atom`\ s
        or an :class:`~sknano.core.atoms.Atoms` object.
    check_operands : :class:`~python:bool`, optional
    degrees : :class:`~python:bool`, optional

    Returns
    -------
    :class:`~python:float`

    Raises
    ------
    :class:`~python:TypeError`
        if `atoms` is not a list of :class:`~sknano.core.atoms.Atom` objects
        or an :class:`~sknano.core.atoms.Atoms` object.
    :class:`~python:ValueError`
        if len(atoms) != 4.

    """
    if check_operands:
        atoms = check_operands_(*atoms, size=4)

    atom1, atom2, atom3, atom4 = atoms
    b12 = Bond(atom1, atom2).vector
    b23 = Bond(atom2, atom3).vector
    b34 = Bond(atom3, atom4).vector
    m = b12.cross(b23)
    n = b23.cross(b34)
    sina = vec.dot(m, b34) * b23.norm
    cosa = vec.dot(m, n)
    angle = np.arctan2(sina, cosa)
    if degrees:
        angle = np.degrees(angle)
    return angle


class Dihedral(AtomTopology):
    """Class representation of dihedral angle between 4 `Atom` objects.

    Parameters
    ----------
    *atoms : {:class:`~python:list`, :class:`~sknano.core.atoms.Atoms`}
        :class:`~python:list` of :class:`~sknano.core.atoms.Atom`\ s
        or an :class:`~sknano.core.atoms.Atoms` object.
    check_operands : :class:`~python:bool`, optional
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.
    id : :class:`~python:int`

    Raises
    ------
    :class:`~python:TypeError`
        if `atoms` is not a list of :class:`~sknano.core.atoms.Atom` objects
        or an :class:`~sknano.core.atoms.Atoms` object.
    :class:`~python:ValueError`
        if len(atoms) != 4.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if self.check_operands and len(self.atoms) != 4:
                raise ValueError('Expected 4 atoms')

        self.fmtstr = "{atom1!r}, {atom2!r}, {atom3!r}, {atom4!r} " + \
            super().fmtstr

    def __str__(self):
        """Return nice string representation of `Dihedral`."""
        return "Dihedral({!r}->{!r}->{!r}->{!r})".format(self.atom1.id,
                                                         self.atom2.id,
                                                         self.atom3.id,
                                                         self.atom4.id)

    @property
    def atom1(self):
        """:class:`~sknano.core.atoms.Atom` 1 in `Dihedral`."""
        return self.atoms[0]

    @property
    def atom2(self):
        """:class:`~sknano.core.atoms.Atom` 2 in `Dihedral`."""
        return self.atoms[1]

    @property
    def atom3(self):
        """:class:`~sknano.core.atoms.Atom` 3 in `Dihedral`."""
        return self.atoms[2]

    @property
    def atom4(self):
        """:class:`~sknano.core.atoms.Atom` 4 in `Dihedral`."""
        return self.atoms[3]

    @property
    def angle(self):
        """An alias for :attr:`AtomTopology.measure`."""
        return self.measure

    def compute_measure(self):
        """Compute the bond angle, which is the measure of an :class:`Angle`.

        Returns
        -------
        :class:`~python:float`

        """
        return compute_dihedral(*self.atoms)

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        super_dict = super().todict()
        super_dict.update(dict(atom1=self.atom1, atom2=self.atom2,
                               atom3=self.atom3, atom4=self.atom4))
        return super_dict


class Dihedrals(AtomsTopology):
    """Base class for collection of atom `Dihedral`\ s.

    Parameters
    ----------
    topolist : {None, sequence, `Dihedrals`}, optional
        if not `None`, then a list of `Dihedral` objects
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.

    """
    @property
    def __item_class__(self):
        return Dihedral

    @property
    def Ndihedrals(self):
        """Number of `Dihedral`\ s in `Dihedrals`."""
        return len(self)

    @property
    def angles(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`~Dihedral.angles`\ s."""
        return self.measures

    @property
    def mean_angle(self):
        """Mean dihedral angle."""
        return self.mean_measure
