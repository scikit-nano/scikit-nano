# -*- coding: utf-8 -*-
"""
=======================================================================================
Class representations of improper angles (:mod:`sknano.core.atoms.mixins._impropers`)
=======================================================================================

.. currentmodule:: sknano.core.atoms.mixins._impropers

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import namedtuple
import numpy as np

# from sknano.core.math import vector as vec
from sknano.core.math import abs_cap
from ._bonds import Bond
from ._topology_base import AngularTopology, AngularTopologyCollection, \
    TopologyStats, check_operands as check_operands_

__all__ = ['compute_improper', 'Improper', 'Impropers', 'ImproperStats']

ImproperStats = namedtuple('ImproperStats', TopologyStats._fields)


def compute_improper(*atoms, check_operands=True, degrees=False):
    """Compute improper angle.

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
    b21 = Bond(atom2, atom1).vector
    b23 = Bond(atom2, atom3).vector
    b34 = Bond(atom3, atom4).vector

    m = b21.cross(b23)
    n = b23.cross(b34)
    cosa = -m.dot(n) / (m.norm * n.norm)
    cosa = abs_cap(cosa, 1)
    angle = np.arccos(cosa)
    if degrees:
        angle = np.degrees(angle)

    return angle


class Improper(AngularTopology):
    """Class representation of improper angle between 4 `Atom` objects.

    Parameters
    ----------
    *atoms : {:class:`~python:list`, :class:`~sknano.core.atoms.Atoms`}
        :class:`~python:list` of :class:`~sknano.core.atoms.Atom`\ s
        or an :class:`~sknano.core.atoms.Atoms` object.
    size : :class:`~python:int`
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.
    id : :class:`~python:int`
    check_operands : :class:`~python:bool`, optional
    degrees : :class:`~python:bool`, optional

    Raises
    ------
    :class:`~python:TypeError`
        if `atoms` is not a list of :class:`~sknano.core.atoms.Atom` objects
        or an :class:`~sknano.core.atoms.Atoms` object.
    :class:`~python:ValueError`
        if len(atoms) != 4.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, size=4, **kwargs)
        self.fmtstr = "{atom1!r}, {atom2!r}, {atom3!r}, {atom4!r} " + \
            super().fmtstr

    def __str__(self):
        """Return nice string representation of `Improper`."""
        return "Improper({!r}->{!r}->{!r}->{!r})".format(self.atom1.id,
                                                         self.atom2.id,
                                                         self.atom3.id,
                                                         self.atom4.id)

    @property
    def atom1(self):
        """:class:`~sknano.core.atoms.Atom` 1 in `Improper`."""
        return self.atoms[0]

    @property
    def atom2(self):
        """:class:`~sknano.core.atoms.Atom` 2 in `Improper`."""
        return self.atoms[1]

    @property
    def atom3(self):
        """:class:`~sknano.core.atoms.Atom` 3 in `Improper`."""
        return self.atoms[2]

    @property
    def atom4(self):
        """:class:`~sknano.core.atoms.Atom` 4 in `Improper`."""
        return self.atoms[3]

    @property
    def improper(self):
        """An alias for :attr:`Topology.measure`."""
        return self.measure

    def compute_measure(self):
        """Compute the bond angle, which is the measure of an :class:`Angle`.

        Returns
        -------
        :class:`~python:float`

        """
        return compute_improper(*self.atoms, degrees=self.degrees)

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        super_dict = super().todict()
        super_dict.update(dict(atom1=self.atom1, atom2=self.atom2,
                               atom3=self.atom3, atom4=self.atom4))
        return super_dict


class Impropers(AngularTopologyCollection):
    """Base class for collection of atom `Improper`\ s.

    Parameters
    ----------
    topolist : {None, sequence, `Impropers`}, optional
        if not `None`, then a list of `Improper` objects
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.
    degrees : :class:`~python:bool`, optional

    """
    @property
    def __item_class__(self):
        return Improper

    @property
    def Nimpropers(self):
        """Number of `Improper`\ s in `Impropers`."""
        return len(self)

    @property
    def impropers(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`~Improper.angles`\ s."""
        return self.measures

    @property
    def mean_improper(self):
        """Mean improper angle."""
        return self.mean_measure

    @property
    def statistics(self):
        """Improper stats."""
        return ImproperStats(**super().statistics._asdict())
