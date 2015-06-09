# -*- coding: utf-8 -*-
"""
===============================================================================
Atoms class for POAV analysis (:mod:`sknano.core.atoms._poav_atoms`)
===============================================================================

An `Atoms` class for POAV analysis

.. currentmodule:: sknano.core.atoms._poav_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numpy as np

# from sknano.core import timethis
from sknano.core.math import vector as vec
from ._kdtree_atoms import KDTAtoms
from ._poav_atom import POAV1, POAV2, POAVR, POAVAtom

__all__ = ['POAVAtoms']


class POAVAtoms(KDTAtoms):
    """An `Atoms` sub-class for POAV analysis.

    Sub-class of `KDTAtoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.POAVAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `POAVAtoms`}, optional
        if not `None`, then a list of `POAVAtom` instance objects or an
        existing `POAVAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return POAVAtom

    # @timethis
    def compute_POAVs(self):
        """Compute `POAV1`, `POAV2`, `POAVR`."""
        super().update_attrs()

        POAV_classes = {'POAV1': POAV1, 'POAV2': POAV2, 'POAVR': POAVR}

        for atom in self:
            # the central atom must have 3 bonds for POAV analysis.
            if atom.bonds.Nbonds == 3:
                for POAV_name, POAV_class in list(POAV_classes.items()):
                    setattr(atom, POAV_name, POAV_class(atom.bonds))

        for atom in self:
            # the central atom must have 3 bonds for POAV analysis.
            if atom.bonds.Nbonds == 3:
                for POAV_name in ('POAV1', 'POAV2', 'POAVR'):
                    POAV = getattr(atom, POAV_name)
                    sigma_pi_angles = []
                    pyramidalization_angles = []
                    misalignment_angles = []
                    for bond, NN in zip(atom.bonds, atom.NN):
                        # first compute the pyramidalization angle
                        sigma_pi_angle = vec.angle(POAV.Vpi, bond.vector)
                        if sigma_pi_angle < np.pi / 2:
                            sigma_pi_angle = np.pi - sigma_pi_angle
                        sigma_pi_angles.append(sigma_pi_angle)
                        pyramidalization_angles.append(
                            sigma_pi_angle - np.pi / 2)

                        # the bonded atom must have a POAV to compute the
                        # misalignment angles
                        if getattr(NN, POAV_name) is not None:
                            NN_POAV = getattr(NN, POAV_name)
                            # compute vector that is orthogonal to the plane
                            # defined by the bond vector and the POAV of the
                            # center atom.
                            nvec = vec.cross(bond.vector, POAV.Vpi)

                            # the misalignment angle is the angle between the
                            # nearest neighbor's POAV and the plane defined by
                            # the bond vector and the POAV of the center atom,
                            # which is pi/2 minus the angle between
                            # the NN POAV and the normal vector to the plane
                            # computed above.
                            misalignment_angles.append(np.abs(
                                np.pi / 2 - vec.angle(NN_POAV.Vpi, nvec)))
                        else:
                            misalignment_angles.append(np.nan)

                    POAV.pyramidalization_angles = pyramidalization_angles
                    POAV.misalignment_angles = misalignment_angles
                    POAV.sigma_pi_angles = sigma_pi_angles

    @property
    def POAV1(self):
        """List of :class:`~sknano.core.atoms.POAVAtom` :class:`POAV1` \
            :attr:`~sknano.core.atoms.POAVAtom.POAV1` attribute."""
        return [atom.POAV1 for atom in self if atom.POAV1 is not None]

    @property
    def POAV2(self):
        """List of :class:`~sknano.core.atoms.POAVAtom` :class:`POAV2` \
            :attr:`~sknano.core.atoms.POAVAtom.POAV2` attribute."""
        return [atom.POAV2 for atom in self if atom.POAV2 is not None]

    @property
    def POAVR(self):
        """List of :class:`~sknano.core.atoms.POAVAtom` :class:`POAVR` \
            :attr:`~sknano.core.atoms.POAVAtom.POAVR` attribute."""
        return [atom.POAVR for atom in self if atom.POAVR is not None]

    def get_POAV_attr(self, POAV_class, attr):
        """Return list of :class:`~sknano.core.atoms.POAVAtom` :class:`POAV1` \
            :class:`POAV2` or :class:`POAVR` attribute.

        Parameters
        ----------
        POAV_class : :class:`~python:str`
        attr : :class:`~python:str`

        Returns
        -------
        :class:`~python:list`

        """
        return [getattr(getattr(atom, POAV_class), attr) for atom in self
                if getattr(atom, POAV_class) is not None]
