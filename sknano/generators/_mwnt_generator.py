# -*- coding: utf-8 -*-
"""
===============================================================================
MWNT structure generator (:mod:`sknano.generators._mwnt_generator`)
===============================================================================

.. currentmodule:: sknano.generators._mwnt_generator

.. todo::

   Add methods to perform fractional translation and cartesian translation
   before structure generation.

.. todo::

   Handle different units in output coordinates.

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

#import copy

#import numpy as np

#from sknano.core import pluralize
#from sknano.core.math import Vector
from sknano.structures import MWNT
#from sknano.utils.geometric_shapes import Cuboid
from ._base import Atoms, GeneratorBase
from ._swnt_generator import SWNTGenerator

__all__ = ['MWNTGenerator']


class MWNTGenerator(MWNT, GeneratorBase):
    """Class for generating single, multi-walled nanotubes.

    .. versionchanged:: 0.2.20

       `MWNTGenerator` no longer generates MWNT *bundles*, only *single*
       MWNTs. To generate bundled MWNT structure data, use the
       `MWNTBundleGenerator` class.

    .. versionadded:: 0.2.8

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.Atom` 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    Lx, Ly, Lz : float, optional
        length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.
    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.
    add_outer_shells : bool, optional
        Build the MWNT by adding outer shells

        .. versionadded:: 0.2.23

    add_inner_shells : bool, optional
        Build the MWNT by adding inner shells

        .. versionadded:: 0.2.23

    max_shells : int, optional
        Maximum number of shells per MWNT.
    max_shell_diameter : float, optional
        Maximum shell diameter, in units of **Angstroms**.
    min_shells : int, optional
        Minimum number of shells per MWNT.
    min_shell_diameter : float, optional
        Minimum shell diameter, in units of **Angstroms**.
    new_shell_type : {None, 'AC', 'ZZ', 'achiral', 'chiral'}, optional
        If `None`, the chiralities of the new shells are constrained only
        by their diameter and will be chosen randomly if more than one
        candidate chirality exists. If not `None`, then the
        `new_shell_type` will be added as a constraint.
    shell_spacing : float, optional
        Shell spacing in units of **Angstroms**. Default
        value is the van der Waals interaction distance of 3.4 Angstroms.
    autogen : bool, optional
        if `True`, automatically call
        :meth:`~MWNTGenerator.generate_unit_cell`,
        followed by :meth:`~MWNTGenerator.generate_structure_data`.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.generators import MWNTGenerator
    >>> mwnt = MWNTGenerator(n=40, m=40, max_shells=5, Lz=1.0, fix_Lz=True)
    >>> mwnt.save_data()

    .. image:: /images/5shell_mwnt_4040_outer_Ch_1cellx1cellx4.06cells-01.png

    """
    def __init__(self, autogen=True, **kwargs):

        super(MWNTGenerator, self).__init__(**kwargs)

        if autogen:
            self.generate_unit_cell()
            self.generate_structure_data()

    def generate_unit_cell(self):
        """Generate the `MWNT` unit cell."""
        self.unit_cell = Atoms()
        for swnt in self.shells:
            self.unit_cell.extend(SWNTGenerator(**swnt.todict()).unit_cell)

    def generate_structure_data(self):
        """Generate structure data.

        .. todo::

           Load the diameter and chirality data from file instead of
           generating it every time.

        """

        #self.atoms = Atoms(atoms=self.unit_cell)
        #self.atoms = Atoms()
        self.structure_data.clear()
        for swnt in self.shells:
            self.atoms.extend(SWNTGenerator(**swnt.todict()).atoms)
        #for n, m in self.Ch:
        #for nz in xrange(int(np.ceil(self.nz))):
        #    dr = Vector([0.0, 0.0, nz * self.T])
        #    for uc_atom in self.unit_cell:
        #        nt_atom = Atom(element=uc_atom.symbol)
        #        nt_atom.r = uc_atom.r + dr
        #        self.atoms.append(nt_atom)

        #Lzmin = self.atoms
        #pmin = [-np.inf, -np.inf, -np.inf]
        #pmax = [np.inf, np.inf, np.inf]
        #region_bounds = Cuboid(pmin=pmin, pmax=pmax)
        #if self.L0 is not None and self.fix_Lz:
        #    region_bounds.zmax = (10 * self.L0 + 1) / 2
        #else:
        #    region_bounds.zmax = (10 * Lzmin + 1) / 2

        #region_bounds.zmin = -region_bounds.zmax
        #region_bounds.update_region_limits()

        #self.atoms.clip_bounds(region_bounds,
        #                                 center_before_clipping=True)

    def save_data(self, fname=None, outpath=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  anchor_point=None, center_CM=True, savecopy=True, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save_data` method
        for documentation.

        """
        if fname is None:
            Nshells = '{}shell_mwnt'.format(self.Nshells)
            chiralities = '@'.join([str(Ch).replace(' ', '') for
                                    Ch in self.Ch])

            fname_wordlist = None
            #if self._assert_integer_nz:
            #    nz = ''.join(('{}'.format(self.nz),
            #                  pluralize('cell', self.nz)))
            #else:
            #    nz = ''.join(('{:.2f}'.format(self.nz),
            #                  pluralize('cell', self.nz)))
            #fname_wordlist = (Nshells, chiralities, nz)
            fname_wordlist = (Nshells, chiralities)
            fname = '_'.join(fname_wordlist)

        super(MWNTGenerator, self).save_data(
            fname=fname, outpath=outpath, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            anchor_point=anchor_point, deg2rad=deg2rad, center_CM=center_CM,
            savecopy=savecopy, **kwargs)
