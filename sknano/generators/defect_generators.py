# -*- coding: utf-8 -*-
"""
===============================================================================
Defect generator base classes (:mod:`sknano.generators.defect_generators`)
===============================================================================

.. currentmodule:: sknano.generators.defect_generators

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from os.path import basename, splitext

import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage

from sknano.core import pluralize
from .base import DefectGeneratorBase

__all__ = ['CrossLinkedDefectGenerator', 'StoneWalesDefectGenerator',
           'VacancyGenerator']


class CrossLinkedDefectGenerator(DefectGeneratorBase):
    """Cross-linked defect generator class."""
    def generate(self, finalize=True):
        pass


class StoneWalesDefectGenerator(DefectGeneratorBase):
    """Stone-Wales defect generator class.

    Parameters
    ----------

    """
    def __init__(self, *args, **kwargs):
        super(StoneWalesDefectGenerator, self).__init__(*args, **kwargs)

    def generate(self, finalize=True):
        pass


class VacancyGenerator(DefectGeneratorBase):
    """Class for generating vacancies in structure data.

    Parameters
    ----------
    """
    def __init__(self, *args, Nsites=None, Nvacs_per_site=None,
                 uniform=False, region=None, vmd_selection_radius=None,
                 **kwargs):
        if Nsites is None:
            raise TypeError('Expected an integer for `Nsites`.')
        self.Nsites = Nsites

        self.Nvacs_per_site = Nvacs_per_site

        self.uniform = uniform
        self.region = region

        if vmd_selection_radius is None:
            vmd_selection_radius = np.sqrt(10.5)
        self.vmd_selection_radius = vmd_selection_radius

        super(VacancyGenerator, self).__init__(*args, **kwargs)

    def generate(self):
        """Generate vacancy structure."""
        atoms = self.atoms
        atoms.update_neighbors()

        if self.uniform:
            vac_atoms = self.uniform_selection()
        else:
            vac_atoms = self.random_selection()

        Nvacs_per_site = self.Nvacs_per_site
        if Nvacs_per_site is not None:
            for atom in vac_atoms[:]:
                vac_atoms.extend(atom.get_n_neighbors(Nvacs_per_site - 1,
                                                      exclude=vac_atoms))

        # self.removed_atoms = \
        #     self.atoms.filtered_ids(vac_atoms.ids, invert=False)

        # if self.echo_vmd_selection_cmd:
        #     self._generate_vmd_selection_cmd()

        # self.remaining_atoms = \
        #     self.atoms.filtered_ids(self.removed_atoms.ids, invert=True)

    def random_selection(self):
        """Select random distribution of atoms."""
        atoms = self.atoms
        ids = np.random.choice(atoms.ids, size=self.Nsites, replace=False)
        return atoms.filtered_ids(ids, invert=False)

    def uniform_selection(self):
        """Select uniform distribution of atoms."""
        pass

    def _generate_vmd_selection_cmd(self):

        r = self.vmd_selection_radius
        selstrs = []
        for atom in self.removed_atoms:
            selstr = "(((x-{:.4f})^2 + ".format(atom.x) + \
                "(y-{:.4f})^2 + ".format(atom.y) + \
                "(z-{:.4f})^2) <= {:.2f})".format(atom.z, r**2)
            selstrs.append(selstr)

        selection_cmd = ' or '.join(selstrs)
        print('copy and paste the following VMD command to select\n'
              'the atoms surrounding the vacancies:\n\n'
              '{}\n'.format(selection_cmd))

    def generate_fname(self, Nvacs=None):
        fname = None
        structure_format = None
        if isinstance(self.structure_or_file, str):
            fname, structure_format = \
                splitext(basename(self.structure_or_file))
        else:
            try:
                fname = self.structure.__class__.__qualname__
            except AttributeError:
                fname = 'VacancyGenerator.structure'

                '+{}_vacancies'.format(self._Nvacs)

        if Nvacs is not None:
            Nvacs = '_'.join((str(Nvacs), pluralize('vacancy', Nvacs)))
            fname = '+'.join((fname, Nvacs))

        if structure_format is not None:
            fname = ''.join((fname, structure_format))

        return fname

    def save(self, fname=None, **kwargs):
        if fname is None:
            fname = self.generate_fname(Nvacs=self.Nvacs)
        super().save(fname, atoms=self.remaining_atoms, **kwargs)
