# -*- coding: utf-8 -*-
"""
===============================================================================
Defect generator base classes (:mod:`sknano.generators._defect_generators`)
===============================================================================

.. currentmodule:: sknano.generators._defect_generators

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import os

import numpy as np

from ._base import Atoms, GeneratorBase

__all__ = ['DefectGenerator']
# __all__ += ['CrossLinkedDefectGenerator']
# __all__ += ['StoneWalesDefectGenerator']
# __all__ += ['VacancyDefectGenerator']


class DefectGenerator(GeneratorBase):
    """Base class for generating structure data with defects.

    Parameters
    ----------
    fname : str
        structure data filename
    outpath : str, optional
        Output path for structure data file.
    structure_format : {None, str}, optional
        chemical file format of saved structure data.
        If `None`, then guess based on `fname` file extension.
        Otherwise, must be one of:

            - `xyz`
            - `data`

    verbose : bool, optional
        Verbose output

    """
    pass
#     def __init__(self, fname, outpath=None, structure_format=None,
#                  verbose=False):

#         if fname.endswith(supported_structure_formats) and \
#                 structure_format is None:
#             for ext in supported_structure_formats:
#                 if fname.endswith(ext):
#                     structure_format = ext
#                     break
#         else:
#             if (not fname.endswith(supported_structure_formats) and
#                     structure_format is None) or \
#                     (structure_format is not None and
#                         structure_format not in supported_structure_formats):
#                 raise ValueError(
#                     '{} is not a supported structure format'.format(
#                         structure_format))

#         self._fname = fname
#         self._outpath = outpath
#         self._structure_format = structure_format
#         self._verbose = verbose

#         # parse structure data
#         self._structure_data = \
#             StructureReader.read(fname, structure_format=self._structure_format)
#         self._atoms = self._structure_data.atoms
#         self._atom_ids = self._atoms.atom_ids
#         self._atom_coords = self._atoms.get_coords(asdict=True)

#     @property
#     def atoms(self):
#         return self._atoms

#     @property
#     def atom_ids(self):
#         return self._atom_ids

#     @property
#     def atom_coords(self):
#         return self._atom_coords


# class CrossLinkedDefectGenerator(DefectGenerator):
#     """Cross-linked defect generator class."""
#     pass


# class StoneWalesDefectGenerator(DefectGenerator):
#     """Stone-Wales defect generator class.

#     Parameters
#     ----------
#     fname : str
#         structure data filename
#     outpath : str, optional
#         Output path for structure data file.
#     structure_format : {None, str}, optional
#         chemical file format of saved structure data.
#         If `None`, then guess based on `fname` file extension.
#         Otherwise, must be one of:

#             - `xyz`
#             - `data`

#     verbose : bool, optional
#         Verbose output

#     """
#     def __init__(self, fname=str, outpath=None, structure_format=None,
#                  verbose=False):
#         super(StoneWalesDefectGenerator, self).__init__(
#             fname=fname, outpath=outpath, structure_format=structure_format,
#             verbose=verbose)

#         self._defect_bounds = None

#         self._vmd_selection_radius = np.sqrt(10.5)
#         self._show_vmd_selection_cmd = True


# class VacancyGenerator(DefectGenerator):
#     """Class for generating vacancies in structure data.

#     Parameters
#     ----------
#     fname : str
#         structure data filename
#     outpath : str, optional
#         Output path for structure data file.
#     structure_format : {None, str}, optional
#         chemical file format of saved structure data.
#         If `None`, then guess based on `fname` file extension.
#         Otherwise, must be one of:

#             - `xyz`
#             - `data`

#     verbose : bool, optional
#         Verbose output

#     """
#     def __init__(self, fname=str, outpath=None, structure_format=None,
#                  verbose=False):
#         super(VacancyGenerator, self).__init__(
#             fname=fname, outpath=outpath, structure_format=structure_format,
#             verbose=verbose)

#         self._Nvacs = 0
#         self._Nvac_clusters = 0
#         self._Nvac_sites = 0

#         self._vac_bounds = None

#         self._vac_ids = np.empty(0, dtype=int)
#         self._vac_type = 'single'
#         self._cluster_size = 1

#         self._vmd_selection_radius = np.sqrt(10.5)
#         self._show_vmd_selection_cmd = True

#     @property
#     def Nvacs(self):
#         return self._Nvacs

#     @Nvacs.setter
#     def Nvacs(self, value=int):
#         self._Nvacs = value

#     @property
#     def Nvac_clusters(self):
#         return self._Nvac_clusters

#     @Nvac_clusters.setter
#     def Nvac_clusters(self, value=int):
#         self._Nvac_clusters = value

#     @property
#     def Nvac_sites(self):
#         return self._Nvac_sites

#     @Nvac_sites.setter
#     def Nvac_sites(self, value=int):
#         self._Nvac_sites = value

#     @property
#     def cluster_size(self):
#         return self._cluster_size

#     @cluster_size.setter
#     def cluster_size(self, value):
#         self._cluster_size = value

#     @property
#     def vac_bounds(self):
#         return self._vac_bounds

#     @vac_bounds.setter
#     def vac_bounds(self, value):
#         self._vac_bounds = value

#     @property
#     def vac_ids(self):
#         return self._vac_ids

#     @property
#     def vac_type(self):
#         return self._vac_type

#     @vac_type.setter
#     def vac_type(self, value):
#         self._vac_type = value

#     @property
#     def vmd_selection_radius(self):
#         return self._vmd_selection_radius

#     @vmd_selection_radius.setter
#     def vmd_selection_radius(self, value):
#         self._vmd_selection_radius = value

#     @property
#     def show_vmd_selection_cmd(self):
#         return self._show_vmd_selection_cmd

#     @show_vmd_selection_cmd.setter
#     def show_vmd_selection_cmd(self, value=bool):
#         self._show_vmd_selection_cmd = value

#     def _random_vacancy_generator(self):
#         """Generate random vacancies in structure data."""
#         self._vac_ids = \
#             np.random.choice(self._atom_ids,
#                              size=self._Nvac_sites,
#                              replace=False)

#     def _generate_vmd_selection_cmd(self):

#         selection_radius = self._vmd_selection_radius
#         selections = []
#         for atom in self._removed_atoms:
#             selection_cmd = \
#                 "(((x-{:.4f})^2 + ".format(atom.x) + \
#                 "(y-{:.4f})^2 + ".format(atom.y) + \
#                 "(z-{:.4f})^2) <= {:.2f})".format(atom.z,
#                                                   selection_radius**2)
#             selections.append(selection_cmd)

#         vmd_selection_cmd = ' or '.join(selections)
#         print('copy and paste the following VMD command to select\n'
#               'the atoms surrounding the vacancies:\n\n'
#               '{}\n'.format(vmd_selection_cmd))

#     def _generate_single_vacancies(self):
#         self._removed_atoms = \
#             self._atoms.filter_ids(self._vac_ids, invert=False)

#     def _generate_multi_vacancies(self):
#         vac_type_properties = {'double': {'cluster_size': 2,
#                                           'NN_cutoff': 1.5},
#                                'triple': {'cluster_size': 3,
#                                           'NN_cutoff': 1.5}}
#         vac_props = vac_type_properties[self._vac_type]
#         self._cluster_size = vac_props['cluster_size']
#         self._atoms.NN_cutoff = vac_props['NN_cutoff']
#         self._atoms.update_neighbors()

#         vac_atoms = Atoms()
#         for vac_id in self._vac_ids:
#             vac_atom = self._atoms.get_atom(vac_id)
#             vac_atoms.append(vac_atom)
#             vac_atoms.extend(np.random.choice(vac_atom.NN,
#                                               size=self._cluster_size-1,
#                                               replace=False).tolist())
#         self._removed_atoms = \
#             self._atoms.filter_ids(vac_atoms.atom_ids, invert=False)

#     def _generate_vacancy_structure(self):
#         """Generate vacancy structure."""
#         if self._vac_type in ('double', 'triple'):
#             self._generate_multi_vacancies()
#         else:
#             self._generate_single_vacancies()

#         if self._show_vmd_selection_cmd:
#             self._generate_vmd_selection_cmd()

#         self._Nvacs = self._removed_atoms.Natoms

#         self._remaining_atoms = \
#             self._atoms.filter_ids(self._removed_atoms.atom_ids, invert=True)
#         #remaining_atoms.assign_unique_ids()

#         self._save_vacancy_structure_data()

#     def _generate_output_fname(self):
#         self._output_fname = \
#             os.path.splitext(os.path.basename(self._fname))[0] + \
#             '+{}_vacancies'.format(self._Nvacs)

#     def _save_vacancy_structure_data(self):
#         self._generate_output_fname()
#         DATAWriter.write(fname=self._output_fname, outpath=self._outpath,
#                          atoms=self._remaining_atoms,
#                          boxbounds=self._structure_data.boxbounds,
#                          comment_line=self._structure_data.comment_line)
#         XYZWriter.write(fname=self._output_fname, outpath=self._outpath,
#                         atoms=self._remaining_atoms,
#                         comment_line=self._structure_data.comment_line)
