# -*- coding: utf-8 -*-
"""
=============================================================================
Generate vacancies in structures (:mod:`sknano.nanogen.vacancygenerator`)
=============================================================================

Module for generating vacancies in nano-structures.

.. currentmodule:: sknano.nanogen.vacancygenerator

.. autoclass:: VacancyGenerator
   :members:
   :inherited-members:
   :show-inheritance:
   :undoc-members:

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

import itertools
import os

import numpy as np

from ..structure_io import DATAReader, DATAWriter, XYZWriter, \
    XYZ2DATAConverter, StructureFormatError, supported_structure_formats

__all__ = ['GrapheneVacancyGenerator',
           'NanotubeVacancyGenerator',
           'VacancyGenerator']


class VacancyGenerator(object):
    """Generate vacancies in structure data.

    Parameters
    ----------
    fname : str
        structure data filename
    Nvac : int
        total number of vacancies to "add" to structure data.
    structure_format : {None, str}, optional
        chemical file format of saved structure data. Must be one of:

            - xyz
            - data

        If ``None``, then guess based on ``fname`` file extension or
        default to ``xyz`` format.
    random_vacancies : bool, optional
        Generate random vacancies in structure data.
    uniform_vacancies : bool, optional
        Generate uniform vacancies in structure data.
    rotate_structure : bool, optional
        rotate structure data about specified ``rotation_axis`` by
        ``rotation_angle``.
    rotation_angle : float, optional
        Angle of rotation to rotate structure data in **degrees**.
    rotation_axis : {'x', 'y', 'z'}, optional
        if ``rotate_data is True``, ``rotation_axis`` must be set.

    """
    def __init__(self, fname=str, structure_format=None):

        if fname.endswith(supported_structure_formats) and \
                structure_format is None:
            for ext in supported_structure_formats:
                if fname.endswith(ext):
                    structure_format = ext
                    break
        else:
            if structure_format is None or \
                    structure_format not in supported_structure_formats:
                structure_format = 'xyz'

        self.fname = fname
        self.structure_format = structure_format

        # parse structure data
        self.atoms = None
        if self.structure_format == 'data':
            self.data = DATAReader(fname)
        elif self.structure_format == 'xyz':
            self.data = XYZ2DATAConverter(fname).convert(return_reader=True)
        else:
            raise StructureFormatError(
                '{} is not a supported structure format'.format(
                    structure_format))

        self.atoms = self.data.atoms
        self.atom_ids = self.atoms.atom_ids
        self.atom_coords = self.atoms.get_coords(as_dict=True)
        self.Nvac = 0
        self.vac_ids = None

    def _random_vacancy_generator(self):
        """Generate random vacancies in structure data."""
        self.vac_ids = \
            np.random.choice(self.atom_ids, size=self.Nvac, replace=False)

    def generate_vacancy_structure(self, show_vmd_selection_cmd=True):
        """Generate vacancy structure."""

        removed_atoms = self.atoms.filter_atoms(self.vac_ids, invert=False)

        if show_vmd_selection_cmd:
            vmd_vac_selection_cmd = []
            for atom in removed_atoms:
                selection_cmd = "(((x - {:.6f})^2 + ".format(atom.x) + \
                                "(y - {:.6f})^2 + ".format(atom.y) + \
                                "(z - {:.6f})^2) <= 9)".format(atom.z)
                vmd_vac_selection_cmd.append(selection_cmd)

            vmd_vac_selection_cmd = ' or '.join(vmd_vac_selection_cmd)
            print('copy and paste the following VMD command to select\n'
                  'the atoms surrounding the vacancies:\n\n'
                  '{}\n'.format(vmd_vac_selection_cmd))

        remaining_atoms = self.atoms.filter_atoms(self.vac_ids, invert=True)
        vacancy_structure_fname = \
            os.path.splitext(os.path.basename(self.fname))[0] + \
            '+{}_vacancies'.format(self.Nvac)
        DATAWriter.write(fname=vacancy_structure_fname,
                         atoms=remaining_atoms,
                         boxbounds=self.data.boxbounds,
                         comment_line=self.data.comment_line)
        XYZWriter.write(fname=vacancy_structure_fname,
                        atoms=remaining_atoms,
                        comment_line=self.data.comment_line)


class GrapheneVacancyGenerator(VacancyGenerator):
    """Generate vacancy structure.

    Parameters
    ----------
    fname : str
        structure data filename
    Nvac : int
        total number of vacancies to "add" to structure data.
    structure_format : str
        chemical file format of structure data. Must be one of:

            - xyz
            - data

    random_vacancies : bool, optional
        Generate random vacancies in structure data.
    uniform_vacancies : bool, optional
        Generate uniform vacancies in structure data.
    Nlayers : int, optional
        number of graphene layers.
    Nvac_per_layer : int, optional
        number of vacancies per layer.
    bin_dim : {'x', 'y', 'z'}, str
        axis along which to generate uniform vacancies
    nbin_sets : int
        number of bin sets to cycle through when picking the vacancy
        coordinates.
    crop_data : bool, optional
        Crop the structure data for each dimension in ``bounds`` by
        limits specified in ``bounds``.
    bounds : dict, optional
        if ``crop_data is True``, ``bounds`` must be provided as a dict
        of the form::

            {'dim': {'min': float, 'max': float}}

        for any dim in ``('x', 'y', 'z')''
    rotate_structure : bool, optional
        rotate structure data about specified ``rotation_axis`` by
        ``rotation_angle``.
    rotation_angle : float, optional
        Angle of rotation to rotate structure data in **degrees**.
    rotation_axis : {'x', 'y', 'z'}, optional
        if ``rotate_data is True``, ``rotation_axis`` must be set.

    Examples
    --------

    First, we'll generate a single layer sheet of graphene using
    the :py:class:`~sknano.nanogen.GrapheneGenerator` class and then
    "add" some vacancies to it using the
    :py:class:`~sknano.nanogen.GrapheneVacancyGenerator` class.

    >>> from sknano.nanogen import GrapheneGenerator, \
            GrapheneVacancyGenerator
    >>> ZZgraphene_1layer = GrapheneGenerator(width=5, length=10, edge='ZZ')
    >>> ZZgraphene_1layer.save_data(fname='5nmx10nm_ZZ_1layer.data',
    ...                             structure_format='data')

    Now we create an instance of the
    :py:class:`~sknano.nanogen.GrapheneVacancyGenerator` class
    and use the saved structure data output above as input to
    :py:class:`~sknano.nanogen.GrapheneVacancyGenerator`.

    >>> gvacgen = GrapheneVacancyGenerator(fname='5nmx10nm_ZZ_1layer.data',
    ...                                    structure_format='data')

    Now let's add 30 *random* vacancies to the graphene layer.

    .. code-block:: python

       >>> gvacgen.generate_vacancy_structure(Nvac=30, random=True)
       copy and paste the following VMD command to select
       the atoms surrounding the vacancies:

       (((x - -24.157000)^2 + (y - 0.000000)^2 + (z - 37.533974)^2) <= 9) or (((x - -24.157000)^2 + (y - 0.000000)^2 + (z - 22.766509)^2) <= 9) or (((x - -22.025500)^2 + (y - 0.000000)^2 + (z - 6.768422)^2) <= 9) or (((x - -22.025500)^2 + (y - 0.000000)^2 + (z - -7.999044)^2) <= 9) or (((x - -24.867500)^2 + (y - 0.000000)^2 + (z - -12.921532)^2) <= 9) or (((x - -17.762500)^2 + (y - 0.000000)^2 + (z - 46.148329)^2) <= 9) or (((x - -19.894000)^2 + (y - 0.000000)^2 + (z - 30.150241)^2) <= 9) or (((x - -18.473000)^2 + (y - 0.000000)^2 + (z - -43.687085)^2) <= 9) or (((x - -15.631000)^2 + (y - 0.000000)^2 + (z - 47.378951)^2) <= 9) or (((x - -16.341500)^2 + (y - 0.000000)^2 + (z - 16.613398)^2) <= 9) or (((x - -14.210000)^2 + (y - 0.000000)^2 + (z - 15.382776)^2) <= 9) or (((x - -14.210000)^2 + (y - 0.000000)^2 + (z - -28.919619)^2) <= 9) or (((x - -9.947000)^2 + (y - 0.000000)^2 + (z - 32.611486)^2) <= 9) or (((x - -9.236500)^2 + (y - 0.000000)^2 + (z - 26.458375)^2) <= 9) or (((x - -11.368000)^2 + (y - 0.000000)^2 + (z - 15.382776)^2) <= 9) or (((x - -12.078500)^2 + (y - 0.000000)^2 + (z - 4.307177)^2) <= 9) or (((x - -7.815500)^2 + (y - 0.000000)^2 + (z - 19.074643)^2) <= 9) or (((x - -4.973500)^2 + (y - 0.000000)^2 + (z - -47.378951)^2) <= 9) or (((x - -3.552500)^2 + (y - 0.000000)^2 + (z - 19.074643)^2) <= 9) or (((x - -1.421000)^2 + (y - 0.000000)^2 + (z - -16.613398)^2) <= 9) or (((x - -2.842000)^2 + (y - 0.000000)^2 + (z - -19.074643)^2) <= 9) or (((x - -0.710500)^2 + (y - 0.000000)^2 + (z - -39.995218)^2) <= 9) or (((x - 12.078500)^2 + (y - 0.000000)^2 + (z - 6.768422)^2) <= 9) or (((x - 15.631000)^2 + (y - 0.000000)^2 + (z - 25.227753)^2) <= 9) or (((x - 13.499500)^2 + (y - 0.000000)^2 + (z - 16.613398)^2) <= 9) or (((x - 15.631000)^2 + (y - 0.000000)^2 + (z - 12.921532)^2) <= 9) or (((x - 14.210000)^2 + (y - 0.000000)^2 + (z - -21.535887)^2) <= 9) or (((x - 16.341500)^2 + (y - 0.000000)^2 + (z - -35.072730)^2) <= 9) or (((x - 18.473000)^2 + (y - 0.000000)^2 + (z - 30.150241)^2) <= 9) or (((x - 19.894000)^2 + (y - 0.000000)^2 + (z - 12.921532)^2) <= 9)

    The rendered structure, with vacancies highlighted, looks like:

    .. image:: /image/5nmx10nm_ZZ_1layer+30_vacancies.png

    """
    def __init__(self, fname=str, structure_format='xyz',
                 rotate_structure=False, rotation_angle=None,
                 rotation_axis=None):

        super(GrapheneVacancyGenerator, self).__init__(
            fname=fname, structure_format=structure_format)

    def generate_vacancy_structure(self, Nvac=int, random=False, uniform=False,
                                   Nlayers=None, Nvac_per_layer=None,
                                   bin_dim=None, nbin_sets=None,
                                   show_vmd_selection_cmd=True):

        self.Nvac = Nvac

        if random:
            super(GrapheneVacancyGenerator, self)._random_vacancy_generator()
        else:
            self.vac_ids = np.empty(0, dtype=int)

            # find the coords of each layer
            y_coords = self.atom_coords['y']
            layer_coords = np.asarray(sorted(list(set(y_coords))))
            bin_sets = []
            nbins = int(Nvac / Nvac_per_layer)
            if nbin_sets is None:
                bin_sets.append(np.arange(0, nbins, 1))
            else:
                for n in xrange(nbin_sets):
                    bin_sets.append(np.arange(n, nbins, nbin_sets))

            bin_set_iter = itertools.cycle((bin_sets))

            vac_bin_edges = np.linspace(self.atom_coords[bin_dim].min(),
                                        self.atom_coords[bin_dim].max(),
                                        num=nbins+1)
            vac_coords_along_bin_dim = \
                vac_bin_edges[:-1] + np.diff(vac_bin_edges) / 2

            for layer_pos in layer_coords:
                bin_set = bin_set_iter.next()
                for vac_pos in vac_coords_along_bin_dim[bin_set]:
                    candidate_vac_atom_indices = \
                        np.where(
                            (self.atom_coords['x'] >=
                             (self.atom_coords['x'].min() + 2.5)) &
                            (self.atom_coords['x'] <=
                             (self.atom_coords['x'].max() - 2.5)) &
                            (np.abs(self.atom_coords['y'] - layer_pos)
                                <= 0.5) &
                            (np.abs(self.atom_coords[bin_dim] - vac_pos)
                                <= 0.75))
                    candidate_vac_atom_ids = \
                        self.atom_ids[candidate_vac_atom_indices]
                    print('candidate_vac_atom_ids: '
                          '{}\n'.format(candidate_vac_atom_ids))

                    self.vac_ids = \
                        np.r_[self.vac_ids, np.random.choice(
                            candidate_vac_atom_ids)]

        super(GrapheneVacancyGenerator, self).generate_vacancy_structure(
            show_vmd_selection_cmd=show_vmd_selection_cmd)


class NanotubeVacancyGenerator(VacancyGenerator):

    def __init__(self, fname=str, structure_format='xyz',
                 rotate_structure=False, rotation_angle=None,
                 rotation_axis=None):

        super(NanotubeVacancyGenerator, self).__init__(
            fname=fname, structure_format=structure_format)

    def generate_vacancy_structure(self, Nvac=int, random=False, uniform=False,
                                   Ntubes=None, Nvac_per_tube=None,
                                   bin_dim=None, nbin_sets=None,
                                   show_vmd_selection_cmd=True):

        self.Nvac = Nvac

        if random:
            super(NanotubeVacancyGenerator, self)._random_vacancy_generator()
        else:
            self.vac_ids = np.empty(0, dtype=int)

            Natoms = len(self.atom_ids)
            print('Natoms: {}'.format(Natoms))
            Natoms_per_tube = int(Natoms / Ntubes)
            print('Natoms_per_tube: {}\n'.format(Natoms_per_tube))

            bin_sets = []
            nbins = Nvac / Nvac_per_tube
            if nbin_sets is None:
                bin_sets.append(np.arange(0, nbins, 1))
            else:
                #tube_array = np.arange(Ntubes)
                #np.random.shuffle(tube_array)
                for n in xrange(nbin_sets):
                    bin_sets.append(np.arange(n, nbins, nbin_sets))

            bin_set_iter = itertools.cycle((bin_sets))

            for i in xrange(Ntubes):
                print('tube n: {:d}'.format(i+1))
                tube_atom_indices = \
                    np.where((self.atom_ids > (Natoms_per_tube * i)) &
                             (self.atom_ids <= (Natoms_per_tube * (i + 1))))
                tube_atom_ids = self.atom_ids[tube_atom_indices]
                print('tube_atom_ids:\n{}\n'.format(tube_atom_ids))

                tube_coords = self.atoms.get_filtered_coords(tube_atom_ids,
                                                             as_dict=True,
                                                             invert=False)

                vac_bin_edges = np.linspace(tube_coords[bin_dim].min(),
                                            tube_coords[bin_dim].max(),
                                            num=nbins+1)
                vac_coords_along_bin_dim = \
                    vac_bin_edges[:-1] + np.diff(vac_bin_edges) / 2

                #if nbins < Nvac_per_tube
                #vac_coords_along_bin_dim = \
                bin_set = bin_set_iter.next()

                for vac_pos in vac_coords_along_bin_dim[bin_set]:
                    candidate_vac_atom_indices = \
                        np.where(np.abs(tube_coords[bin_dim] - vac_pos) <= 1)
                    candidate_vac_atom_ids = \
                        tube_atom_ids[candidate_vac_atom_indices]

                    print(u'vac_pos along {}-axis: {:.2f} \u00c5'.format(
                        bin_dim, vac_pos))
                    print('N candidate_vac_atom_ids: '
                          '{}\n'.format(len(candidate_vac_atom_ids)))

                    self.vac_ids = \
                        np.r_[self.vac_ids,
                              np.random.choice(candidate_vac_atom_ids)]
