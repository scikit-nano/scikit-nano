# -*- coding: utf-8 -*-
"""
=================================================================================
Graphene defect generators (:mod:`sknano.generators._graphene_defect_generators`)
=================================================================================

.. currentmodule:: sknano.generators._graphene_defect_generators

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import itertools

import numpy as np

from sknano.core import xyz
from sknano.core.refdata import CCbond

from ._defect_generators import DefectGenerator, VacancyGenerator
from ._graphene_generator import GrapheneGenerator
from ._unrolled_swnt_generator import UnrolledSWNTGenerator

__all__ = ['GrapheneDefectGenerator', 'GrapheneVacancyGenerator']


class GrapheneDefectGenerator(DefectGenerator):
    """Graphene defect generator class."""
    pass


class GrapheneVacancyGenerator(VacancyGenerator):
    """Generate vacancies in graphene.

    .. versionadded:: 0.2.6

    Parameters
    ----------
    length : float, optional
        Length of graphene sheet in **nanometers**
    width : float, optional
        Width of graphene sheet in **nanometers**
    edge : {'AC', 'armchair', 'ZZ', 'zigzag'}, optional
        **A**\ rm\ **C**\ hair or **Z**\ ig\ **Z**\ ag edge along
        the `length` of the sheet.
    fname : str, optional
        Structure data filename. If you don't provide a structure data file,
        the you **must** provide the structure data parameters to
        generate the structure data file.
    outpath : str, optional
        Output path for structure data file.
    structure_format : {None, str}, optional
        Chemical file format of saved structure data.
        If `None`, then guess based on `fname` file extension.
        Otherwise, must be one of:

            - `xyz`
            - `data`

    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2.
    bond : float, optional
        bond length between nearest-neighbor atoms in **Angstroms**.
    nlayers : int, optional
        Number of graphene layers.
    layer_spacing : float, optional
        Distance between layers in **Angstroms**.
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers
    verbose : bool, optional
        Verbose output

    Raises
    ------
    `TypeError`
        If `fname` is `None` and `n` and `m` are not integers or
        `length` and `width` are not floats.

    Examples
    --------

    Import the `GrapheneVacancyGenerator` class.

    >>> from sknano.generators import GrapheneVacancyGenerator

    You can supply any existing structure data file
    (as long as its a supported format) to
    `GrapheneVacancyGenerator` using the
    `fname` keyword argument. For example, if you have an existing
    LAMMPS structure data file for graphene named
    *10nmx5nm_ZZ_1layer.data*, you could load that into your
    vacancy generator by typing:

    >>> gvacgen = GrapheneVacancyGenerator(fname='10nmx5nm_ZZ_1layer.data')

    If say, you had your LAMMPS data file saved with a `.lammps`
    extension, then you would need to specify the `structure_format`
    keyword argument as well:

    >>> gvacgen = GrapheneVacancyGenerator(fname='10nmx5nm_ZZ_1layer.lammps',
    ...                                    structure_format='data')

    If you don't have an existing structure data file and want to generate
    your graphene structure on-the-fly, then you need to provide
    structure parameters to the `GrapheneVacancyGenerator`
    constructor. These arguments are passed to an instance of
    `GrapheneGenerator` "under the hood" (see the `GrapheneGenerator` docs
    for examples and more detailed documentation). The `length` and `width`
    parameters are required.

    In the following example, we'll use the `GrapheneVacancyGenerator` class
    to generate a *10nm x 2nm*, single layer graphene with `armchair` edges,
    and then poke some holes in it.

    >>> gvacgen = GrapheneVacancyGenerator(length=10, width=2, edge='AC')

    Now let's add 20 *random* vacancies to the graphene layer by
    calling the `~GrapheneVacancyGenerator.generate_vacancy_structure`
    method:

    .. code-block:: python

       >>> gvacgen.generate_vacancy_structure(Nvac_sites=20)
       copy and paste the following VMD command to select
       the atoms surrounding the vacancies:

       (((x - -38.7646)^2 + (y - 0.0000)^2 + (z - 2.8420)^2) <= 12.00) or
       (((x - -38.7646)^2 + (y - 0.0000)^2 + (z - -5.6840)^2) <= 12.00) or
       (((x - -36.3034)^2 + (y - 0.0000)^2 + (z - -2.8420)^2) <= 12.00) or
       (((x - -36.3034)^2 + (y - 0.0000)^2 + (z - -7.1050)^2) <= 12.00) or
       (((x - -10.4603)^2 + (y - 0.0000)^2 + (z - 7.8155)^2) <= 12.00) or
       (((x - -6.7684)^2 + (y - 0.0000)^2 + (z - -2.8420)^2) <= 12.00) or
       (((x - -4.3072)^2 + (y - 0.0000)^2 + (z - -7.1050)^2) <= 12.00) or
       (((x - -3.0766)^2 + (y - 0.0000)^2 + (z - -3.5525)^2) <= 12.00) or
       (((x - 6.7684)^2 + (y - 0.0000)^2 + (z - -4.9735)^2) <= 12.00) or
       (((x - 6.7684)^2 + (y - 0.0000)^2 + (z - -9.2365)^2) <= 12.00) or
       (((x - 19.0746)^2 + (y - 0.0000)^2 + (z - 4.9735)^2) <= 12.00) or
       (((x - 22.7665)^2 + (y - 0.0000)^2 + (z - -9.9470)^2) <= 12.00) or
       (((x - 25.2278)^2 + (y - 0.0000)^2 + (z - -1.4210)^2) <= 12.00) or
       (((x - 27.6890)^2 + (y - 0.0000)^2 + (z - -7.1050)^2) <= 12.00) or
       (((x - 32.6115)^2 + (y - 0.0000)^2 + (z - 5.6840)^2) <= 12.00) or
       (((x - 33.8421)^2 + (y - 0.0000)^2 + (z - 0.7105)^2) <= 12.00) or
       (((x - 35.0727)^2 + (y - 0.0000)^2 + (z - -1.4210)^2) <= 12.00) or
       (((x - 37.5340)^2 + (y - 0.0000)^2 + (z - 5.6840)^2) <= 12.00) or
       (((x - 38.7646)^2 + (y - 0.0000)^2 + (z - -9.2365)^2) <= 12.00) or
       (((x - 39.9952)^2 + (y - 0.0000)^2 + (z - -9.9470)^2) <= 12.00)

    This will save the vacancy structure data in both
    `LAMMPS data` and `xyz` formats. You can load the `xyz` file
    into VMD and then copy the VMD selection command from the output above
    and paste it into VMD to highlight the atoms surrounding the vacancies.

    As you know, removing atoms from a material will change the equilibrium
    atomic arrangement of the atoms. Obtaining the new equilibrium structure
    in the vacancy structure data requires running molecular dynamics
    to compute the new forces on the atoms and move them to a minimum energy
    configuration. Currently this is not supported within this toolkit and
    requires using an external molecular dynamics package.  This toolkit
    outputs the structure data in the native `LAMMPS data` file format and is
    the molecular dynamics engine I used to relax the structure data after
    adding the vacancies.

    The rendered, relaxed structure, with vacancies highlighted, looks like:

    .. image:: /images/10nmx2nm_1layer_AC_CC_graphene+20_vacancies-01.png

    In addition adding single-vacancies, you can generate di-vacancies
    and tri-vacancies by setting `vac_type` to `double` and `triple`,
    respectively.

    Here's an example of generating single layer graphene with a
    `di-vacancy` at each of the 20 `Nvac_sites`.

    >>> gvacgen = GrapheneVacancyGenerator(length=10, width=5, edge='ZZ')
    >>> gvacgen.generate_vacancy_structure(Nvac_sites=20, vac_type='double')

    I'm not going to show the VMD selection command since it will not
    be the same if you run this example on your computer. After generating
    the vacancy structure, I relaxed the structure data using LAMMPS.

    Here's the rendered structure, with the 40 vacancies highlighted:

    .. image:: /images/10nmx5nm_1layer_ZZ_CC_graphene+40_vacancies-01.png

    """
    def __init__(self, n=None, m=None, nx=1, ny=1, nz=1,
                 length=None, width=None, edge=None,
                 element1='C', element2='C', bond=CCbond,
                 nlayers=1, layer_spacing=3.35, stacking_order='AB',
                 fname=None, outpath=None, structure_format=None,
                 rotate_structure=False, rotation_angle=None,
                 rotation_axis=None, verbose=False, **kwargs):

        if fname is None:
            if isinstance(length, (int, float)) and \
                    isinstance(width, (int, float)):
                gg = GrapheneGenerator(
                    length=length, width=width, edge=edge, element1=element1,
                    element2=element2, bond=bond, nlayers=nlayers,
                    layer_spacing=layer_spacing, stacking_order=stacking_order,
                    verbose=verbose)
            elif isinstance(n, int) and isinstance(m, int):
                gg = UnrolledSWNTGenerator(
                    n=n, m=m, nx=nx, ny=ny, nz=nz,
                    element1=element1, element2=element2, bond=bond,
                    nlayers=nlayers, layer_spacing=layer_spacing,
                    stacking_order=stacking_order, verbose=verbose)
            else:
                raise TypeError('Either `fname` must be set or `n` and `m` '
                                'must be specified as integers or `length` '
                                'and `width` must be specified as floats.')
            gg.save(outpath=outpath, structure_format='data', **kwargs)
            fname = gg.fname
            structure_format = gg.structure_format

        super(GrapheneVacancyGenerator, self).__init__(
            fname=fname, outpath=outpath,
            structure_format=structure_format, verbose=verbose)

    def generate_vacancy_structure(self, Nvac_sites=None, Nvac_clusters=None,
                                   cluster_size=None, vac_type='single',
                                   uniform=False, bin_axis='z',
                                   distribute_evenly=False,
                                   vmd_selection_radius=np.sqrt(10.5),
                                   show_vmd_selection_cmd=True):
        """Generate vacancy structure.

        Parameters
        ----------
        Nvac_sites : int, optional
            total number of vacancy sites to "add" to structure data.
        Nvac_clusters : int, optional
            total number of vacancy cluster sites to "add" to structure data.
        vac_type : {'single', 'double', 'triple'}, optional
        uniform : bool, optional
            Generate vacancies uniformly distributed along `bin_axis`.
        bin_axis : {'x', 'y', 'z'}, optional
            axis along which to generate uniform vacancies
        distribute_evenly : bool, optional
        show_vmd_selection_cmd : bool, optional
            Generate a VMD selection string that can be used to
            select the atoms surrounding the vacancies.
        vmd_selection_radius : float, optional
            Cutoff radius for VMD selection command in units of **Angstroms**.

        Raises
        ------
        `TypeError`
            If `Nvac_sites` is `None` and `Nvac_clusters` is `None`

        """
        if Nvac_sites is None and Nvac_clusters is None:
            raise TypeError('`Nvac_sites` or `Nvac_clusters` must '
                            'be an integer.')
        elif Nvac_sites is None and Nvac_clusters is not None:
            Nvac_sites = Nvac_clusters
            if cluster_size is None:
                cluster_size = 1

        self._Nvac_sites = Nvac_sites
        self._Nvac_clusters = Nvac_clusters
        self._cluster_size = cluster_size

        self._vac_type = vac_type
        self._vmd_selection_radius = vmd_selection_radius
        self._show_vmd_selection_cmd = show_vmd_selection_cmd
        bin_axis_index = xyz.index(bin_axis)

        if uniform:
            # find the coords of each layer
            y_coords = self._atom_coords['y']
            Nlayers = len(set(y_coords))
            layer_coords = np.asarray(sorted(list(set(y_coords))))

            Nvac_sites_per_layer = int(Nvac_sites / Nlayers)
            extra = Nvac_sites % Nlayers

            nbins = Nvac_sites_per_layer * Nlayers

            bin_sets = []
            for n in range(Nlayers):
                bin_sets.append(np.arange(n, nbins, Nlayers))
            bin_set_iter = itertools.cycle((bin_sets))

            bin_edges = np.linspace(self._atom_coords.T[bin_axis_index].min(),
                                    self._atom_coords.T[bin_axis_index].max(),
                                    num=nbins+1)

            bin_mid_pts = bin_edges[:-1] + np.diff(bin_edges) / 2

            if self._verbose:
                print('Nlayers: {}'.format(Nlayers))
                print('Nvac_sites: {}'.format(Nvac_sites))
                print('Nvac_sites_per_layer: {}'.format(Nvac_sites_per_layer))
                print('extra vacancies: {}'.format(extra))
                print('nbins: {}'.format(nbins))
                print('bin_sets:\n{}'.format(bin_sets))
                print('bin_edges:\n{}'.format(bin_edges))
                print('bin_mid_pts:\n{}'.format(bin_mid_pts))

            ortho_axis = 'x' if bin_axis == 'z' else 'z'
            ortho_mid_pt = self._atom_coords[ortho_axis].min() + \
                (self._atom_coords[ortho_axis].max() -
                 self._atom_coords[ortho_axis].min()) / 2
            atoms_along_ortho_axis = {'+': [], '-': []}

            for y in layer_coords:
                bin_set = next(bin_set_iter)
                for vac_pos in bin_mid_pts[bin_set]:
                    candidate_vac_atom_indices = \
                        np.where(
                            (self._atom_coords[ortho_axis] >=
                             (self._atom_coords[ortho_axis].min() + 2.5)) &
                            (self._atom_coords[ortho_axis] <=
                             (self._atom_coords[ortho_axis].max() - 2.5)) &
                            (np.abs(self._atom_coords['y'] - y) <= 0.5) &
                            (np.abs(self._atom_coords.T[bin_axis_index] -
                                    vac_pos)
                                <= 1))
                    candidate_vac_atom_ids = \
                        self._atom_ids[candidate_vac_atom_indices]

                    if self._verbose:
                        print('candidate_vac_atom_ids: '
                              '{}\n'.format(candidate_vac_atom_ids))

                    rand_vac_atom_id = None
                    while True:
                        rand_vac_atom_id = \
                            np.random.choice(candidate_vac_atom_ids)
                        if distribute_evenly:
                            rand_vac_atom = \
                                self._atoms.get_atoms(asarray=True)[
                                    self._atom_ids == rand_vac_atom_id][0]
                            ortho_pos = getattr(rand_vac_atom, ortho_axis)
                            if ortho_pos >= ortho_mid_pt and \
                                    len(atoms_along_ortho_axis['+']) < \
                                    Nvac_sites_per_layer / 2:
                                atoms_along_ortho_axis['+'].append(
                                    rand_vac_atom)
                                break
                            elif ortho_pos < ortho_mid_pt and \
                                    len(atoms_along_ortho_axis['-']) < \
                                    Nvac_sites_per_layer / 2:
                                atoms_along_ortho_axis['-'].append(
                                    rand_vac_atom)
                                break
                            else:
                                continue
                        else:
                            break

                    self._vac_ids = np.r_[self._vac_ids, rand_vac_atom_id]
        else:
            super(GrapheneVacancyGenerator, self)._random_vacancy_generator()

        super(GrapheneVacancyGenerator, self)._generate_vacancy_structure()
