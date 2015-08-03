# -*- coding: utf-8 -*-
"""
=================================================================================
Nanotube defect generators (:mod:`sknano.generators._nanotube_defect_generators`)
=================================================================================

.. currentmodule:: sknano.generators._nanotube_defect_generators

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import itertools

import numpy as np

from sknano.core import xyz
from sknano.core.refdata import CCbond

from ._defect_generators import CrossLinkedDefectGenerator, \
    StoneWalesDefectGenerator, VacancyGenerator
from ._nanotube_bundle_generators import NanotubeBundleGenerator

__all__ = ['NanotubeStoneWalesDefectGenerator',
           'NanotubeVacancyGenerator',
           'CrossLinkedNanotubeBundleGenerator']


class NanotubeStoneWalesDefectGenerator(StoneWalesDefectGenerator):
    """Generate Stone-Wales defects in nanotube structure data.

    .. versionadded:: 0.2.26

    Parameters
    ----------
    n, m : int, optional
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions.
    fname : str, optional
        Structure data filename. If you don't provide a structure data file,
        the you **must** provide the structure data parameters to
        generate the structure data file.
    structure_format : {None, str}, optional
        Chemical file format of saved structure data.
        If `None`, then guess based on `fname` file extension.
        Otherwise, must be one of:

            - xyz
            - data

    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    bundle_packing : {None, 'hexagonal', 'cubic'}, optional
        close packing arrangement of bundles
    bundle_geometry : {None, 'triangle', 'hexagon', 'square', 'rectangle',
                       'rhombus', 'rhomboid'}, optional
    Lx, Ly, Lz : float, optional
        length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.
    verbose : bool, optional
        Verbose output

    Raises
    ------
    `TypeError`
        If `fname` is `None` and `n` and `m` are not integers.

    Examples
    --------

    You can provide any existing structure data file in any supported
    format or else provide structure parameters to generate a new structure
    in real-time. The `NanotubeStoneWalesDefectGenerator` uses
    the `NanotubeBundleGenerator` class to generate
    nanotube structure data. See the structure generator classes provided
    by the :mod:`~sknano.generators` module for more detailed docs.

    In the next example, we'll generate a nanotube bundle and add
    some SW defects to it.

    >>> from sknano.generators import NanotubeStoneWalesDefectGenerator
    >>> swntbundle = NanotubeStoneWalesDefectGenerator(
    ...     n=10, m=5, Lz=10, fix_Lz=True, bundle_geometry='hexagon')

    Notice that I used the `Lz` keyword argument to specify the
    **length** of the nanotube bundle in **nanometers**, and I also
    set `fix_Lz=True`. If you don't set `fix_Lz=True`, then the
    length of the nanotube will get truncated such that it's an integer
    multiple of the unit cell length, which, depending on the chirality, may
    end up being very different than the desired length set by `Lz`.

    Next, I add 20 SW defects distributed uniformly along the
    :math:`z`-axis. For the `NanotubeStoneWalesDefectGenerator`,
    if `uniform` is `True`, the code is written such it will divide the total
    number of defects evenly among all nanotubes and distribute those
    defects uniformly along the specified `bin_axis`. The circumferential
    location of defects is selected randomly.  Future revisions will allow
    more precise control of the positioning.

    >>> swntbundle.generate_structure_defects(Nvac_sites=20, uniform=True,
    ...                                       bin_axis='z', defect_type='575')

    Running that command left me with the saved defect structure data which
    I loaded into VMD. I used the VMD selection command that was printed to the
    terminal window to highlight the atoms and bonds surrounding the defects.

    The rendered structure looks like:

    After relaxing the defected structure data in LAMMPS, the rendered
    structure looks like:

    """
    def __init__(self, n=None, m=None, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond,
                 bundle_packing=None, bundle_geometry=None,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False,
                 rotate_structure=False, rotation_angle=None,
                 rotation_axis=None,
                 fname=None, outpath=None, structure_format=None,
                 verbose=False, **kwargs):

        self._Ntubes = None
        self._Natoms_per_tube = None

        if fname is None and isinstance(n, int) and isinstance(m, int):
            generator = NanotubeBundleGenerator(
                n=n, m=m, nx=nx, ny=ny, nz=nz, element1=element1,
                element2=element2, bond=bond,
                bundle_packing=bundle_packing, bundle_geometry=bundle_geometry,
                Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz, verbose=verbose)
            generator.save(structure_format='data', **kwargs)
            fname = generator.fname
            structure_format = generator.structure_format
            self._Ntubes = generator.Ntubes
            self._Natoms_per_tube = generator.Natoms_per_tube
        else:
            raise TypeError('Either `fname` but be provided or '
                            '`n` and `m` must be specified as integers')

        super(NanotubeStoneWalesDefectGenerator, self).__init__(
            fname=fname, outpath=outpath, structure_format=structure_format,
            verbose=verbose)

    def generate_structure_defects(self, N_defect_sites=None, uniform=False,
                                   bin_axis='z', defect_bounds=None,
                                   Ntubes=None,
                                   vmd_selection_radius=np.sqrt(10.5),
                                   show_vmd_selection_cmd=True):
        """Generate structure defects.

        Parameters
        ----------
        uniform : bool, optional
            Generate uniform defects in structure data.
        bin_axis : {'x', 'y', 'z'}, optional
            axis along which to generate uniform defects
        Ntubes : {None, int}, optional
        show_vmd_selection_cmd : bool, optional
            Generate a VMD selection string that can be used to
            select the atoms surrounding the defects.
        vmd_selection_radius : float, optional
            Cutoff radius for VMD selection command in units of **Angstroms**.

        """

        self._vmd_selection_radius = vmd_selection_radius
        self._show_vmd_selection_cmd = show_vmd_selection_cmd

        if uniform:
            pass
        else:
            super(NanotubeStoneWalesDefectGenerator,
                  self)._generate_random_defects()

        super(NanotubeStoneWalesDefectGenerator,
              self)._generate_structure_defects()


class CrossLinkedNanotubeBundleGenerator(CrossLinkedDefectGenerator):
    pass


class NanotubeVacancyGenerator(VacancyGenerator):
    """Generate vacancies in nanotubes.

    .. versionadded:: 0.2.6

    Parameters
    ----------
    n, m : int, optional
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions.
    fname : str, optional
        Structure data filename. If you don't provide a structure data file,
        the you **must** provide the structure data parameters to
        generate the structure data file.
    structure_format : {None, str}, optional
        Chemical file format of saved structure data.
        If `None`, then guess based on `fname` file extension.
        Otherwise, must be one of:

            - xyz
            - data

    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    bundle_packing : {None, 'hexagonal', 'cubic'}, optional
        close packing arrangement of bundles
    bundle_geometry : {None, 'triangle', 'hexagon', 'square', 'rectangle',
                       'rhombus', 'rhomboid'}, optional
    Lx, Ly, Lz : float, optional
        length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.
    verbose : bool, optional
        Verbose output

    Raises
    ------
    `TypeError`
        If `fname` is `None` and `n` and `m` are not integers.

    Examples
    --------

    The `NanotubeVacancyGenerator` class behaves nearly
    identically to the `GrapheneVacancyGenerator` class.
    You can provide any existing structure data file in any supported
    format or else provide structure parameters to generate a structure
    in real-time. The `NanotubeVacancyGenerator` uses
    the `NanotubeBundleGenerator` class to generate
    nanotubes. See the structure generator classes provided by the
    :mod:`~sknano.generators` module for more detailed docs.

    In the next example, we'll generate a nanotube bundle and then
    poke some holes in it.

    >>> from sknano.generators import NanotubeVacancyGenerator
    >>> ntvg = NanotubeVacancyGenerator(n=10, m=5, Lz=10, fix_Lz=True,
    ...                                 bundle_geometry='hexagon')

    Notice that I used the `Lz` keyword argument to specify the
    **length** of the nanotube bundle in **nanometers**, and I also
    set `fix_Lz=True`. If you don't set `fix_Lz=True`, then the
    length of the nanotube will get truncated such that it's an integer
    multiple of the unit cell length, which,
    depending on the chirality, may end up being very different than
    the desired length set by `Lz`.

    Next, I add 35 **double** vacancies distributed uniformly along the
    :math:`z`-axis. For the `NanotubeVacancyGenerator`, if `uniform` is `True`,
    the code is written such it will divide the total number of vacancies
    evenly among all nanotubes and distribute those vacancies uniformly along
    the specified `bin_axis`. The circumferential location of vacancies is
    selected randomly.  Future revisions will allow more precise control of the
    positioning.

    >>> ntvg.generate_vacancy_structure(Nvac_sites=35, uniform=True,
    ...                                 bin_axis='z', vac_type='double')

    Running that command left me with the saved vacancy structure data which
    I loaded into VMD. I used the VMD selection command that was printed to the
    terminal window to highlight the atoms and bonds surrounding the vacancies.

    The rendered structure looks like:

    .. image:: /images/1005_hcp_7tube_hexagon+70_vacancies-01.png

    Here's a final example which generates a *tri-vacancy* at each
    `Nvac_sites`.

    >>> ntvg = NanotubeVacancyGenerator(n=10, m=10, nz=50)
    >>> ntvg.generate_vacancy_structure(Nvac_sites=20, uniform=True,
    ...                                 bin_axis='z', vac_type='triple')

    After relaxing the vacancy structure data in LAMMPS, the rendered
    structure looks like:

    .. image:: /images/1010_hcp_1cellx1cellx50cells+60_vacancies-01.png

    This poor nanotube has taken quite a beating.

    """
    def __init__(self, n=None, m=None, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond,
                 bundle_packing=None, bundle_geometry=None,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False,
                 rotate_structure=False, rotation_angle=None,
                 rotation_axis=None,
                 fname=None, outpath=None, structure_format=None,
                 verbose=False, **kwargs):

        self._Ntubes = None
        self._Natoms_per_tube = None

        if fname is None and isinstance(n, int) and isinstance(m, int):
            ntbg = NanotubeBundleGenerator(n=n, m=m, nx=nx, ny=ny, nz=nz,
                                           element1=element1,
                                           element2=element2,
                                           bond=bond,
                                           bundle_packing=bundle_packing,
                                           bundle_geometry=bundle_geometry,
                                           Lx=Lx, Ly=Ly, Lz=Lz,
                                           fix_Lz=fix_Lz,
                                           verbose=verbose)
            ntbg.save(structure_format='data', **kwargs)
            fname = ntbg.fname
            structure_format = ntbg.structure_format
            self._Ntubes = ntbg.Ntubes
            self._Natoms_per_tube = ntbg.Natoms_per_tube
        else:
            raise TypeError('Either `fname` but be provided or '
                            '`n` and `m` must be specified as integers')

        super(NanotubeVacancyGenerator, self).__init__(
            fname=fname, outpath=outpath, structure_format=structure_format,
            verbose=verbose)

    def generate_vacancy_structure(self, Nvac_sites=None, Nvac_clusters=None,
                                   cluster_size=None, vac_type='single',
                                   uniform=False, bin_axis='z',
                                   vac_bounds=None, Ntubes=None,
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
            Generate uniform vacancies in structure data.
        bin_axis : {'x', 'y', 'z'}, optional
            axis along which to generate uniform vacancies
        Ntubes : {None, int}, optional
        show_vmd_selection_cmd : bool, optional
            Generate a VMD selection string that can be used to
            select the atoms surrounding the vacancies.
        vmd_selection_radius : float, optional
            Cutoff radius for VMD selection command in units of **Angstroms**.

        Raises
        ------
        TypeError
            If `Nvac_sites` is `None` and `Nvac_clusters` is `None` or
            if `uniform` is `True` and the `structure_data` was read in
            from a file and `Ntubes` is `None`.

        """

        if Nvac_sites is None and Nvac_clusters is None:
            raise TypeError('`Nvac_sites` or `Nvac_clusters` must be an '
                            'integer.')
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
            if Ntubes is None:
                Ntubes = self._Ntubes
            if not isinstance(Ntubes, int):
                raise TypeError('`Ntubes` must be specified as an integer.')

            Natoms = len(self._atom_ids)
            if self._Natoms_per_tube is not None:
                Natoms_per_tube = self._Natoms_per_tube
            else:
                Natoms_per_tube = int(Natoms / Ntubes)

            Nvac_sites_per_tube = int(Nvac_sites / Ntubes)
            extra = Nvac_sites % Ntubes
            nbins = Ntubes * Nvac_sites_per_tube

            bin_sets = []
            for n in range(Ntubes):
                bin_sets.append(np.arange(n, Nvac_sites, Ntubes))
            bin_set_iter = itertools.cycle((bin_sets))

            vac_bin_edges = np.linspace(self._atom_coords.T[bin_axis].min(),
                                        self._atom_coords.T[bin_axis].max(),
                                        num=nbins+1)
            vac_coords_along_bin_axis = \
                vac_bin_edges[:-1] + np.diff(vac_bin_edges) / 2

            if self._verbose:
                print('Ntubes: {}'.format(Ntubes))
                print('Natoms: {}'.format(Natoms))
                print('Natoms_per_tube: {}'.format(Natoms_per_tube))
                print('Nvac_sites: {}'.format(Nvac_sites))
                print('Nvac_sites_per_tube: {}'.format(Nvac_sites_per_tube))
                print('extra vacancies: {}'.format(extra))
                print('nbins: {}'.format(nbins))
                print('bin_sets:\n{}'.format(bin_sets))
                print('vac_bin_edges:'
                      '\n{}'.format(vac_bin_edges))
                print('vac_coords_along_bin_axis:'
                      '\n{}'.format(vac_coords_along_bin_axis))

            for i in range(Ntubes):
                tube_atom_indices = \
                    np.where((self._atom_ids > (Natoms_per_tube * i)) &
                             (self._atom_ids <= (Natoms_per_tube * (i + 1))))

                tube_atom_ids = self._atom_ids[tube_atom_indices]

                tube_coords = self._atoms.filter_ids(tube_atom_ids,
                                                     invert=False).coords

                if self._verbose:
                    print('tube n: {:d}'.format(i + 1))
                    print('tube_atom_ids:\n{}'.format(tube_atom_ids))
                    print('tube_atom_coords:\n{}'.format(tube_coords))

                bin_set = next(bin_set_iter)
                for vac_pos in vac_coords_along_bin_axis[bin_set]:
                    candidate_vac_atom_indices = \
                        np.where(np.abs(tube_coords.T[bin_axis_index] -
                                        vac_pos) <= 1)
                    candidate_vac_atom_ids = \
                        tube_atom_ids[candidate_vac_atom_indices]

                    if self._verbose:
                        print('vac_pos along {}-axis: {:.2f} \u00c5'.format(
                            bin_axis, vac_pos))
                        print('N candidate_vac_atom_ids: '
                              '{}\n'.format(len(candidate_vac_atom_ids)))

                    self._vac_ids = \
                        np.r_[self._vac_ids,
                              np.random.choice(candidate_vac_atom_ids)]
        else:
            super(NanotubeVacancyGenerator, self)._random_vacancy_generator()

        super(NanotubeVacancyGenerator, self)._generate_vacancy_structure()
