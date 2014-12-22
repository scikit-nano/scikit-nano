# -*- coding: utf-8 -*-
"""
==============================================================================
SWNT structure class (:mod:`sknano.structures._swnt`)
==============================================================================

.. currentmodule:: sknano.structures._swnt

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from ._base import StructureBase
from ._compute_funcs import compute_tube_mass, compute_linear_mass_density
from ._extras import attr_strfmt, attr_symbols, attr_units
from ._mixins import SWNTMixin

__all__ = ['SWNT', 'Nanotube']


class SWNT(SWNTMixin, StructureBase):
    """SWNT structure class.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nz : int, optional
        Number of repeat unit cells along the :math:`z` dimension
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b
    Lz : float, optional
        Length of the nanotube in **nanometers**.
        Overrides the :math:`n_z` cell values.

        .. versionadded:: 0.2.5

    tube_length : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides `nz` value.

        .. deprecated:: 0.2.5
           Use `Lz` instead

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.

    Examples
    --------

    >>> from sknano.structures import SWNT

    Create a SWNT with :math:`\\mathbf{C}_{h} = (10, 10)` chirality.

    >>> swnt = SWNT(n=10, m=10, verbose=True)
    >>> print(unicode(swnt))
    SWNT(n=10, m=10, element1='C', element2='C', bond=1.42, nz=1)
    n: 10
    m: 10
    t₁: 1
    t₂: -1
    d: 10
    dR: 30
    N: 20
    R: (1, 0)
    θc: 30.00°
    Ch: 42.60 Å
    T: 2.46 Å
    dt: 13.56 Å
    rt: 6.78 Å
    electronic_type: metallic

    Change the chirality to :math:`\\mathbf{C}_{h} = (20, 10)`.

    >>> swnt.n = 20
    >>> print(unicode(swnt))
    SWNT(n=20, m=10, element1='C', element2='C', bond=1.42, nz=1)
    n: 20
    m: 10
    t₁: 4
    t₂: -5
    d: 10
    dR: 10
    N: 140
    R: (1, -1)
    θc: 19.11°
    Ch: 65.07 Å
    T: 11.27 Å
    dt: 20.71 Å
    rt: 10.36 Å
    electronic_type: semiconducting, type 2

    Change the chirality to :math:`\\mathbf{C}_{h} = (20, 0)`.

    >>> swnt.m = 0
    >>> print(unicode(swnt))
    SWNT(n=20, m=0, element1='C', element2='C', bond=1.42, nz=1)
    n: 20
    m: 0
    t₁: 1
    t₂: -2
    d: 20
    dR: 20
    N: 40
    R: (1, -1)
    θc: 0.00°
    Ch: 49.19 Å
    T: 4.26 Å
    dt: 15.66 Å
    rt: 7.83 Å
    electronic_type: semiconducting, type 1

    """
    # add each attribute in the order I want them to appear in
    # verbose output mode
    _structure_attrs = ['n', 'm', 't1', 't2', 'd', 'dR', 'N', 'R',
                        'chiral_angle', 'Ch', 'T', 'dt', 'rt',
                        'electronic_type']

    def __init__(self, *Ch, nz=None, tube_length=None, Lz=None,
                 fix_Lz=False, **kwargs):

        try:
            n, m = Ch
        except ValueError:
            try:
                n, m = Ch[0]
            except IndexError:
                n = kwargs['n']
                del kwargs['n']
                m = kwargs['m']
                del kwargs['m']

        super(SWNT, self).__init__(**kwargs)

        if tube_length is not None and Lz is None:
            Lz = tube_length

        self.L0 = Lz  # store initial value of Lz

        self.n = n
        self.m = m

        self.fix_Lz = fix_Lz
        if Lz is not None:
            self.nz = 10 * float(Lz) / self.T
        elif nz is not None:
            self.nz = nz
        else:
            self.nz = 1

    def __str__(self):
        """Return nice string representation of `SWNT`."""
        strrep = repr(self)
        if self.verbose:
            strrep += '\n'
            for attr in self._structure_attrs:
                var = attr
                if attr in attr_symbols:
                    var = attr_symbols[attr]
                if attr in attr_strfmt:
                    if attr in attr_units:
                        strrep += \
                            "{}: {}{}\n".format(
                                var, attr_strfmt[attr].format(
                                    getattr(self, attr)), attr_units[attr])
                    else:
                        strrep += "{}: {}\n".format(
                            var, attr_strfmt[attr].format(getattr(self, attr)))
                else:
                    if attr in attr_units:
                        strrep += "{}: {}{}\n".format(
                            var, getattr(self, attr), attr_units[attr])
                    else:
                        strrep += "{}: {}\n".format(
                            var, getattr(self, attr))

        return strrep

    def __repr__(self):
        """Return canonical string representation of `SWNT`."""
        strrep = "SWNT({!r}, element1={!r}, element2={!r}, bond={!r}"
        Ch = (self.n, self.m)
        if self.fix_Lz:
            strrep += ", Lz={!r}, fix_Lz={!r})"
            return strrep.format(Ch, self.element1, self.element2,
                                 self.bond, self.Lz, self.fix_Lz)
        else:
            strrep += ", nz={!r})"
            return strrep.format(Ch, self.element1, self.element2,
                                 self.bond, self.nz)

    def todict(self):
        """Return :class:`~python:dict` of `SWNT` attributes."""
        return dict(n=self.n, m=self.m, nz=self.nz, Lz=self.Lz,
                    fix_Lz=self.fix_Lz, element1=self.element1,
                    element2=self.element2, bond=self.bond)

Nanotube = SWNT
