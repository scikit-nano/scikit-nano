# -*- coding: utf-8 -*-
"""
==============================================================================
SWNT structure class (:mod:`sknano.core.structures.swnt`)
==============================================================================

.. currentmodule:: sknano.core.structures.swnt

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

try:
    from math import gcd
except ImportError:
    from fractions import gcd

import numbers
import numpy as np

from sknano.core import deprecated, deprecate_kwarg
from sknano.core.atoms import Atom, BasisAtom, BasisAtoms
from sknano.core.crystallography import Crystal3DLattice, UnitCell
from sknano.core.refdata import aCC, grams_per_Da
from .base import NanoStructureBase, r_CC_vdw
from .extras import attr_strfmt, attr_symbols, attr_units, \
    get_chiral_indices, get_chiral_type
from .nanotube_bundle import NanotubeBundleBase


__all__ = ['compute_d', 'compute_dR', 'compute_N', 'compute_t1', 'compute_t2',
           'compute_Ch', 'compute_chiral_angle', 'compute_T', 'compute_dt',
           'compute_rt', 'compute_M', 'compute_R', 'compute_R_chiral_angle',
           'compute_symmetry_operation', 'compute_unit_cell_symmetry_params',
           'compute_psi', 'compute_tau',
           'compute_l1', 'compute_l2', 'compute_L',
           'compute_Lx', 'compute_Ly', 'compute_Lz',
           'compute_electronic_type', 'compute_Natoms',
           'compute_Natoms_per_tube', 'compute_Natoms_per_unit_cell',
           'compute_unit_cell_mass', 'compute_linear_mass_density',
           'compute_symmetry_chiral_angle', 'compute_tube_diameter',
           'compute_tube_radius', 'compute_tube_length', 'compute_tube_mass',
           'NanotubeUnitCell', 'SWNTMixin', 'NanotubeMixin', 'SWNTBase',
           'SWNT', 'Nanotube']


class NanotubeUnitCell(UnitCell):
    """Primitive graphene unit cell with 2 atom basis.

    Parameters
    ----------
    bond : :class:`~python:float`, optional
    a : :class:`~python:float`, optional
    gamma : {60, 120}, optional
    basis : {:class:`~python:list`, :class:`~sknano.core.atoms.BasisAtoms`}, \
    optional
    coords : {:class:`~python:list`}, optional
    cartesian : {:class:`~python:bool`}, optional

    """
    def __init__(self, *args, lattice=None, **kwargs):
        basis = kwargs.get('basis')
        if lattice is None or not isinstance(basis, BasisAtoms):
            lattice, basis = \
                self._generate_unit_cell_parameters(*args, **kwargs)
        super().__init__(lattice=lattice, basis=basis)

    def _generate_unit_cell_parameters(self, *Ch, bond=aCC, basis=['C', 'C'],
                                       gutter=r_CC_vdw, wrap_coords=False,
                                       eps=0.01, **kwargs):
        n, m, kwargs = get_chiral_indices(*Ch, **kwargs)
        e1, e2 = basis
        N = compute_N(*Ch)
        T = compute_T(*Ch, bond=bond, length=True)
        rt = compute_rt(*Ch, bond=bond)
        psi, tau, dpsi, dtau = \
            compute_unit_cell_symmetry_params(*Ch, bond=bond)

        a = compute_dt(*Ch, bond=bond) + 2 * gutter
        c = T
        lattice = Crystal3DLattice.hexagonal(a, c)
        basis = BasisAtoms()

        verbose = kwargs.get('verbose', False)
        debug = kwargs.get('debug', False)
        if verbose:
            print('dpsi: {}'.format(dpsi))
            print('dtau: {}\n'.format(dtau))

        for i in range(N):
            for j, element in enumerate((e1, e2), start=1):
                theta = i * psi
                h = i * tau

                if j == 2:
                    theta += dpsi
                    h -= dtau

                x = rt * np.cos(theta)
                y = rt * np.sin(theta)
                z = h

                while z > T - eps:
                    z -= T

                if z < 0:
                    z += T

                xs, ys, zs = \
                    lattice.cartesian_to_fractional([x, y, z])

                if wrap_coords:
                    xs, ys, zs = \
                        lattice.wrap_fractional_coordinate([xs, ys, zs])

                if debug:
                    print('i={}: x, y, z = ({:.6f}, {:.6f}, {:.6f})'.format(
                        i, x, y, z))

                    print('xs, ys, zs = ({:.6f}, {:.6f}, {:.6f})'.format(
                        xs, ys, zs))

                atom = BasisAtom(element, lattice=lattice, xs=xs, ys=ys, zs=zs)
                atom.rezero()

                if verbose:
                    print('Basis Atom:\n{}'.format(atom))

                basis.append(atom)

        return lattice, basis


def compute_d(*Ch):
    """Compute :math:`d=\\gcd{(n, m)}`

    :math:`d` is the **G**\ reatest **C**\ ommon **D**\ ivisor of
    :math:`n` and :math:`m`.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

    Returns
    -------
    gcd : :class:`python:int`
        Greatest Common Divisor of :math:`n` and :math:`m`

    """
    n, m, _ = get_chiral_indices(*Ch)

    return gcd(n, m)


def compute_dR(*Ch):
    """Compute :math:`d_R=\\gcd{(2n + m, 2m + n)}`

    :math:`d_R` is the **G**\ reatest **C**\ ommon **D**\ ivisor of
    :math:`2n + m` and :math:`2m + n`.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

    Returns
    -------
    int
        greatest common divisor of :math:`2n+m` and :math:`2m+n`

    """
    n, m, _ = get_chiral_indices(*Ch)

    return gcd(2 * m + n, 2 * n + m)


def compute_N(*Ch):
    """Compute :math:`N = \\frac{2(n^2+m^2+nm)}{d_R}`.

    :math:`N` is the number of graphene hexagons mapped to a nanotube
    *unit cell*.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

    Returns
    -------
    int
        Number of hexagons per nanotube unit cell:
        :math:`N = \\frac{2(n^2+m^2+nm)}{d_R}`.

    """
    n, m, _ = get_chiral_indices(*Ch)

    dR = compute_dR(n, m)
    try:
        return int(2 * (n ** 2 + m ** 2 + n * m) / dR)
    except ZeroDivisionError:
        return 0


def compute_t1(*Ch):
    """Compute :math:`t_1 = \\frac{2m + n}{d_R}`

    where :math:`d_R = \\gcd{(2n + m, 2m + n)}`.

    The component of the translation vector :math:`\\mathbf{T}`
    along :math:`\\mathbf{a}_1`:

    .. math::

        \\mathbf{T} = t_1\\mathbf{a}_{1} + t_2\\mathbf{a}_2

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

    Returns
    -------
    int
        :math:`t_1`

    """
    n, m, _ = get_chiral_indices(*Ch)

    dR = compute_dR(n, m)
    try:
        return int((2 * m + n) / dR)
    except ZeroDivisionError:
        return 0


def compute_t2(*Ch):
    """Compute :math:`t_2 = -\\frac{2n + m}{d_R}`

    where :math:`d_R = \\gcd{(2n + m, 2m + n)}`.

    The component of the translation vector :math:`\\mathbf{T}`
    along :math:`\\mathbf{a}_2`:

    .. math::

        \\mathbf{T} = t_1\\mathbf{a}_1 + t_2\\mathbf{a}_2

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

    Returns
    -------
    int
        :math:`t_2`

    """
    n, m, _ = get_chiral_indices(*Ch)

    dR = compute_dR(n, m)
    try:
        return -int((2 * n + m) / dR)
    except ZeroDivisionError:
        return 0


def compute_Ch(*Ch, bond=None, **kwargs):
    """Compute nanotube circumference :math:`|\\mathbf{C}_{h}|` in \
    **\u212b**.

    .. math::

       |\\mathbf{C}_h| = a\\sqrt{n^2 + m^2 + nm} =
       \\sqrt{3}a_{\\mathrm{CC}}\\sqrt{n^2 + m^2 + nm}

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

    Returns
    -------
    float
        Nanotube circumference :math:`|\\mathbf{C}_h|` in \u212b.

    """
    n, m, _ = get_chiral_indices(*Ch)

    if bond is None:
        bond = aCC

    return bond * np.sqrt(3 * (n ** 2 + m ** 2 + n * m))


def compute_chiral_angle(*Ch, degrees=True, **kwargs):
    """Compute chiral angle :math:`\\theta_c`.

    .. math::

       \\theta_c = \\tan^{-1}\\left(\\frac{\\sqrt{3} m}{2n + m}\\right)

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    degrees : bool, optional
        If `True`, return angle in degrees.

    Returns
    -------
    float
        Chiral angle :math:`\\theta_{c}` in
        degrees (default) or radians (if `degrees=False`).

    """
    n, m, _ = get_chiral_indices(*Ch)

    theta = np.arctan(np.sqrt(3) * m / (2 * n + m))
    # return np.arccos((2*n + m) / (2 * np.sqrt(n**2 + m**2 + n*m)))
    if degrees or kwargs.get('rad2deg', False):
        return np.degrees(theta)
    else:
        return theta


def compute_T(*Ch, bond=None, length=True, **kwargs):
    """Compute length of nanotube unit cell :math:`|\\mathbf{T}|` in \
    \u212b.

    .. math::

       |\\mathbf{T}| = \\frac{\\sqrt{3} |\\mathbf{C}_{h}|}{d_{R}}
       = \\frac{\\sqrt{3}a\\sqrt{n^2 + m^2 + nm}}{d_{R}}
       = \\frac{3a_{\\mathrm{CC}}\\sqrt{n^2 + m^2 + nm}}{d_{R}}

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b
    length : bool, optional
        Compute the magnitude (i.e., length) of the translation vector.

    Returns
    -------
    float or 2-tuple of ints
        If `length` is `True`, then
        return the length of unit cell in \u212b.

        If `length` is `False`, return the componets of the
        translation vector as a 2-tuple of ints
        (:math:`t_1`, :math:`t_2`).

    """
    n, m, _ = get_chiral_indices(*Ch)

    if length:
        if bond is None:
            bond = aCC

        Ch = compute_Ch(n, m, bond=bond)
        dR = compute_dR(n, m)

        try:
            return np.sqrt(3) * Ch / dR
        except (FloatingPointError, ZeroDivisionError):
            return 0
    else:
        t1 = compute_t1(n, m)
        t2 = compute_t2(n, m)

        return (t1, t2)


def compute_dt(*Ch, bond=None, **kwargs):
    """Compute nanotube diameter :math:`d_t` in \u212b.

    .. math::

       d_t = \\frac{|\\mathbf{C}_h|}{\\pi}

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

    Returns
    -------
    float
        Nanotube diameter :math:`d_t` in \u212b.

    """
    n, m, _ = get_chiral_indices(*Ch)

    Ch = compute_Ch(n, m, bond=bond)
    return Ch / np.pi


def compute_rt(*Ch, bond=None, **kwargs):
    """Compute nanotube radius :math:`r_t` in \u212b.

    .. math::

        r_t = \\frac{|\\mathbf{C}_h|}{2\\pi}


    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

    Returns
    -------
    float
        Nanotube radius :math:`r_t` in \u212b.

    """
    n, m, _ = get_chiral_indices(*Ch)

    Ch = compute_Ch(n, m, bond=bond)
    return Ch / (2 * np.pi)


def compute_M(*Ch):
    """Compute :math:`M = mp - nq`

    :math:`M` is the number of multiples of the translation vector
    :math:`\\mathbf{T}` in the vector :math:`N\\mathbf{R}`.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

    Returns
    -------
    int
        :math:`M = mp - nq`

    """
    n, m, _ = get_chiral_indices(*Ch)

    p, q = compute_R(n, m)
    return m * p - n * q


def compute_R(*Ch, bond=None, length=False, **kwargs):
    """Compute symmetry vector :math:`\\mathbf{R} = (p, q)`.

    The *symmetry vector* is any lattice vector of the unfolded graphene
    layer that represents a *symmetry operation* of the nanotube. The
    symmetry vector :math:`\\mathbf{R}` can be written as:

    .. math::

        \\mathbf{R} = p\\mathbf{a}_1 + q\\mathbf{a}_2

    where :math:`p` and :math:`q` are integers.
    The *symmetry vector* represents a *symmetry operation* of the nanotube
    which arises as a *screw translation*, which is a combination of
    a rotation :math:`\\psi` and translation :math:`\\tau`. The symmetry
    operation of the nanotube can be written as:

    .. math::

        R = (\\psi|\\tau)

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b
    length : bool, optional
        If `True`, return :math:`|\\mathbf{R}|`.

    Returns
    -------
    (p, q) : tuple
        2-tuple of ints -- components of :math:`\\mathbf{R}`.
    float
        Length of :math:`\\mathbf{R}` (:math:`|\\mathbf{R}|`) if `length`
        is `True` in units of **\u212b**.

    """
    n, m, _ = get_chiral_indices(*Ch)

    N = compute_N(n, m)
    t1 = compute_t1(n, m)
    t2 = compute_t2(n, m)

    p = q = 0
    for i in range(0, t1 + n + 1):
        for j in range(t2, m + 1):
            R = t1 * j - t2 * i
            if R == 1:
                M = m * i - n * j
                if M > 0 and M <= N:
                    p = i
                    q = j

    if length:
        if bond is None:
            bond = aCC

        return bond * np.sqrt(3 * (p ** 2 + q ** 2 + p * q))
    else:
        return (p, q)


def compute_R_chiral_angle(*Ch, degrees=True, **kwargs):
    """Compute "chiral angle" of symmetry vector :math:`\\theta_R`.

    .. math::

        \\theta_R = \\tan^{-1}\\left(\\frac{\\sqrt{3}q}{2p + q}\\right)

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    degrees : bool, optional
        If `True`, return angle in degrees

    Returns
    -------
    float
        Chiral angle of *symmetry vector* :math:`\\theta_R` in
        degrees (default) or radians (if `degrees=False`).


    """
    n, m, _ = get_chiral_indices(*Ch)

    p, q = compute_R(n, m)
    theta = np.arctan((np.sqrt(3) * q) / (2 * p + q))
    if degrees or kwargs.get('rad2deg', False):
        return np.degrees(theta)
    else:
        return theta


def compute_symmetry_operation(*Ch, bond=None):
    """Compute symmetry operation :math:`(\\psi|\\tau)`.

    The *symmetry vector* `R` represents a *symmetry
    operation* of the nanotube which arises as a
    *screw translation*--a combination of a rotation
    :math:`\\psi` and translation :math:`\\tau`.
    The symmetry operation of the nanotube can be written as:

    .. math::

        R = (\\psi|\\tau)

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

    Returns
    -------
    (psi, tau) : tuple
        2-tuple of floats -- :math:`\\psi` in radians and
        :math:`\\tau` in \u212b.

    """
    n, m, _ = get_chiral_indices(*Ch)

    psi = compute_psi(n, m)
    tau = compute_tau(n, m, bond=bond)
    return (psi, tau)


def compute_unit_cell_symmetry_params(*Ch, bond=None):
    """Tuple of `SWNT` unit cell *symmetry parameters*."""
    psi, tau = compute_symmetry_operation(*Ch, bond=bond)
    aCh = compute_chiral_angle(*Ch, degrees=False)
    rt = compute_rt(*Ch, bond=bond)
    dpsi = bond * np.cos(np.pi / 6 - aCh) / rt
    dtau = bond * np.sin(np.pi / 6 - aCh)
    return psi, tau, dpsi, dtau


def compute_psi(*Ch):
    """Compute rotation component of symmetry operation \
    :math:`\\psi` in **radians**.

    .. math::

        \\psi = \\frac{2\\pi}{N}

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

    Returns
    -------
    float
        Rotation component of symmetry operation :math:`\\psi`
        in **radians**.

    """
    n, m, _ = get_chiral_indices(*Ch)

    N = compute_N(n, m)
    try:
        return 2 * np.pi / N
    except (FloatingPointError, ZeroDivisionError):
        return 0


def compute_tau(*Ch, bond=None, **kwargs):
    """Compute translation component of symmetry operation \
    :math:`\\tau` in **\u212b**.

    .. math::

        \\tau = \\frac{M|\\mathbf{T}|}{N}

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

    Returns
    -------
    float
        Translation component of symmetry operation :math:`\\tau`
        in **\u212b**.

    """
    n, m, _ = get_chiral_indices(*Ch)

    M = compute_M(n, m)
    N = compute_N(n, m)
    T = compute_T(n, m, bond=bond)
    try:
        return M * T / N
    except ZeroDivisionError:
        return 0


def compute_electronic_type(*Ch):
    """Compute nanotube electronic type.

    .. versionadded:: 0.2.7

    The electronic type is determined as follows:

    if :math:`(2n + m)\\,\\mathrm{mod}\\,3=0`, the nanotube is
    **metallic**.

    if :math:`(2n + m)\\,\\mathrm{mod}\\,3=1`, the nanotube is
    **semiconducting, type 1**.

    if :math:`(2n + m)\\,\\mathrm{mod}\\,3=2`, the nanotube is
    **semiconducting, type 2**.

    The :math:`x\\,\\mathrm{mod}\\,y` notation is mathematical
    shorthand for the *modulo* operation, which computes the
    **remainder** of the division of :math:`x` by :math:`y`.
    So, for example, all *armchair* nanotubes must be metallic
    since the chiral indices satisfy: :math:`2n + m = 2n + n = 3n` and
    therefore :math:`3n\\,\\mathrm{mod}\\,3` i.e. the remainder of the
    division of :math:`3n/3=n` is always zero.

    .. note::
       Mathematically, :math:`(2n + m)\\,\\mathrm{mod}\\,3` is equivalent
       to :math:`(n - m)\\,\\mathrm{mod}\\,3` when distinguishing
       between metallic and semiconducting. However, when
       distinguishing between semiconducting types,
       one must be careful to observe the following convention:

       * Semiconducting, **type 1** means:

         * :math:`(2n + m)\\,\\mathrm{mod}\\,3=1`
         * :math:`(n - m)\\,\\mathrm{mod}\\,3=2`

       * Semiconducting, **type 2** means:

         * :math:`(2n + m)\\,\\mathrm{mod}\\,3=2`
         * :math:`(n - m)\\,\\mathrm{mod}\\,3=1`

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

    Returns
    -------
    :class:`~python:str`

    """
    n, m, _ = get_chiral_indices(*Ch)

    if (2 * n + m) % 3 == 1:
        return 'semiconducting, type 1'
    elif (2 * n + m) % 3 == 2:
        return 'semiconducting, type 2'
    else:
        return 'metallic'


@deprecate_kwarg(kwarg='nz', since='0.4.0', alternative='n3')
def compute_Natoms(*Ch, n3=1):
    """Compute :math:`N_{\\mathrm{atoms/tube}}`

    .. math::

        N_{\\mathrm{atoms/tube}} = N_{\\mathrm{atoms/cell}} \\times
        n_{3-\\mathrm{cells}}

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    n3 : {int, float}
        Number of nanotube unit cells

    Returns
    -------
    int
        :math:`N_{\\mathrm{atoms/tube}}`
    """
    n, m, _ = get_chiral_indices(*Ch)

    Natoms_per_unit_cell = compute_Natoms_per_unit_cell(n, m)
    return int(Natoms_per_unit_cell * n3)


@deprecate_kwarg(kwarg='nz', since='0.4.0', alternative='n3')
def compute_Natoms_per_tube(*Ch, n3=1):
    """Compute :math:`N_{\\mathrm{atoms/tube}}`

    .. math::

        N_{\\mathrm{atoms/tube}} = N_{\\mathrm{atoms/cell}} \\times
        n_{3-\\mathrm{cells}}

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    n3 : {int, float}
        Number of nanotube unit cells

    Returns
    -------
    int
        :math:`N_{\\mathrm{atoms/tube}}`
    """
    return compute_Natoms(*Ch, n3=n3)


def compute_Natoms_per_unit_cell(*Ch):
    """Compute :math:`N_{\mathrm{atoms/cell}} = 2N`.

    .. math::

        N_{\\mathrm{atoms}} = 2N = \\frac{4(n^2 + m^2 + nm)}{d_R}

    where :math:`N` is the number of graphene hexagons mapped to the
    nanotube unit cell.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.

    Returns
    -------
    `2N` : int
        Number of atoms in nanotube *unit cell*:
        N_{\\mathrm{atoms}} = 2N = \\frac{4(n^2 + m^2 + nm)}{d_R}


    """
    n, m, _ = get_chiral_indices(*Ch)

    N = compute_N(n, m)
    return 2 * N


def compute_unit_cell_mass(*Ch, element1=None, element2=None, **kwargs):
    """Compute nanotube unit cell mass in *Daltons*/*atomic mass units* (amu) \
        units.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2

    Returns
    -------
    float
        Unit cell mass in **Daltons**.

    Notes
    -----

    .. todo::

        Handle different elements and perform accurate calculation by
        determining number of atoms of each element.

    """
    n, m, _ = get_chiral_indices(*Ch)

    N = compute_N(n, m)

    if element1 is None:
        element1 = 'C'
    if element2 is None:
        element2 = 'C'

    mass = N * (Atom(element1).mass + Atom(element2).mass)
    return mass


def compute_linear_mass_density(*Ch, bond=None, element1=None, element2=None,
                                **kwargs):
    """Compute nanotube linear mass density (mass per unit length) in \
    **grams/Å**.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

    Returns
    -------
    float
        Linear mass density in units of **g/Å**.

    """
    n, m, _ = get_chiral_indices(*Ch)

    mass = compute_unit_cell_mass(n, m, element1=element1, element2=element2,
                                  **kwargs)
    T = compute_T(n, m, bond=bond, length=True, **kwargs)

    try:
        linear_mass_density = mass / T
        # there are 1.6605e-24 grams / Da
        linear_mass_density *= grams_per_Da
        return linear_mass_density
    except ZeroDivisionError:
        return 0


@deprecated(since='0.4.0', alternative='compute_l1')
def compute_Lx(*Ch, **kwargs):
    """Compute the axis-aligned length along the `x`-axis in **Angstroms**.

    Calculated as:

    .. math::

       L_x = n_x * (d_t + 2 r_{\\mathrm{vdW}})

    """
    return compute_l1(*Ch, **kwargs)


@deprecated(since='0.4.0', alternative='compute_l2')
def compute_Ly(*Ch, **kwargs):
    """Compute the axis-aligned length along the `y`-axis in **Angstroms**.

    Calculated as:

    .. math::

       L_y = n_y * (d_t + 2 r_{\\mathrm{vdW}})

    """
    return compute_l2(*Ch, **kwargs)


@deprecated(since='0.4.0', alternative='compute_L')
def compute_Lz(*Ch, **kwargs):
    """Compute the axis-aligned length along the `z`-axis in **Angstroms**.

    :math:`L_z = L_{\\mathrm{tube}}` in **Angstroms**.

    .. versionchanged:: 0.4.0

       Changed units from nanometers to **Angstroms**

    .. math::

        L_z = n_z |\\mathbf{T}|


    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nz : {int, float}
        Number of nanotube unit cells
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

    Returns
    -------
    float
        :math:`L_z = L_{\\mathrm{tube}}` in **Angstroms**

    """
    return compute_L(*Ch, **kwargs)


def compute_l1(*Ch, n1=1, bond=None, gutter=r_CC_vdw):
    """Compute the length along the :math:`\\mathbf{a}_1` lattice vector in \
        **Angstroms**.

    Calculated as:

    .. math::

       \\ell_1 = n_1 * (d_t + 2 r_{\\mathrm{vdW}})

    """
    return n1 * (compute_dt(*Ch, bond=bond) + 2 * gutter)


def compute_l2(*Ch, n2=1, bond=None, gutter=r_CC_vdw):
    """Compute the length along the :math:`\\mathbf{a}_2` lattice vector in \
        **Angstroms**.

    Calculated as:

    .. math::

       \\ell_2 = n_2 * (d_t + 2 r_{\\mathrm{vdW}})

    """
    return n2 * (compute_dt(*Ch, bond=bond) + 2 * gutter)


def compute_L(*Ch, n3=1, bond=None, **kwargs):
    """Compute the length along the :math:`\\mathbf{a}_3` lattice vector in \
        **Angstroms**.

    :math:`\\ell_3 = \\ell_{\\mathrm{tube}}` in **Angstroms**.

    .. math::

       \\ell_3 = n_3 |\\mathbf{T}|


    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    n3 : {int, float}
        Number of nanotube unit cells all :math:`\\mathbf{a}_3`
    bond : float, optional
        Distance between nearest neighbor atoms (i.e., bond length).
        Must be in units of **\u212b**. Default value is
        the carbon-carbon bond length in graphite, approximately
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.42` \u212b

    Returns
    -------
    float
        :math:`\\ell_3 = \\ell_{\\mathrm{tube}}` in **Angstroms**

    """
    n, m, _ = get_chiral_indices(*Ch)

    if not (isinstance(n3, numbers.Real) or n3 > 0):
        raise TypeError('Expected a real, positive number')

    T = compute_T(n, m, bond=bond, **kwargs)
    return n3 * T


def compute_symmetry_chiral_angle(*Ch, degrees=True):
    """Alias for :func:`compute_R_chiral_angle`."""
    return compute_R_chiral_angle(*Ch, degrees=degrees)


def compute_tube_diameter(*Ch, bond=None, **kwargs):
    """Alias for :func:`compute_dt`"""
    return compute_dt(*Ch, bond=bond, **kwargs)


def compute_tube_radius(*Ch, bond=None, **kwargs):
    """Alias for :func:`compute_rt`"""
    return compute_rt(*Ch, bond=bond, **kwargs)


@deprecate_kwarg(kwarg='nz', since='0.4.0', alternative='n3')
def compute_tube_length(*Ch, n3=1, bond=None, **kwargs):
    """Alias for :func:`compute_L`"""
    return compute_L(*Ch, n3=n3, bond=bond, **kwargs)


@deprecate_kwarg(kwarg='nz', since='0.4.0', alternative='n3')
def compute_tube_mass(*Ch, n3=1, element1=None, element2=None, **kwargs):
    """Compute nanotube mass in **grams**.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    n3 : {int, float}
        Number of nanotube unit cells
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2

    Returns
    -------
    float
        Nanotube mass in **grams**.

    Notes
    -----

    .. todo::

        Handle different elements and perform accurate calculation by
        determining number of atoms of each element.

    """
    n, m, _ = get_chiral_indices(*Ch)

    if not (isinstance(n3, numbers.Real) or n3 > 0):
        raise TypeError('Expected a real, positive number')

    Natoms = compute_Natoms(n, m, n3=n3)

    if element1 is None:
        element1 = 'C'
    if element2 is None:
        element2 = 'C'

    atom1 = Atom(element1)
    atom2 = Atom(element2)

    mass = Natoms * (atom1.mass + atom2.mass) / 2
    # there are 1.6605e-24 grams / Da
    mass *= grams_per_Da
    return mass


class SWNTMixin:
    """Mixin class for nanotube classes."""
    @property
    def n(self):
        """Chiral index :math:`n`.

        The component of the chiral vector :math:`\\mathbf{C}_h`
        along :math:`\\mathbf{a}_1`:

        .. math::

           \\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)

        """
        return self._n

    @n.setter
    def n(self, value):
        """Set chiral index :math:`n`"""
        if not (isinstance(value, numbers.Real) or value >= 0):
            raise TypeError('Expected an integer.')
        self._n = int(value)

    @n.deleter
    def n(self):
        del self._n

    @property
    def m(self):
        """Chiral index :math:`m`.

        The component of the chiral vector :math:`\\mathbf{C}_h`
        along :math:`\\mathbf{a}_2`:

        .. math::

           \\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)

        """
        return self._m

    @m.setter
    def m(self, value):
        """Set chiral index :math:`m`"""
        if not (isinstance(value, numbers.Real) or value >= 0):
            raise TypeError('Expected an integer.')
        self._m = int(value)

    @m.deleter
    def m(self):
        del self._m

    @property
    def d(self):
        """:math:`d=\\gcd{(n, m)}`

        :math:`d` is the **G**\ reatest **C**\ ommon **D**\ ivisor of
        :math:`n` and :math:`m`.

        """
        return compute_d(self.n, self.m)

    @property
    def dR(self):
        """:math:`d_R=\\gcd{(2n + m, 2m + n)}`

        :math:`d_R` is the **G**\ reatest **C**\ ommon **D**\ ivisor of
        :math:`2n + m` and :math:`2m + n`.

        """
        return compute_dR(self.n, self.m)

    @property
    def N(self):
        """Number of graphene hexagons in nanotube *unit cell*.

        .. math::

           N = \\frac{4(n^2 + m^2 + nm)}{d_R}

        """
        return compute_N(self.n, self.m)

    @property
    def t1(self):
        """:math:`t_{1} = \\frac{2m + n}{d_{R}}`

        where :math:`d_R = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_1`:

        .. math::

           \\mathbf{T} = t_1\\mathbf{a}_{1} + t_2\\mathbf{a}_2

        """
        return compute_t1(self.n, self.m)

    @property
    def t2(self):
        """:math:`t_2 = -\\frac{2n + m}{d_R}`

        where :math:`d_R = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_2`:

        .. math::

           \\mathbf{T} = t_1\\mathbf{a}_1 + t_2\\mathbf{a}_2

        """
        return compute_t2(self.n, self.m)

    @property
    def Ch_vec(self):
        """SWNT chiral vector."""
        return (self.n, self.m)

    @property
    def Ch(self):
        """SWNT circumference :math:`|\\mathbf{C}_h|` in **\u212b**"""
        return compute_Ch(self.n, self.m, bond=self.bond)

    @property
    def dt(self):
        """Nanotube diameter :math:`d_t = \\frac{|\\mathbf{C}_h|}{\\pi}` \
            in \u212b."""
        return compute_dt(self.n, self.m, bond=self.bond)

    @property
    def rt(self):
        """Nanotube radius :math:`r_t = \\frac{|\\mathbf{C}_h|}{2\\pi}` \
            in \u212b."""
        return compute_rt(self.n, self.m, bond=self.bond)

    @property
    def chiral_angle(self):
        """Chiral angle :math:`\\theta_c` in **degrees**.

        .. math::

           \\theta_c = \\tan^{-1}\\left(\\frac{\\sqrt{3} m}{2n + m}\\right)

        """
        return compute_chiral_angle(self.n, self.m)

    @property
    def chiral_type(self):
        """`SWNT` chiral type."""
        return get_chiral_type((self.n, self.m))

    @property
    def Tvec(self):
        """`SWNT` translation vector."""
        return (self.t1, self.t2)

    @property
    def T(self):
        """Length of nanotube unit cell :math:`|\\mathbf{T}|` in \u212b.

        .. math::

           |\\mathbf{T}| = \\frac{\\sqrt{3} |\\mathbf{C}_{h}|}{d_{R}}

        """
        return compute_T(self.n, self.m, bond=self.bond, length=True)

    @property
    def M(self):
        """:math:`M = np - nq`

        :math:`M` is the number of multiples of the translation vector
        :math:`\\mathbf{T}` in the vector :math:`N\\mathbf{R}`.

        """
        return compute_M(self.n, self.m)

    @property
    def R(self):
        """Symmetry vector :math:`\\mathbf{R} = (p, q)`.

        .. math::

           \\mathbf{R} = p\\mathbf{a}_1 + q\\mathbf{a}_2

        """
        return compute_R(self.n, self.m, bond=self.bond, length=False)

    @property
    @deprecated(since='0.4.0', alternative='n3', obj_type='attribute')
    def nz(self):
        """Number of nanotube unit cells along the :math:`z`-axis."""
        return self.n3

    @nz.setter
    def nz(self, value):
        self.n3 = value

    def _update_nz(self):
        self._update_n3()

    @property
    def n3(self):
        """Number of nanotube unit cells along the :math:`\\mathbf{a}_3` \
            lattice vector."""
        return self._n3

    @n3.setter
    def n3(self, value):
        if not (isinstance(value, numbers.Real) or value > 0):
            raise TypeError('Expected a real, positive number.')
        self._n3 = value
        if self._integral_n3:
            self._n3 = int(np.ceil(value))

    def _update_n3(self):
        try:
            self.n3 = self.n3
        except AttributeError:
            pass

    @property
    def electronic_type(self):
        """SWNT electronic type.

        .. versionadded:: 0.2.7

        The electronic type is determined as follows:

        if :math:`(2n + m)\\,\\mathrm{mod}\\,3=0`, the nanotube is
        **metallic**.

        if :math:`(2n + m)\\,\\mathrm{mod}\\,3=1`, the nanotube is
        **semiconducting, type 1**.

        if :math:`(2n + m)\\,\\mathrm{mod}\\,3=2`, the nanotube is
        **semiconducting, type 2**.

        The :math:`x\\,\\mathrm{mod}\\,y` notation is mathematical
        shorthand for the *modulo* operation, which computes the
        **remainder** of the division of :math:`x` by :math:`y`.
        So, for example, all *armchair* nanotubes must be metallic
        since the chiral indices satisfy: :math:`2n + m = 2n + n = 3n` and
        therefore :math:`3n\\,\\mathrm{mod}\\,3` i.e. the remainder of the
        division of :math:`3n/3=n` is always zero.

        .. note::
           Mathematically, :math:`(2n + m)\\,\\mathrm{mod}\\,3` is equivalent
           to :math:`(n - m)\\,\\mathrm{mod}\\,3` when distinguishing
           between metallic and semiconducting. However, when
           distinguishing between semiconducting types,
           one must be careful to observe the following convention:

           * Semiconducting, **type 1** means:

             * :math:`(2n + m)\\,\\mathrm{mod}\\,3=1`
             * :math:`(n - m)\\,\\mathrm{mod}\\,3=2`

           * Semiconducting, **type 2** means:

             * :math:`(2n + m)\\,\\mathrm{mod}\\,3=2`
             * :math:`(n - m)\\,\\mathrm{mod}\\,3=1`

        """
        return compute_electronic_type(self.n, self.m)

    @property
    @deprecated(since='0.4.0', alternative='L', obj_type='attribute')
    def Lz(self):
        """SWNT length :math:`L_z = L_{\\mathrm{tube}}` in **Angstroms**."""
        return self.L

    @property
    @deprecated(since='0.4.0', alternative='fix_L', obj_type='attribute')
    def fix_Lz(self):
        """:class:`~python:bool` indicating whether :attr:`SWNTMixin.Lz` is \
            fixed or calculated."""
        return self.fix_L

    @fix_Lz.setter
    def fix_Lz(self, value):
        self.fix_L = value

    @property
    def L(self):
        """SWNT length :math:`L = L_{\\mathrm{tube}}` in **Angstroms**."""
        return self.n3 * self.T

    @property
    def fix_L(self):
        """:class:`~python:bool` indicating whether :attr:`SWNTMixin.L` is \
            fixed or calculated."""
        return self._fix_L

    @fix_L.setter
    def fix_L(self, value):
        if not isinstance(value, bool):
            raise TypeError('Expected `True` or `False`')
        self._fix_L = value
        self._integral_n3 = False if self.fix_L else True
        self._update_n3()
        self._update_fmtstr()

    def _update_fmtstr(self):
        if self.fix_L:
            self.fmtstr = self.fmtstr.replace("n3={n3!r}",
                                              "L={L!r}, fix_L=True")
        else:
            self.fmtstr = self.fmtstr.replace("L={L!r}, fix_L=True",
                                              "n3={n3!r}")

    @property
    def Natoms(self):
        """Number of atoms in nanotube.

        .. versionchanged:: 0.3.0

           **Returns total number of atoms per nanotube.**
           Use :attr:`~SWNT.Natoms_per_unit_cell` to get the number of
           atoms per unit cell.

        .. math::

           N_{\\mathrm{atoms}} = 2N\\times n_z =
           \\frac{4(n^2 + m^2 + nm)}{d_R}\\times n_z

        where :math:`N` is the number of graphene hexagons mapped to the
        nanotube unit cell and :math:`n_3` is the number of unit cells.

        """
        return compute_Natoms(self.n, self.m, n3=self.n3)

    @property
    def Natoms_per_unit_cell(self):
        """Number of atoms in nanotube unit cell.

        .. math::

           N_{\\mathrm{atoms}} = 2N = \\frac{4(n^2 + m^2 + nm)}{d_R}

        where :math:`N` is the number of graphene hexagons mapped to the
        nanotube unit cell.

        """
        return compute_Natoms_per_unit_cell(self.n, self.m)

    @property
    def Natoms_per_tube(self):
        """Number of atoms in nanotube :math:`N_{\\mathrm{atoms/tube}}`."""
        return self.Natoms

    @property
    def linear_mass_density(self):
        """Linear mass density of nanotube in g/Å."""
        return compute_linear_mass_density(self.n, self.m, bond=self.bond,
                                           element1=self.element1,
                                           element2=self.element2)

    @property
    def tube_length(self):
        """Alias for :attr:`SWNT.L`"""
        return self.L

    @property
    def mass(self):
        """SWNT mass in **grams**."""
        return compute_tube_mass(self.n, self.m, n3=self.n3,
                                 element1=self.element1,
                                 element2=self.element2)

    @property
    def tube_mass(self):
        """An alias for :attr:`~SWNTMixin.mass`."""
        return self.mass

    @property
    def unit_cell_mass(self):
        """Unit cell mass in atomic mass units."""
        return compute_unit_cell_mass(self.n, self.m,
                                      element1=self.element1,
                                      element2=self.element2)

    @property
    def unit_cell_symmetry_params(self):
        """Tuple of `SWNT` unit cell *symmetry parameters*."""
        psi, tau = compute_symmetry_operation(self.n, self.m, bond=self.bond)
        aCh = compute_chiral_angle(self.n, self.m, degrees=False)
        dpsi = self.bond * np.cos(np.pi / 6 - aCh) / self.rt
        dtau = self.bond * np.sin(np.pi / 6 - aCh)

        return psi, tau, dpsi, dtau

NanotubeMixin = SWNTMixin


class SWNTBase(SWNTMixin, NanoStructureBase):
    """Base SWNT structure class."""
    # add each attribute in the order I want them to appear in
    # verbose output mode
    _structure_attrs = ['n', 'm', 't1', 't2', 'd', 'dR', 'N', 'R',
                        'chiral_angle', 'Ch', 'T', 'dt', 'rt',
                        'electronic_type']

    @deprecate_kwarg(kwarg='nz', since='0.4.0', alternative='n3')
    @deprecate_kwarg(kwarg='Lz', since='0.4.0', alternative='L')
    @deprecate_kwarg(kwarg='fix_Lz', since='0.4.0', alternative='fix_L')
    def __init__(self, *Ch, n3=None, gutter=None, L=None, fix_L=False,
                 wrap_coords=False, **kwargs):

        n, m, kwargs = get_chiral_indices(*Ch, **kwargs)
        L = kwargs.pop('tube_length', L)

        super().__init__(**kwargs)

        if gutter is None:
            gutter = self.vdw_radius
        self.gutter = gutter
        self.wrap_coords = wrap_coords

        self.L0 = L  # store initial value of L

        self.n = n
        self.m = m

        self.fix_L = fix_L
        if L is not None:
            self.n3 = float(L) / self.T
        elif n3 is not None:
            self.n3 = n3
        else:
            self.n3 = 1

        # self.generate_unit_cell()
        self.unit_cell = \
            NanotubeUnitCell((self.n, self.m), bond=self.bond,
                             basis=self.basis, gutter=self.gutter,
                             wrap_coords=wrap_coords)
        self.scaling_matrix = [1, 1, int(np.ceil(self.n3))]
        fmtstr = ", ".join(("{Ch!r}", super().fmtstr))
        if self.fix_L:
            fmtstr = ", ".join((fmtstr, "L={L!r}", "fix_L=True"))
        else:
            fmtstr = ", ".join((fmtstr, "n3={n3!r}"))
        self.fmtstr = ", ".join((fmtstr, "gutter={gutter!r}",
                                 "wrap_coords={wrap_coords!r}"))

    def __str__(self):
        """Return nice string representation of `SWNT`."""
        fmtstr = repr(self)
        if self.verbose:
            fmtstr += '\n'
            for attr in self._structure_attrs:
                var = attr
                if attr in attr_symbols:
                    var = attr_symbols[attr]
                if attr in attr_strfmt:
                    if attr in attr_units:
                        fmtstr += \
                            "{}: {}{}\n".format(
                                var, attr_strfmt[attr].format(
                                    getattr(self, attr)), attr_units[attr])
                    else:
                        fmtstr += "{}: {}\n".format(
                            var, attr_strfmt[attr].format(getattr(self, attr)))
                else:
                    if attr in attr_units:
                        fmtstr += "{}: {}{}\n".format(
                            var, getattr(self, attr), attr_units[attr])
                    else:
                        fmtstr += "{}: {}\n".format(
                            var, getattr(self, attr))

        return fmtstr

    def todict(self):
        """Return :class:`~python:dict` of `SWNT` attributes."""
        attr_dict = super().todict()
        attr_dict.update(dict(Ch=(self.n, self.m),
                              n3=self.n3, L=self.L, fix_L=self.fix_L,
                              gutter=self.gutter,
                              wrap_coords=self.wrap_coords))
        return attr_dict


class SWNT(NanotubeBundleBase, SWNTBase):
    """SWNT structure class.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of integers (i.e., *Ch = ((n, m)) or
        2 integers (i.e., *Ch = (n, m) specifying the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    n1 : :class:`python:int`, optional
        Number of nanotubes along the :math:`x` axis
    n2 : :class:`python:int`, optional
        Number of nanotubes along the :math:`y` axis
    n3 : :class:`python:int`, optional
        Number of repeat unit cells in the :math:`z` direction, along
        the *length* of the nanotube.
    basis : {:class:`python:list`}, optional
        List of :class:`python:str`\ s of element symbols or atomic number
        of the two atom basis (default: ['C', 'C'])

        .. versionadded:: 0.3.10

    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.Atom` 1 and 2

        .. deprecated:: 0.3.10
           Use `basis` instead

    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    L : float, optional
        Length of nanotube in units of **Angstroms**.
        Overrides the `n3` value.

        .. versionadded:: 0.2.5

        .. versionchanged:: 0.4.0

           Changed units from nanometers to **Angstroms**

    tube_length : float, optional
        Length of nanotube in units of **Angstroms**.
        Overrides the `n3` value.

        .. deprecated:: 0.2.5
           Use `L` instead

    fix_L : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.

        .. versionadded:: 0.2.6

    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.core.structures import SWNT

    Create a SWNT with :math:`\\mathbf{C}_{h} = (10, 10)` chirality.

    >>> swnt = SWNT((10, 10), verbose=True)
    >>> print(swnt)
    SWNT((10, 10), n3=1)
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
    >>> print(swnt)
    SWNT((20, 10), n3=1)
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
    >>> print(swnt)
    SWNT((20, 0), n3=1)
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

    The next example defines a :math:`\\mathbf{C}_{h} = (10, 10)`
    hexagonally close packed (*hcp*)
    :math:`5\\times 3\\times 10` :class:`SWNT` bundle.

    >>> swnt_bundle = SWNT((10, 10), n1=5, n2=3, n3=10, bundle_packing='hcp')
    >>> print(swnt_bundle)
    SWNT((10, 10), n1=5, n2=3, n3=10, basis=['C', 'C'], bond=1.42,
    bundle_packing='hcp', bundle_geometry=None)

    """
    pass

Nanotube = SWNT
