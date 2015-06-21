# -*- coding: utf-8 -*-
"""
===============================================================================
Compute functions (:mod:`sknano.structures._compute_funcs`)
===============================================================================

.. currentmodule:: sknano.structures._compute_funcs

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from fractions import gcd
import numbers
import numpy as np

from sknano.core.atoms import Atom
from sknano.core.refdata import CCbond, dVDW, grams_per_Da
from ._extras import get_chiral_indices

__all__ = ['compute_d', 'compute_dR', 'compute_N', 'compute_t1', 'compute_t2',
           'compute_Ch', 'compute_chiral_angle', 'compute_T', 'compute_dt',
           'compute_rt', 'compute_M', 'compute_R', 'compute_R_chiral_angle',
           'compute_symmetry_operation', 'compute_psi', 'compute_tau',
           'compute_Lx', 'compute_Ly', 'compute_Lz',
           'compute_electronic_type', 'compute_Natoms',
           'compute_Natoms_per_tube', 'compute_Natoms_per_unit_cell',
           'compute_unit_cell_mass', 'compute_linear_mass_density',
           'compute_bundle_density', 'compute_symmetry_chiral_angle',
           'compute_tube_diameter', 'compute_tube_radius',
           'compute_tube_length', 'compute_tube_mass']


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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')
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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')
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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')
    dR = compute_dR(n, m)
    try:
        return int(2 * (n**2 + m**2 + n * m) / dR)
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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')
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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')
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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')
    if bond is None:
        bond = CCbond

    return bond * np.sqrt(3 * (n**2 + m**2 + n * m))


def compute_chiral_angle(*Ch, rad2deg=True):
    """Compute chiral angle :math:`\\theta_c`.

    .. math::

       \\theta_c = \\tan^{-1}\\left(\\frac{\\sqrt{3} m}{2n + m}\\right)

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    rad2deg : bool, optional
        If `True`, return angle in degrees.

    Returns
    -------
    float
        Chiral angle :math:`\\theta_{c}` in
        degrees (default) or radians (if `rad2deg=False`).

    """
    n, m, _ = get_chiral_indices(*Ch)

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')
    theta = np.arctan(np.sqrt(3) * m / (2 * n + m))
    #return np.arccos((2*n + m) / (2 * np.sqrt(n**2 + m**2 + n*m)))
    if rad2deg:
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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

    if length:
        if bond is None:
            bond = CCbond

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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

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
            bond = CCbond

        return bond * np.sqrt(3 * (p**2 + q**2 + p * q))
    else:
        return (p, q)


def compute_R_chiral_angle(*Ch, rad2deg=True):
    """Compute "chiral angle" of symmetry vector :math:`\\theta_R`.

    .. math::

        \\theta_R = \\tan^{-1}\\left(\\frac{\\sqrt{3}q}{2p + q}\\right)

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    rad2deg : bool, optional
        If `True`, return angle in degrees

    Returns
    -------
    float
        Chiral angle of *symmetry vector* :math:`\\theta_R` in
        degrees (default) or radians (if `rad2deg=False`).


    """
    n, m, _ = get_chiral_indices(*Ch)

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

    p, q = compute_R(n, m)
    theta = np.arctan((np.sqrt(3) * q) / (2 * p + q))
    if rad2deg:
        return np.degrees(theta)
    else:
        return theta


def compute_symmetry_operation(*Ch, bond=None):
    """Compute symmetry operation :math:`(\\psi|\\tau)`.

    The *symmetry vector* `R` represents a *symmetry
    operation* of the nanotube which arises as a *screw translation*, which
    is a combination of a rotation :math:`\\psi` and translation
    :math:`\\tau`.
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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

    psi = compute_psi(n, m)
    tau = compute_tau(n, m, bond=bond)
    return (psi, tau)


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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

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
    str

    """
    n, m, _ = get_chiral_indices(*Ch)

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

    if (2 * n + m) % 3 == 1:
        return 'semiconducting, type 1'
    elif (2 * n + m) % 3 == 2:
        return 'semiconducting, type 2'
    else:
        return 'metallic'


def compute_Natoms(*Ch, nz=1):
    """Compute :math:`N_{\\mathrm{atoms/tube}}`

    .. math::

        N_{\\mathrm{atoms/tube}} = N_{\\mathrm{atoms/cell}} \\times
        n_{z-\\mathrm{cells}}

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nz : {int, float}
        Number of nanotube unit cells

    Returns
    -------
    int
        :math:`N_{\\mathrm{atoms/tube}}`
    """
    n, m, _ = get_chiral_indices(*Ch)

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

    Natoms_per_unit_cell = compute_Natoms_per_unit_cell(n, m)
    return int(Natoms_per_unit_cell * nz)


def compute_Natoms_per_tube(*Ch, nz=1):
    """Compute :math:`N_{\\mathrm{atoms/tube}}`

    .. math::

        N_{\\mathrm{atoms/tube}} = N_{\\mathrm{atoms/cell}} \\times
        n_{z-\\mathrm{cells}}

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nz : {int, float}
        Number of nanotube unit cells

    Returns
    -------
    int
        :math:`N_{\\mathrm{atoms/tube}}`
    """
    return compute_Natoms(*Ch, nz=nz)


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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

    N = compute_N(n, m)

    if element1 is None:
        element1 = 'C'
    if element2 is None:
        element2 = 'C'

    mass = N * (Atom(element1).m + Atom(element2).m)
    return mass


def compute_linear_mass_density(*Ch, bond=None, element1=None, element2=None,
                                **kwargs):
    """Compute nanotube linear mass density (mass per unit length) in \
    **grams/nm**.

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
        Linear mass density in units of **g/nm**.

    """
    n, m, _ = get_chiral_indices(*Ch)

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')

    mass = compute_unit_cell_mass(n, m, element1=element1, element2=element2,
                                  **kwargs)
    T = compute_T(n, m, bond=bond, length=True, **kwargs)

    try:
        linear_mass_density = mass / T
        # there are 1.6605e-24 grams / Da and 10 angstroms / nm
        linear_mass_density *= 10 * grams_per_Da
        return linear_mass_density
    except ZeroDivisionError:
        return 0


def compute_bundle_density(*Ch, d_vdw=None, bond=None,
                           element1=None, element2=None):
    """Compute nanotube bundle mass density \
    :math:`\\rho_{\\mathrm{bundle}}(n, m)` in :math:`\\mathrm{g/cm^3}`.

    .. math::

        \\rho_{\\mathrm{bundle}}(n, m) = \\frac{8\\pi^2 m_{\\mathrm{C}}
        \\sqrt{n^2 + m^2 + nm}}{9\\sqrt{3}a_{\\mathrm{CC}}^3 \\times
        \\left(\\sqrt{n^2 + m^2 + nm} +
        \\frac{\\pi d_{\\mathrm{vdW}}}{\\sqrt{3}a_{\\mathrm{CC}}}\\right)^2}

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    d_vdw : int
        van der Waals distance between nearest-neighbor tubes
    bond : float, optional
        Bond length.

    Returns
    -------
    float
        :math:`\\rho_{\\mathrm{bundle}}` in units of
        :math:`\\mathrm{\\frac{g}{cm^3}}`

    """
    n, m, _ = get_chiral_indices(*Ch)

    if bond is None:
        bond = CCbond

    if d_vdw is None:
        if n == m:
            d_vdw = 3.38
        elif (m == 0) or (n == 0):
            d_vdw = 3.41
        else:
            d_vdw = 3.39

    if element1 is None:
        element1 = 'C'
    if element2 is None:
        element2 = 'C'

    if element1 == element2:
        bundle_density = 8 * np.pi**2 * Atom(element1).m * \
            np.sqrt(n**2 + m**2 + n*m) / \
            (9 * np.sqrt(3) * bond**3 *
                (np.sqrt(n**2 + m**2 + n*m) +
                    np.pi * d_vdw / (np.sqrt(3) * bond))**2)
    else:
        bundle_density = 0

    # there are 1.6605e-24 grams / Da and 1e-8 cm / angstrom
    bundle_density *= grams_per_Da / (1e-8)**3
    return bundle_density


def compute_Lx(*Ch, nx=1, bond=None, gutter=dVDW):
    return nx * (compute_dt(*Ch, bond=bond) + gutter) / 10


def compute_Ly(*Ch, ny=1, bond=None, gutter=dVDW):
    return ny * (compute_dt(*Ch, bond=bond) + gutter) / 10


def compute_Lz(*Ch, nz=1, bond=None, **kwargs):
    """Compute :math:`L_z = L_{\\mathrm{tube}}` in **nanometers**.

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
        :math:`L_z = L_{\\mathrm{tube}}` in **nanometers**

    """
    n, m, _ = get_chiral_indices(*Ch)

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(nz, numbers.Real) or nz > 0):
        raise TypeError('Expected a real, positive number')

    T = compute_T(n, m, bond=bond, **kwargs)
    return nz * T / 10


def compute_symmetry_chiral_angle(*Ch, rad2deg=True):
    """Alias for :func:`compute_R_chiral_angle`."""
    return compute_R_chiral_angle(*Ch, rad2deg=rad2deg)


def compute_tube_diameter(*Ch, bond=None, **kwargs):
    """Alias for :func:`compute_dt`"""
    return compute_dt(*Ch, bond=bond, **kwargs)


def compute_tube_radius(*Ch, bond=None, **kwargs):
    """Alias for :func:`compute_rt`"""
    return compute_rt(*Ch, bond=bond, **kwargs)


def compute_tube_length(*Ch, nz=1, bond=None, **kwargs):
    """Alias for :func:`compute_Lz`"""
    return compute_Lz(*Ch, nz=nz, bond=bond, **kwargs)


def compute_tube_mass(*Ch, nz=1, element1=None, element2=None):
    """Compute nanotube mass in **grams**.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nz : {int, float}
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

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(nz, numbers.Real) or nz > 0):
        raise TypeError('Expected a real, positive number')

    Natoms = compute_Natoms(n, m, nz=nz)

    if element1 is None:
        element1 = 'C'
    if element2 is None:
        element2 = 'C'

    atom1 = Atom(element1)
    atom2 = Atom(element2)

    mass = Natoms * (atom1.m + atom2.m) / 2
    # there are 1.6605e-24 grams / Da
    mass *= grams_per_Da
    return mass
