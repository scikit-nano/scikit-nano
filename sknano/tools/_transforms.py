# -*- coding: utf-8 -*-
"""
===============================================================================
Linear algebra functions for transformations (:mod:`sknano.tools._transforms`)
===============================================================================

.. currentmodule:: sknano.tools._transforms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numpy as np

from ._luts import xyz_axes

__all__ = ['rotation_matrix', 'transformation_matrix', 'rotate_point']

I = np.identity(4)

_str2array = {}
for i, axis in enumerate(xyz_axes):
    _str2array[axis] = np.asarray(I[:3, i])


def Rx(angle=None, deg2rad=False):
    pass


def Ry(angle=None, deg2rad=False):
    pass


def Rz(angle=None, deg2rad=False):
    pass


def rotation_matrix(angle=None, rot_axis=None, deg2rad=False, R4x4=False):
    """Generate a :math:`3\\times 3` or :math:`4\\times 4` rotation matrix.

    Parameters
    ----------
    angle : float
        Rotation angle in **radians** unless `deg2rad` is `True`.
        The sense of the rotation is defined by the *right hand rule*:
        If you point your thumb alo
    rot_axis : {array_like, str}
        A 3-element list or ndarray defining the 3 components,
        :math:`u, v, w`, of the vector defining the axis of rotation.
    deg2rad : bool, optional
        Angle is in degrees and needs to be converted to radians
    R4x4 : bool, optional
        If `True`, return a :math:`4\\times 4` matrix.
        If `False`, return a :math:`3\\times 3` matrix.

    Returns
    -------
    Rmat : ndarray
        rotation matrix

    """
    if angle is None or rot_axis is None:
        raise TypeError('`angle` and `rot_axis` are required')
    if isinstance(rot_axis, (str, unicode)):
        try:
            rot_axis = _str2array[rot_axis]
        except KeyError:
            raise ValueError('Invalid `rot_axis` string: {}'.format(rot_axis))
    elif isinstance(rot_axis, (list, np.ndarray)) and len(rot_axis) != 3:
        raise ValueError('`rot_axis` must be a 3-element list or ndarray')

    if deg2rad:
        angle = np.radians(angle)

    cosa = np.cos(angle)
    sina = np.sin(angle)

    Rmat = None

    if np.allclose(rot_axis, I[:3, 0]):
        Rmat = np.array([[1, 0, 0, 0],
                         [0, cosa, -sina, 0],
                         [0, sina, cosa, 0],
                         [0, 0, 0, 1]])
    elif np.allclose(rot_axis, I[:3, 1]):
        Rmat = np.array([[cosa, 0, sina, 0],
                         [0, 1, 0, 0],
                         [-sina, 0, cosa, 0],
                         [0, 0, 0, 1]])
    elif np.allclose(rot_axis, I[:3, 2]):
        Rmat = np.array([[cosa, -sina, 0, 0],
                         [sina, cosa, 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])
    else:
        u, v, w = rot_axis
        uu, vv, ww = u**2, v**2, w**2
        uv, uw, vw = u * v, u * w, v * w
        l = np.sqrt(uu + vv + ww)
        ll = l**2

        r11 = (uu + (vv + ww) * cosa) / ll
        r22 = (vv + (uu + ww) * cosa) / ll
        r33 = (ww + (uu + vv) * cosa) / ll

        r12 = (uv * (1 - cosa) - w * l * sina) / ll
        r21 = (uv * (1 - cosa) + w * l * sina) / ll

        r13 = (uw * (1 - cosa) + v * l * sina) / ll
        r31 = (uw * (1 - cosa) - v * l * sina) / ll

        r23 = (vw * (1 - cosa) - u * l * sina) / ll
        r32 = (vw * (1 - cosa) + u * l * sina) / ll

        Rmat = np.array([[r11, r12, r13, 0],
                         [r21, r22, r23, 0],
                         [r31, r32, r33, 0],
                         [0, 0, 0, 1]])

        Rmat[np.where(np.abs(Rmat) <= np.finfo(float).eps)] = 0.0

    if not R4x4:
        Rmat = Rmat[:3,:3]

    return Rmat


def rotate_point(r0=None, angle=None, p0=None, p1=None, deg2rad=False):
    """Rotate point about arbitrary axis.

    Parameters
    ----------
    axis : `Vector`

    Returns
    -------
    ndarray
        3-element ndarray of (:math:`x,y,z`) coordinates of rotated
        point.

    """
    if deg2rad:
        angle = np.radians(angle)

    a, b, c = p0
    d, e, f = p1
    u, v, w = d - a, e - b, f - c

    Tmat = transformation_matrix(angle=angle, P0=[a, b, c],
                                 rot_axis=[u, v, w], deg2rad=deg2rad)

    x, y, z = r0
    t0 = np.array([x, y, z, 1])
    t = np.dot(Tmat, t0)
    r = t[:3]

    return r


def transformation_matrix(angle=None, rot_axis=None, rvec=None,
                          P0=None, P1=None, deg2rad=False, T4x4=True):
    """Generate a :math:`3\\times 3` or :math:`4\\times 4` rotation matrix.

    Parameters
    ----------
    angle : float
        Rotation angle in **radians** unless `deg2rad` is `True`.
        The sense of the rotation is defined by the *right hand rule*:
        If you point your thumb alo
    rot_axis : {array_like, str}
        A 3-element list or ndarray defining the 3 components,
        :math:`u, v, w`, of the vector defining the axis of rotation.
    deg2rad : bool, optional
        Angle is in degrees and needs to be converted to radians
    T4x4 : bool, optional
        If `True`, return a :math:`4\\times 4` matrix.
        If `False`, return a :math:`3\\times 3` matrix.

    Returns
    -------
    Tmat : ndarray
        Transformation matrix

    """
    if angle is None or P0 is None or (rot_axis is None and P1 is None):
        raise TypeError('`angle`, `P0`, and `rot_axis` or `P1` are required')
    if isinstance(rot_axis, (str, unicode)):
        try:
            rot_axis = _str2array[rot_axis]
        except KeyError:
            raise ValueError('Invalid `rot_axis` string: {}'.format(rot_axis))
    elif isinstance(rot_axis, (list, np.ndarray)) and len(rot_axis) != 3:
        raise ValueError('`rot_axis` must be a 3-element list or ndarray')
    elif rot_axis is None and isinstance(P1, (list, np.ndarray)) and \
            len(P1) != 3:
        raise ValueError('`P1` must be a 3-element list or ndarray')

    if rot_axis is None:
        rot_axis = np.asarray(P1) - np.asarray(P0)

    Tmat = rotation_matrix(angle=angle, rot_axis=rot_axis, deg2rad=deg2rad,
                           R4x4=True)

    if deg2rad:
        angle = np.radians(angle)

    cosa = np.cos(angle)
    sina = np.sin(angle)

    a, b, c = P0
    u, v, w = rot_axis
    uu, vv, ww = u**2, v**2, w**2
    au, bu, cu = a * u, b * u, c * u
    av, bv, cv = a * v, b * v, c * v
    aw, bw, cw = a * w, b * w, c * w
    l = np.sqrt(uu + vv + ww)
    ll = l**2

    t14 = ((a * (vv + ww) - u * (bv + cw)) * (1 - cosa) +
           (bw - cv) * l * sina) / ll
    t24 = ((b * (uu + ww) - v * (au + cw)) * (1 - cosa) +
           (cu - aw) * l * sina) / ll
    t34 = ((c * (uu + vv) - w * (au + bv)) * (1 - cosa) +
           (av - bu) * l * sina) / ll

    Tmat[:3, 3] = np.array([t14, t24, t34])

    Tmat[np.where(np.abs(Tmat) <= np.finfo(float).eps)] = 0.0

    if not T4x4:
        Tmat = Tmat[:3,:3]

    return Tmat
