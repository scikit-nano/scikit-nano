# -*- coding: utf-8 -*-
"""
===============================================================================
Linear algebra functions for transformations (:mod:`sknano.core._transforms`)
===============================================================================

.. currentmodule:: sknano.core._transforms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numpy as np

__all__ = ['rotate_point', 'rotation_matrix', 'transformation_matrix']

I = np.identity(4)

_str2array = {}
for i, axis in enumerate(('x', 'y', 'z')):
    _str2array[axis] = np.asarray(I[:3, i])


def Rx(angle=None, deg2rad=False):
    pass


def Ry(angle=None, deg2rad=False):
    pass


def Rz(angle=None, deg2rad=False):
    pass


def rotate_point(point=None, angle=None, rot_axis=None, axis_origin=None,
                 deg2rad=False):
    """Rotate point about arbitrary axis.

    .. versionadded:: 0.2.26

    Parameters
    ----------
    point : array_like or :class:`~sknano.core.Point`
        3-element list or ndarray or :class:`~sknano.core.Point`
        defining the :math:`x, y, z` coordinates of point to rotate.
    angle : float
        Rotation angle in **radians** unless `deg2rad` is `True`.
        The *sense* of the rotation is defined by the *right hand rule*:
        If your right-hand's thumb points along the `rot_axis`,
        then your fingers wrap around the axis in the *positive sense* of
        the rotation angle.
    rot_axis : {array_like, str, :class:`~sknano.core.Vector`}
        3-element list or ndarray or :class:`~sknano.core.Vector` defining
        the 3 components, :math:`u, v, w`, of the vector defining the axis
        of rotation.
    axis_origin : {array_like, :class:`~sknano.core.Point`}
        3-element list or ndarray or :class:`~sknano.core.Point` defining
        origin point of rotation axis.
    deg2rad : bool, optional
        Angle is in degrees and needs to be converted to radians

    Returns
    -------
    rotated_point : ndarray
        3-element ndarray of (:math:`x,y,z`) coordinates of rotated
        point.

    """
    if deg2rad:
        angle = np.radians(angle)

    Tmat = transformation_matrix(angle=angle, rot_axis=rot_axis,
                                 axis_origin=axis_origin,
                                 deg2rad=deg2rad)

    if point is None:
        point = np.zeros(3)

    x, y, z = point
    point = np.array([x, y, z, 1])
    rotated_transformation = np.dot(Tmat, point)
    rotated_point = rotated_transformation[:3]

    return rotated_point


def rotation_matrix(angle=None, rot_axis=None, deg2rad=False, R4x4=False):
    """Generate a :math:`3\\times 3` or :math:`4\\times 4` rotation matrix.

    Parameters
    ----------
    angle : float
        Rotation angle in **radians** unless `deg2rad` is `True`.
        The *sense* of the rotation is defined by the *right hand rule*:
        If your right-hand's thumb points along the `rot_axis`,
        then your fingers wrap around the axis in the *positive sense* of
        the rotation angle.
    rot_axis : {array_like, str, :class:`~sknano.core.Vector`}
        3-element list or ndarray or :class:`~sknano.core.Vector` defining
        the 3 components, :math:`u, v, w`, of the vector defining the axis
        of rotation.
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
    Rmat = transformation_matrix(angle=angle, rot_axis=rot_axis,
                                 deg2rad=deg2rad)

    if not R4x4:
        Rmat = Rmat[:3,:3]

    return Rmat


def transformation_matrix(angle=None, rot_axis=None, axis_origin=None,
                          deg2rad=False):
    """Generate a :math:`4\\times 4` transformation matrix.

    .. versionadded:: 0.2.26

    Parameters
    ----------
    angle : float
        Rotation angle in **radians** unless `deg2rad` is `True`.
        The *sense* of the rotation is defined by the *right hand rule*:
        If your right-hand's thumb points along the `rot_axis`,
        then your fingers wrap around the axis in the *positive sense* of
        the rotation angle.
    rot_axis : {array_like, str, :class:`~sknano.core.Vector`}
        3-element list or ndarray or :class:`~sknano.core.Vector` defining
        the 3 components, :math:`u, v, w`, of the vector defining the axis
        of rotation.
    axis_origin : {array_like, :class:`~sknano.core.Point`}
        3-element list or ndarray or :class:`~sknano.core.Point` defining
        origin point of axis of rotation
    deg2rad : bool, optional
        Angle is in degrees and needs to be converted to radians

    Returns
    -------
    Tmat : ndarray
        :math:`4\\times 4` transformation matrix

    """
    if angle is None or rot_axis is None:
        raise TypeError('`angle` and `rot_axis` are required')
    if isinstance(rot_axis, (str, unicode)):
        try:
            rot_axis = _str2array[rot_axis]
        except KeyError:
            raise ValueError('Invalid `rot_axis` string: {}'.format(rot_axis))
    elif not (isinstance(rot_axis, (list, np.ndarray)) and len(rot_axis) == 3):
        raise ValueError('`rot_axis` must be a 3-element list or ndarray or '
                         '`Vector`')

    if axis_origin is None:
        axis_origin = np.zeros(3, dtype=float)
    elif not (isinstance(axis_origin, (list, np.ndarray)) and
              len(axis_origin) == 3):
        raise ValueError('`axis_origin` must be a 3-element list or ndarray '
                         'or `Point`')

    if deg2rad:
        angle = np.radians(angle)

    cosa = np.cos(angle)
    sina = np.sin(angle)

    Tmat = None

    if np.allclose(rot_axis, I[:3, 0]):
        Tmat = np.array([[1, 0, 0, 0],
                         [0, cosa, -sina, 0],
                         [0, sina, cosa, 0],
                         [0, 0, 0, 1]])
    elif np.allclose(rot_axis, I[:3, 1]):
        Tmat = np.array([[cosa, 0, sina, 0],
                         [0, 1, 0, 0],
                         [-sina, 0, cosa, 0],
                         [0, 0, 0, 1]])
    elif np.allclose(rot_axis, I[:3, 2]):
        Tmat = np.array([[cosa, -sina, 0, 0],
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

        Tmat = np.array([[r11, r12, r13, 0],
                         [r21, r22, r23, 0],
                         [r31, r32, r33, 0],
                         [0, 0, 0, 1]])

    if not np.allclose(axis_origin, np.zeros(3)):
        a, b, c = axis_origin
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

    return Tmat
