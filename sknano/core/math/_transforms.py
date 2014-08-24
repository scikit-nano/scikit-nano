# -*- coding: utf-8 -*-
"""
===============================================================================
Functions for linear algebra transforms (:mod:`sknano.core.math._transforms`)
===============================================================================

.. currentmodule:: sknano.core.math._transforms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numpy as np
from ._point import Point
from ._vector import Vector

__all__ = ['Rx', 'Ry', 'Rz',
           'rotation_matrix',
           'transformation_matrix',
           'rotation_transform',
           'rotate_point']

I = np.identity(4)

_str2array = {}
for i, axis in enumerate(('x', 'y', 'z')):
    _str2array[axis] = Vector(np.asarray(I[:3, i]))


def Rx(angle, deg2rad=False):
    if deg2rad:
        angle = np.radians(angle)
    cosa = np.cos(angle)
    sina = np.sin(angle)

    return np.array([[1.0, 0.0, 0.0],
                     [0.0, cosa, -sina],
                     [0.0, sina, cosa]])


def Ry(angle, deg2rad=False):
    if deg2rad:
        angle = np.radians(angle)
    cosa = np.cos(angle)
    sina = np.sin(angle)

    return np.array([[cosa, 0.0, sina],
                     [0.0, 1.0, 0.0],
                     [-sina, 0.0, cosa]])


def Rz(angle, deg2rad=False):
    if deg2rad:
        angle = np.radians(angle)
    cosa = np.cos(angle)
    sina = np.sin(angle)

    return np.array([[cosa, -sina, 0.0],
                     [sina, cosa, 0.0],
                     [0.0, 0.0, 1.0]])


def rotation_matrix(angle=None, rot_axis=None, deg2rad=False):
    """Generate an :math:`n\\times n` rotation matrix.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    angle : float
        Rotation angle in **radians**. If `deg2rad` is `True`, `angle` will be
        converted to radians from degrees.  The *sense* of the rotation is
        defined by the *right hand rule*: If your right-hand's thumb points
        along the `rot_axis`, then your fingers wrap around the axis in the
        *positive sense* of the rotation angle.
    rot_axis : {None, array_like, str}, optional
        An :math:`n`-element array_like sequence defining the :math:`n`
        components of the rotation axis or the string `x`, `y`, or `z`
        representing the :math:`x, y, z` axes of a Cartesian coordinate
        system in 3D with unit vectors
        :math:`\\mathbf{v}_x=\\mathbf{\\hat{x}}`,
        :math:`\\mathbf{v}_y=\\mathbf{\\hat{y}}`, and
        :math:`\\mathbf{v}_z=\\mathbf{\\hat{z}}`, respectively.
    deg2rad : bool, optional
        If `True`, convert `angle` from degrees to radians.

    Returns
    -------
    Rmat : ndarray
        If `rot_axis` is `None` then `Rmat` will be a :math:`2D`
        rotation matrix :math:`R(\\theta)` that rotates :math:`2D` vectors
        counterclockwise by `angle` :math:`\\theta`.

        If `rot_axis` is not `None` then `Rmat` will be a rotation matrix that
        gives a rotation around the direction of the vector `rot_axis`.

    """
    Rmat = transformation_matrix(angle, rot_axis=rot_axis, deg2rad=deg2rad)
    return Rmat[:-1,:-1]


def transformation_matrix(angle, rot_axis=None, axis_origin=None,
                          deg2rad=False):
    """Generate an :math:`n+1\\times n+1` transformation matrix for
    an affine transformation in :math:`n` dimensions.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    angle : float
        Rotation angle in **radians**. If `deg2rad` is `True`, `angle` will be
        converted to radians from degrees.  The *sense* of the rotation is
        defined by the *right hand rule*: If your right-hand's thumb points
        along the `rot_axis`, then your fingers wrap around the axis in the
        *positive sense* of the rotation angle.
    rot_axis : {None, array_like, str}, optional
        An :math:`n`-element array_like sequence defining the :math:`n`
        components of the rotation axis or the string `x`, `y`, or `z`
        representing the :math:`x, y, z` axes of a Cartesian coordinate
        system in 3D with unit vectors
        :math:`\\mathbf{v}_x=\\mathbf{\\hat{x}}`,
        :math:`\\mathbf{v}_y=\\mathbf{\\hat{y}}`, and
        :math:`\\mathbf{v}_z=\\mathbf{\\hat{z}}`, respectively.
    axis_origin : {None, array_like}, optional
        An :math:`n`-element list or ndarray or
        :class:`~sknano.core.math.Point` defining
        the origin of the rotation axis.

        If `axis_origin` is not `None` and `rot_axis` is a `Vector` instance,
        then the origin of the vector defined by :attr:`Vector.p0` will be
        changed to `axis_origin`.

        If `axis_origin` is `None`, then it defaults to an
        :math:`n`-element array of zeros.
    deg2rad : bool, optional
        If `True`, convert `angle` from degrees to radians.

    Returns
    -------
    Tmat : ndarray
        :math:`n+1\\times n+1` transformation matrix for an affine transform
        in :math:`n` dimensions.

        If `rot_axis` is `None` and `axis_origin` is `None`,
        then `Tmat` will be a :math:`2D` rotation matrix :math:`R(\\theta)`
        that rotates :math:`2D` vectors counterclockwise by `angle`
        :math:`\\theta`.

        If `rot_axis` is `None` and `axis_origin` is a 2-element sequence,
        then `Rmat` will be a :math:`2D` rotation matrix :math:`R(\\theta)`
        about the :math:`2D` `Point` `axis_origin` by `angle`
        :math:`\\theta`.

        If `rot_axis` is not `None` and `axis_origin` is `None`,
        then `Rmat` will be a rotation matrix that gives a rotation around
        the direction of the vector `rot_axis`.

    Notes
    -----

    """
    if deg2rad:
        angle = np.radians(angle)

    cosa = np.cos(angle)
    sina = np.sin(angle)

    # Handle 2D rotation about origin
    if rot_axis is None and axis_origin is None:
        return Rz(angle)

    # Handle 3D rotation about origin
    # Handle 3D rotation around the rotation vector
    # Handle N-D rotation about origin
    # Handle N-D rotation around the N-D rotation vector anchored at
    # an arbitrary N-D point.
    if rot_axis is not None and isinstance(rot_axis, (str, unicode)):
        try:
            rot_axis = _str2array[rot_axis]
        except KeyError:
            raise ValueError('Invalid `rot_axis` string: {}'.format(rot_axis))
    elif not isinstance(rot_axis, (tuple, list, np.ndarray)):
        raise ValueError('`rot_axis` must be a sequence')

    if axis_origin is not None and \
            isinstance(axis_origin, (tuple, list, np.ndarray)) and \
            len(axis_origin) == len(rot_axis):
        if isinstance(rot_axis, Vector):
            rot_axis.p0 = Point(axis_origin)
        else:
            rot_axis = Vector(rot_axis, p0=axis_origin)

    if not isinstance(rot_axis, Vector):
        rot_axis = Vector(rot_axis)

    if rot_axis.nd == 2:
        Tmat = Rz(angle)
        # Handle 2D rotation about arbitrary 2D point
        if not np.allclose(rot_axis, np.zeros(2)):
            a, b = rot_axis
            m13 = a - a * cosa + b * sina
            m23 = b - b * cosa - a * sina
            Tmat[:, 2] = np.array([m13, m23])
    else:
        Tmat = np.zeros((4, 4))

        if np.allclose(rot_axis, I[:3, 0]):
            Tmat[:3,:3] = Rx(angle)
        elif np.allclose(rot_axis, I[:3, 1]):
            Tmat[:3,:3] = Ry(angle)
        elif np.allclose(rot_axis, I[:3, 2]):
            Tmat[:3,:3] = Rz(angle)
        else:
            u, v, w = rot_axis
            uu, vv, ww = u**2, v**2, w**2
            uv, uw, vw = u * v, u * w, v * w
            l = np.sqrt(uu + vv + ww)
            ll = l**2

            m11 = (uu + (vv + ww) * cosa) / ll
            m22 = (vv + (uu + ww) * cosa) / ll
            m33 = (ww + (uu + vv) * cosa) / ll
            m12 = (uv * (1 - cosa) - w * l * sina) / ll
            m21 = (uv * (1 - cosa) + w * l * sina) / ll
            m13 = (uw * (1 - cosa) + v * l * sina) / ll
            m31 = (uw * (1 - cosa) - v * l * sina) / ll
            m23 = (vw * (1 - cosa) - u * l * sina) / ll
            m32 = (vw * (1 - cosa) + u * l * sina) / ll

            Tmat[:3,:3] = np.array([[m11, m12, m13],
                                    [m21, m22, m23],
                                    [m31, m32, m33]])

        if not np.allclose(rot_axis.p0, np.zeros(3)):
            a, b, c = rot_axis.p0
            u, v, w = rot_axis
            uu, vv, ww = u**2, v**2, w**2
            au, bu, cu = a * u, b * u, c * u
            av, bv, cv = a * v, b * v, c * v
            aw, bw, cw = a * w, b * w, c * w
            l = np.sqrt(uu + vv + ww)
            ll = l**2

            m14 = ((a * (vv + ww) - u * (bv + cw)) * (1 - cosa) +
                   (bw - cv) * l * sina) / ll
            m24 = ((b * (uu + ww) - v * (au + cw)) * (1 - cosa) +
                   (cu - aw) * l * sina) / ll
            m34 = ((c * (uu + vv) - w * (au + bv)) * (1 - cosa) +
                   (av - bu) * l * sina) / ll

            Tmat[:, 3] = np.array([m14, m24, m34, 1.0])

    Tmat[np.where(np.abs(Tmat) <= np.finfo(float).eps)] = 0.0

    return Tmat


def rotation_transform(angle, rot_axis=None, axis_origin=None, deg2rad=False):
    """Generate an :math:`n\\times n` rotation matrix.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    angle : float
        Rotation angle in **radians** unless `deg2rad` is `True`.
        The *sense* of the rotation is defined by the *right hand rule*:
        If your right-hand's thumb points along the `rot_axis`,
        then your fingers wrap around the axis in the *positive sense* of
        the rotation angle.
    rot_axis : {array_like, str, :class:`~sknano.core.math.Vector`}
        3-element list or ndarray or :class:`~sknano.core.math.Vector` defining
        the 3 components, :math:`u, v, w`, of the vector defining the axis
        of rotation.
    axis_origin : {array_like, :class:`~sknano.core.math.Point`}
        3-element list or ndarray or :class:`~sknano.core.math.Point` defining
        origin point of rotation axis.
    deg2rad : bool, optional
        Angle is in degrees and needs to be converted to radians

    Returns
    -------
    transformation_function

    """
    Tmat = transformation_matrix(angle, rot_axis=rot_axis,
                                 axis_origin=axis_origin,
                                 deg2rad=deg2rad)

    def transform_function(p0):
        p1 = np.ones(len(p0) + 1)
        p1[:-1] = p0
        return Point(np.dot(Tmat, p1)[:3])

    return transform_function


def rotate_point(point, angle, rot_axis=None, axis_origin=None,
                 deg2rad=False):
    """Rotate point about arbitrary axis.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    point : array_like or :class:`~sknano.core.math.Point`
    angle : float
    rot_axis : {array_like, str, :class:`~sknano.core.math.Vector`}
    axis_origin : {array_like, :class:`~sknano.core.math.Point`}
    deg2rad : bool, optional

    Returns
    -------
    rotated_point : ndarray

    """
    rotfn = rotation_transform(angle, rot_axis=rot_axis,
                               axis_origin=axis_origin,
                               deg2rad=deg2rad)

    return rotfn(point)
