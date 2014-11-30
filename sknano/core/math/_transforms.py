# -*- coding: utf-8 -*-
"""
===============================================================================
Functions for linear algebra transforms (:mod:`sknano.core.math._transforms`)
===============================================================================

.. currentmodule:: sknano.core.math._transforms

"""
from __future__ import absolute_import, division, print_function
import six
__docformat__ = 'restructuredtext en'

import numpy as np

__all__ = ['transformation_matrix',
           'affine_transform',
           'reflection_transform',
           'rotation_transform',
           'scaling_transform',
           'Rx', 'Ry', 'Rz',
           'reflection_matrix',
           'rotation_matrix',
           'scaling_matrix',
           'rotate', 'scale', 'translate']

I = np.identity(4)

_str2array = {}
for i, axis in enumerate(('x', 'y', 'z')):
    _str2array[axis] = I[:3, i]


def transformation_matrix(angle=None, rot_axis=None, anchor_point=None,
                          rot_point=None, from_vector=None, to_vector=None,
                          deg2rad=False, verbose=False):
    """Generate an :math:`(n+1)\\times(n+1)` transformation matrix for an \
        affine transformation in :math:`n` dimensions.

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
    anchor_point : {None, array_like}, optional
        An :math:`n`-element list or ndarray or
        :class:`~sknano.core.math.Point` defining
        the origin of the rotation axis.

        If `anchor_point` is not `None` and `rot_axis` is a `Vector` instance,
        then the origin of the vector defined by :attr:`Vector.p0` will be
        changed to `anchor_point`.

        If `anchor_point` is `None`, then it defaults to an
        :math:`n`-element array of zeros.
    deg2rad : bool, optional
        If `True`, convert `angle` from degrees to radians.

    Returns
    -------
    Tmat : ndarray
        :math:`n+1\\times n+1` transformation matrix for an affine transform
        in :math:`n` dimensions.

        If `rot_axis` is `None` and `anchor_point` is `None`,
        then `Tmat` will be a :math:`2D` rotation matrix :math:`R(\\theta)`
        that rotates :math:`2D` vectors counterclockwise by `angle`
        :math:`\\theta`.

        If `rot_axis` is `None` and `anchor_point` is a 2-element sequence,
        then `Rmat` will be a :math:`2D` rotation matrix :math:`R(\\theta)`
        about the :math:`2D` `Point` `anchor_point` by `angle`
        :math:`\\theta`.

        If `rot_axis` is not `None` and `anchor_point` is `None`,
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
    if rot_axis is None and anchor_point is None:
        return Rz(angle)

    # Handle 3D rotation about origin
    # Handle 3D rotation around the rotation vector
    # Handle N-D rotation about origin
    # Handle N-D rotation around the N-D rotation vector anchored at
    # an arbitrary N-D point.
    if rot_axis is not None and isinstance(rot_axis, (str, six.text_type)):
        try:
            rot_axis = _str2array[rot_axis]
        except KeyError:
            raise ValueError(
                'Invalid `rot_axis` string: {}'.format(rot_axis))
    elif not isinstance(rot_axis, (tuple, list, np.ndarray)):
        raise ValueError('`rot_axis` must be a sequence')

    if anchor_point is None:
        anchor_point = np.zeros(len(rot_axis))

    nd = len(rot_axis)
    if nd == 2:
        Tmat = Rz(angle)
        # Handle 2D rotation about arbitrary 2D point
        if not np.allclose(rot_axis, np.zeros(2)):
            a, b = rot_axis
            m13 = a - a * cosa + b * sina
            m23 = b - b * cosa - a * sina
            Tmat[:2, 2] = np.array([m13, m23])
    elif nd == 3:
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

        if not np.allclose(anchor_point, np.zeros(3)):
            a, b, c = anchor_point
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


def affine_transform(angle=None, rot_axis=None, anchor_point=None,
                     rot_point=None, from_vector=None, to_vector=None,
                     deg2rad=False, verbose=False):
    raise NotImplementedError("Not implemented")


def reflection_transform():
    raise NotImplementedError("Not implemented")


def rotation_transform(angle=None, rot_axis=None, anchor_point=None,
                       rot_point=None, from_vector=None, to_vector=None,
                       deg2rad=False, verbose=False):
    raise NotImplementedError("Not implemented")


def scaling_transform():
    raise NotImplementedError("Not implemented")


def Rx(angle, deg2rad=False):
    """Generate the :math:`3\\times3` rotation matrix :math:`R_x(\\theta)` \
        for a rotation about the :math:`x` axis by an angle :math:`\\theta`.

    Parameters
    ----------
    angle : float
        The rotation angle :math:`\\theta` in *radians*. If the angle is
        given in *degrees*, then you must set `deg2rad=True` to correctly
        calculate the rotation matrix.
    deg2rad : bool, optional
        if `True`, then `angle` is converted from degrees to radians.

    Returns
    -------
    :class:`~numpy:numpy.ndarray`
        :math:`3\\times3` rotation matrix :math:`R_x(\\theta)` for a
        rotation about the :math:`x` axis by an angle :math:`\\theta`.

    .. math::

       R_x = \\begin{pmatrix}
       1 & 0 & 0\\\\
       0 & \\cos\\theta & -\\sin\\theta\\\\
       0 & \\sin\\theta & \\cos\\theta
       \\end{pmatrix}

    Examples
    --------
    >>> import numpy as np
    >>> from sknano.core.math import Rx
    >>> Rx(np.pi/4)
    array([[ 1.        ,  0.        ,  0.        ],
           [ 0.        ,  0.70710678, -0.70710678],
           [ 0.        ,  0.70710678,  0.70710678]])
    >>> np.alltrue(Rx(np.pi/4) == Rx(45, deg2rad=True))
    True

    """
    if deg2rad:
        angle = np.radians(angle)
    cosa = np.cos(angle)
    sina = np.sin(angle)

    Rmat = np.array([[1.0, 0.0, 0.0], [0.0, cosa, -sina], [0.0, sina, cosa]])
    Rmat[np.where(np.abs(Rmat) <= np.finfo(float).eps)] = 0.0
    return Rmat


def Ry(angle, deg2rad=False):
    """Generate the :math:`3\\times3` rotation matrix :math:`R_y(\\theta)` \
        for a rotation about the :math:`y` axis by an angle :math:`\\theta`.

    Parameters
    ----------
    angle : float
        The rotation angle :math:`\\theta` in *radians*. If the angle is
        given in *degrees*, then you must set `deg2rad=True` to correctly
        calculate the rotation matrix.
    deg2rad : bool, optional
        if `True`, then `angle` is converted from degrees to radians.

    Returns
    -------
    :class:`~numpy:numpy.ndarray`
        :math:`3\\times3` rotation matrix :math:`R_y(\\theta)` for a
        rotation about the :math:`y` axis by an angle :math:`\\theta`:

    .. math::

       R_y = \\begin{pmatrix}
       \\cos\\theta & 0 & \\sin\\theta\\\\
       0 & 1 & 0\\\\
       -\\sin\\theta & 0 & \\cos\\theta
       \\end{pmatrix}

    Examples
    --------
    >>> import numpy as np
    >>> from sknano.core.math import Ry
    >>> Ry(np.pi/4)
    array([[ 0.70710678,  0.        ,  0.70710678],
           [ 0.        ,  1.        ,  0.        ],
           [-0.70710678,  0.        ,  0.70710678]])
    >>> np.alltrue(Ry(np.pi/4) == Ry(45, deg2rad=True))
    True

    """
    if deg2rad:
        angle = np.radians(angle)
    cosa = np.cos(angle)
    sina = np.sin(angle)

    Rmat = np.array([[cosa, 0.0, sina], [0.0, 1.0, 0.0], [-sina, 0.0, cosa]])
    Rmat[np.where(np.abs(Rmat) <= np.finfo(float).eps)] = 0.0
    return Rmat


def Rz(angle, deg2rad=False):
    """Generate the :math:`3\\times3` rotation matrix :math:`R_z(\\theta)` \
        for a rotation about the :math:`z` axis by an angle :math:`\\theta`.

    Parameters
    ----------
    angle : float
        The rotation angle :math:`\\theta` in *radians*. If the angle is
        given in *degrees*, then you must set `deg2rad=True` to correctly
        calculate the rotation matrix.
    deg2rad : bool, optional
        if `True`, then `angle` is converted from degrees to radians.

    Returns
    -------
    :class:`~numpy:numpy.ndarray`
        :math:`3\\times3` rotation matrix :math:`R_z(\\theta)` for a
        rotation about the :math:`z` axis by an angle :math:`\\theta`.

    .. math::

       R_z = \\begin{pmatrix}
       \\cos\\theta & -\\sin\\theta & 0\\\\
       \\sin\\theta & \\cos\\theta & 0\\\\
       0 & 0 & 1
       \\end{pmatrix}

    Examples
    --------
    >>> import numpy as np
    >>> from sknano.core.math import Rz
    >>> Rz(np.pi/4)
    array([[ 0.70710678, -0.70710678,  0.        ],
           [ 0.70710678,  0.70710678,  0.        ],
           [ 0.        ,  0.        ,  1.        ]])
    >>> np.alltrue(Rz(np.pi/4) == Rz(45, deg2rad=True))
    True

    """
    if deg2rad:
        angle = np.radians(angle)
    cosa = np.cos(angle)
    sina = np.sin(angle)

    Rmat = np.array([[cosa, -sina, 0.0], [sina, cosa, 0.0], [0.0, 0.0, 1.0]])
    Rmat[np.where(np.abs(Rmat) <= np.finfo(float).eps)] = 0.0
    return Rmat


def reflection_matrix():
    pass


def rotation_matrix(angle=None, rot_axis=None, anchor_point=None,
                    rot_point=None, from_vector=None, to_vector=None,
                    deg2rad=False, verbose=False):
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
    anchor_point : :class:`~sknano.core.math.Point`, optional
    rot_point : :class:`~sknano.core.math.Point`, optional
    from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
    deg2rad : bool, optional
        If `True`, convert `angle` from degrees to radians.

    Returns
    -------
    Rmat : :class:`~numpy:numpy.ndarray`
        If `rot_axis` is `None` then `Rmat` will be a :math:`2D`
        rotation matrix :math:`R(\\theta)` that rotates :math:`2D` vectors
        counterclockwise by `angle` :math:`\\theta`.

        If `rot_axis` is not `None` then `Rmat` will be a rotation matrix that
        gives a rotation around the direction of the vector `rot_axis`.

    """
    Rmat = transformation_matrix(angle=angle, rot_axis=rot_axis,
                                 anchor_point=anchor_point, rot_point=rot_point,
                                 from_vector=from_vector, to_vector=to_vector,
                                 deg2rad=deg2rad, verbose=verbose)
    return Rmat[:-1,:-1]


def scaling_matrix(s=None, v=None):
    """Return scaling matrix.

    Parameters
    ----------
    s : {list, float}
    v : {None, :class:`~sknano.core.math.Vector`}, optional

    Returns
    -------
    :class:`~numpy:numpy.ndarray`

    """
    if not isinstance(s, (list, float)):
        raise TypeError('Expected `s` to be a list or a float.')
    if isinstance(s, float) and not isinstance(v, (list, np.ndarray)):
        raise TypeError('Expected `v` to be a list or numpy array')

    if isinstance(s, list):
        return np.diag(s)
    else:
        from ._vector import Vector
        v = Vector(v)

        Smat = np.zeros((v.nd, v.nd))

        if v.nd == 2:
            s11 = s * v.x**2 + v.y**2
            s22 = v.x**2 + s * v.y**2
            s12 = s21 = (s - 1) * v.x * v.y
            Smat[:, :] = 1 / v.norm**2 * np.array([[s11, s12],
                                                   [s21, s22]])
        else:

            s11 = 1 + (s - 1) * v.x**2 / v.norm**2
            s22 = 1 + (s - 1) * v.y**2 / v.norm**2
            s33 = 1 + (s - 1) * v.z**2 / v.norm**2

            s12 = s21 = (s - 1) * v.x * v.y / v.norm**2
            s13 = s31 = (s - 1) * v.x * v.z / v.norm**2
            s23 = s32 = (s - 1) * v.y * v.z / v.norm**2

            Smat[:,:] = np.array([[s11, s12, s13],
                                  [s21, s22, s23],
                                  [s31, s32, s33]])

        Smat[np.where(np.abs(Smat) <= np.finfo(float).eps)] = 0.0
        return Smat


def rotate(obj, angle=None, rot_axis=None, anchor_point=None, rot_point=None,
           from_vector=None, to_vector=None, deg2rad=False,
           transform_matrix=None, verbose=False):
    """Rotate object.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    obj : array_like
    angle : float, optional
    rot_axis : {array_like, str, :class:`~sknano.core.math.Vector`}
    anchor_point : {array_like, :class:`~sknano.core.math.Point`}
    rot_point : :class:`~sknano.core.math.Point`, optional
    from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
    deg2rad : bool, optional
    transform_matrix : array_like, optional

    Returns
    -------
    rotated object : array_like

    """
    if verbose:
        print('In rotate\n'
              'obj: {}\n'.format(obj) +
              'rot_axis: {}\n'.format(rot_axis) +
              'anchor_point: {}\n'.format(anchor_point) +
              'transform_matrix: {}\n'.format(transform_matrix))

    t = np.append(np.asarray(obj), 1)
    try:
        rot_obj = np.dot(transform_matrix, t)[:-1]
    except TypeError:
        if rot_axis is None and len(obj) > 2:
            raise TypeError('`rot_axis` must be a sequence with the same '
                            'shape as `obj`')
        tmatrix = \
            transformation_matrix(angle=angle, rot_axis=rot_axis,
                                  anchor_point=anchor_point,
                                  rot_point=rot_point, from_vector=from_vector,
                                  to_vector=to_vector, deg2rad=deg2rad,
                                  verbose=verbose)

        rot_obj = np.dot(tmatrix, t)[:-1]

    if rot_obj.__class__ != obj.__class__:
        rot_obj = obj.__class__(rot_obj)

    return rot_obj


def scale():
    pass


def translate(obj, t):
    """Translate object points by a vector `t`."""
    pass
