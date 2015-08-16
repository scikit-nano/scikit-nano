# -*- coding: utf-8 -*-
"""
===============================================================================
Functions for linear algebra transforms (:mod:`sknano.core.math._transforms`)
===============================================================================

.. currentmodule:: sknano.core.math._transforms

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

__docformat__ = 'restructuredtext en'

from itertools import combinations

import numpy as np

__all__ = ['rotate', 'Rx', 'Ry', 'Rz', 'rotation_matrix',
           'reflection_matrix', 'scaling_matrix', 'translation_matrix',
           'translation_from_matrix',
           'transformation_matrix', 'axis_angle_from_rotation_matrix']

I = np.identity(4)

_str2array = {axis: I[:3, i] for i, axis in enumerate(('x', 'y', 'z'))}


def rotate(obj, angle=None, axis=None, anchor_point=None, rot_point=None,
           from_vector=None, to_vector=None, degrees=False,
           transform_matrix=None, verbose=False, **kwargs):
    """Rotate object.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    obj : array_like
    angle : float, optional
        Rotation angle in **radians** *unless* `degrees` is `True`.
    axis : {array_like, str, :class:`~sknano.core.math.Vector`}
        Rotation axis.
    anchor_point : {array_like, :class:`~sknano.core.math.Point`}
    rot_point : :class:`~sknano.core.math.Point`, optional
    from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
    degrees : bool, optional
    transform_matrix : array_like, optional

    Returns
    -------
    rotated object : array_like

    """
    if axis is None and 'rot_axis' in kwargs:
        axis = kwargs['rot_axis']
        del kwargs['rot_axis']

    if verbose:
        print('In rotate\n'
              'obj: {}\n'.format(obj) +
              'axis: {}\n'.format(axis) +
              'anchor_point: {}\n'.format(anchor_point) +
              'transform_matrix: {}\n'.format(transform_matrix))

    objarr = np.asarray(obj)
    nd = len(objarr)
    objarr = np.append(objarr, 1)

    if transform_matrix is not None and \
            isinstance(transform_matrix, np.ndarray):

        if transform_matrix.shape[-1] != nd + 1:
            augmented_tmatrix = np.eye(nd + 1)
            augmented_tmatrix[:nd, :nd] = transform_matrix
            transform_matrix = augmented_tmatrix

    try:
        rot_obj = np.dot(transform_matrix, objarr)[:-1]
    except TypeError:
        if axis is None and from_vector is None and to_vector is None and \
                nd > 2:
            raise ValueError('Expected `axis` to be a 3D `Vector`.')
        transform_matrix = \
            transformation_matrix(angle=angle, axis=axis,
                                  anchor_point=anchor_point,
                                  rot_point=rot_point, from_vector=from_vector,
                                  to_vector=to_vector, degrees=degrees,
                                  verbose=verbose, **kwargs)
        return rotate(obj, transform_matrix=transform_matrix)
    else:
        if rot_obj.__class__ != obj.__class__:
            rot_obj = obj.__class__(rot_obj)
        return rot_obj


def Rx(angle, degrees=False):
    """Generate the :math:`3\\times3` rotation matrix :math:`R_x(\\theta)` \
        for a rotation about the :math:`x` axis by an angle :math:`\\theta`.

    Parameters
    ----------
    angle : float
        The rotation angle :math:`\\theta` in *radians*. If the angle is
        given in *degrees*, then you must set `degrees=True` to correctly
        calculate the rotation matrix.
    degrees : bool, optional
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
    >>> np.alltrue(Rx(np.pi/4) == Rx(45, degrees=True))
    True

    """
    if degrees:
        angle = np.radians(angle)
    cosa = np.cos(angle)
    sina = np.sin(angle)

    Rmat = np.array([[1.0, 0.0, 0.0], [0.0, cosa, -sina], [0.0, sina, cosa]])
    Rmat[np.where(np.abs(Rmat) <= np.finfo(float).eps)] = 0.0
    return Rmat


def Ry(angle, degrees=False):
    """Generate the :math:`3\\times3` rotation matrix :math:`R_y(\\theta)` \
        for a rotation about the :math:`y` axis by an angle :math:`\\theta`.

    Parameters
    ----------
    angle : float
        The rotation angle :math:`\\theta` in *radians*. If the angle is
        given in *degrees*, then you must set `degrees=True` to correctly
        calculate the rotation matrix.
    degrees : bool, optional
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
    >>> np.alltrue(Ry(np.pi/4) == Ry(45, degrees=True))
    True

    """
    if degrees:
        angle = np.radians(angle)
    cosa = np.cos(angle)
    sina = np.sin(angle)

    Rmat = np.array([[cosa, 0.0, sina], [0.0, 1.0, 0.0], [-sina, 0.0, cosa]])
    Rmat[np.where(np.abs(Rmat) <= np.finfo(float).eps)] = 0.0
    return Rmat


def Rz(angle, degrees=False):
    """Generate the :math:`3\\times3` rotation matrix :math:`R_z(\\theta)` \
        for a rotation about the :math:`z` axis by an angle :math:`\\theta`.

    Parameters
    ----------
    angle : float
        The rotation angle :math:`\\theta` in *radians*. If the angle is
        given in *degrees*, then you must set `degrees=True` to correctly
        calculate the rotation matrix.
    degrees : bool, optional
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
    >>> np.alltrue(Rz(np.pi/4) == Rz(45, degrees=True))
    True

    """
    if degrees:
        angle = np.radians(angle)
    cosa = np.cos(angle)
    sina = np.sin(angle)

    Rmat = np.array([[cosa, -sina, 0.0], [sina, cosa, 0.0], [0.0, 0.0, 1.0]])
    Rmat[np.where(np.abs(Rmat) <= np.finfo(float).eps)] = 0.0
    return Rmat


def reflection_matrix(v):
    """Generate reflection matrix that represents reflection of points \
        in a mirror plane defined by normal `Vector` `v`.

    Parameters
    ----------
    v : :class:`~sknano.core.math.Vector`

    Returns
    -------
    :class:`~numpy:numpy.ndarray`

    """
    from ._vector import Vector
    v = Vector(v)

    Rmat = np.zeros((v.nd, v.nd))
    for i in range(v.nd):
        Rmat[i, i] = 1 - 2 * v[i] ** 2 / v.norm ** 2

    for i, j in combinations(list(range(v.nd)), 2):
        Rmat[i, j] = Rmat[j, i] = -2 * v[i] * v[j] / v.norm ** 2

    Rmat[np.where(np.abs(Rmat) <= np.finfo(float).eps)] = 0.0
    return Rmat


def rotation_matrix(angle=None, axis=None, anchor_point=None,
                    rot_point=None, from_vector=None, to_vector=None,
                    degrees=False, verbose=False, **kwargs):
    """Generate an :math:`n\\times n` rotation matrix.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    angle : float
        Rotation angle in **radians**. If `degrees` is `True`, `angle` will be
        converted to radians from degrees.  The *sense* of the rotation is
        defined by the *right hand rule*: If your right-hand's thumb points
        along the `axis`, then your fingers wrap around the axis in the
        *positive sense* of the rotation angle.
    axis : {None, array_like, str}, optional
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
    degrees : bool, optional
        If `True`, convert `angle` from degrees to radians.

    Returns
    -------
    Rmat : :class:`~numpy:numpy.ndarray`
        If `axis` is `None` then `Rmat` will be a :math:`2D`
        rotation matrix :math:`R(\\theta)` that rotates :math:`2D` vectors
        counterclockwise by `angle` :math:`\\theta`.

        If `axis` is not `None` then `Rmat` will be a rotation matrix that
        gives a rotation around the direction of the vector `axis`.

    """
    if axis is None and 'rot_axis' in kwargs:
        axis = kwargs['rot_axis']
        del kwargs['rot_axis']

    Rmat = transformation_matrix(angle=angle, axis=axis,
                                 anchor_point=anchor_point,
                                 rot_point=rot_point, from_vector=from_vector,
                                 to_vector=to_vector, degrees=degrees,
                                 verbose=verbose, **kwargs)
    return Rmat[:-1, :-1]


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

        Smat = np.zeros((v.nd + 1, v.nd + 1))
        for i in range(v.nd):
            Smat[i, i] = 1 + (s - 1) * v[i] ** 2 / v.norm ** 2

        for i, j in combinations(list(range(v.nd)), 2):
            Smat[i, j] = Smat[j, i] = \
                (s - 1) * v[i] ** 2 * v[j] ** 2 / v.norm ** 2

        Smat[np.where(np.abs(Smat) <= np.finfo(float).eps)] = 0.0
        return Smat


def translation_matrix(v):
    """Generate :math:`4\\times4` translation matrix given 3D translation \
        `Vector` `v`.

    Parameters
    ----------
    v : `Vector`
        Translation vector.

    Returns
    -------
    M : :class:`~numpy:numpy.matrix`
        Translation matrix
    """
    from . import Vector
    v = Vector(v)
    nd = v.nd
    M = np.identity(nd + 1)
    M[:nd, nd] = v[:nd]
    return np.asmatrix(M)


def translation_from_matrix(M):
    from . import Vector
    return Vector(np.asarray(M)[:-1, -1])


def transformation_matrix(angle=None, axis=None, anchor_point=None,
                          rot_point=None, from_vector=None, to_vector=None,
                          degrees=False, verbose=False, **kwargs):
    """Generate an :math:`(n+1)\\times(n+1)` transformation matrix for an \
        affine transformation in :math:`n` dimensions.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    angle : float
        Rotation angle in **radians**. If `degrees` is `True`, `angle` will be
        converted to radians from degrees.  The *sense* of the rotation is
        defined by the *right hand rule*: If your right-hand's thumb points
        along the `axis`, then your fingers wrap around the axis in the
        *positive sense* of the rotation angle.
    axis : {None, array_like, str}, optional
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

        If `anchor_point` is not `None` and `axis` is a `Vector` instance,
        then the origin of the vector defined by :attr:`Vector.p0` will be
        changed to `anchor_point`.

        If `anchor_point` is `None`, then it defaults to an
        :math:`n`-element array of zeros.
    degrees : bool, optional
        If `True`, convert `angle` from degrees to radians.

    Returns
    -------
    Tmat : ndarray
        :math:`n+1\\times n+1` transformation matrix for an affine transform
        in :math:`n` dimensions.

        If `axis` is `None` and `anchor_point` is `None`,
        then `Tmat` will be a :math:`2D` rotation matrix :math:`R(\\theta)`
        that rotates :math:`2D` vectors counterclockwise by `angle`
        :math:`\\theta`.

        If `axis` is `None` and `anchor_point` is a 2-element sequence,
        then `Rmat` will be a :math:`2D` rotation matrix :math:`R(\\theta)`
        about the :math:`2D` `Point` `anchor_point` by `angle`
        :math:`\\theta`.

        If `axis` is not `None` and `anchor_point` is `None`,
        then `Rmat` will be a rotation matrix that gives a rotation around
        the direction of the vector `axis`.

    Notes
    -----

    """
    if angle is None and from_vector is None and to_vector is None:
        raise TypeError('Expected `angle` or `from_vector` and `to_vector` '
                        'set to values with valid types.')

    if angle is not None and degrees:
        angle = np.radians(angle)

    if axis is None and 'rot_axis' in kwargs:
        axis = kwargs['rot_axis']
        del kwargs['rot_axis']

    from . import vector, Vector, Point, NullVector

    if angle is None and \
            all([v is not None for v in (from_vector, to_vector)]):
        from_vector = Vector(from_vector)
        to_vector = Vector(to_vector)
        if from_vector.nd != to_vector.nd:
            raise ValueError('`from_vector` and `to_vector` must have same '
                             'dimensions.')
        angle = vector.angle(from_vector, to_vector)
        if from_vector.nd == 2:
            if not np.allclose(Vector(np.dot(Rz(angle)[:2, :2], from_vector)),
                               to_vector):
                angle = -angle

            return transformation_matrix(angle=angle, rot_point=rot_point)

        elif from_vector.nd == 3:
            axis = vector.cross(from_vector, to_vector)
            if axis == NullVector:
                return I

            return transformation_matrix(angle=angle, axis=axis,
                                         anchor_point=anchor_point,
                                         rot_point=rot_point)
        else:
            raise ValueError('currently, transformation_matrices can only '
                             'be computed for rotation transformations in '
                             '2D and 3D.')
    else:
        cosa = np.cos(angle)
        sina = np.sin(angle)

        if axis is None:

            if rot_point is not None and \
                    isinstance(rot_point, (list, np.ndarray)) and \
                    len(rot_point) == 2:
                p = Point(rot_point)
            else:
                p = Point(nd=2)

            Tmat = Rz(angle)

            if not np.allclose(p, np.zeros(2)):
                for i in range(2):
                    j, = list(set(range(2)) - {i})
                    Tmat[i, 2] = p[i] - p[i] * cosa + (-1) ** i * p[j] * sina

        else:
            # Handle 3D rotation about origin
            # Handle 3D rotation around the rotation vector
            # Handle N-D rotation about origin
            # Handle N-D rotation around the N-D rotation vector anchored at
            # an arbitrary N-D point.

            if isinstance(axis, str):
                try:
                    axis = _str2array[axis]
                except KeyError:
                    raise ValueError(
                        'Invalid `axis` string: {}'.format(axis))
            elif not isinstance(axis, (tuple, list, np.ndarray)):
                raise ValueError('`axis` must be a sequence')

            if anchor_point is None and rot_point is not None:
                anchor_point = rot_point

            if anchor_point is None and isinstance(axis, Vector):
                anchor_point = axis.p0

            if anchor_point is not None and \
                    isinstance(anchor_point, (list, np.ndarray)) and \
                    len(anchor_point) == 3:
                anchor_point = Point(anchor_point)
            else:
                anchor_point = Point()

            v = Vector(np.asarray(axis), p0=anchor_point)

            Tmat = np.zeros((4, 4))

            if np.allclose(v, I[:3, 0]):
                Tmat[:3, :3] = Rx(angle)
            elif np.allclose(v, I[:3, 1]):
                Tmat[:3, :3] = Ry(angle)
            elif np.allclose(v, I[:3, 2]):
                Tmat[:3, :3] = Rz(angle)
            else:
                for i in range(3):
                    Tmat[i, i] = (v[i] ** 2 +
                                  (v[(i + 1) % 3] ** 2 +
                                   v[(i + 2) % 3] ** 2) * cosa) / v.norm ** 2

                for i, j in combinations(list(range(3)), 2):
                    k = list(set(range(3)) - {i, j})[0]
                    Tmat[i, j] = \
                        (v[i] * v[j] * (1 - cosa) +
                         (-1) ** (i + j) * v[k] * v.norm * sina) / v.norm ** 2
                    Tmat[j, i] = \
                        (v[i] * v[j] * (1 - cosa) -
                         (-1) ** (i + j) * v[k] * v.norm * sina) / v.norm ** 2

            if not np.allclose(v.p0, np.zeros(3)):

                p = v.p0

                for i in range(3):
                    j, k = list(set(range(3)) - {i})
                    Tmat[i, 3] = ((p[i] * (v[j] ** 2 + v[k] ** 2) -
                                   v[i] * (p[j] * v[j] + p[k] * v[k])) *
                                  (1 - cosa) + (-1) ** i *
                                  (p[j] * v[k] - p[k] * v[j]) *
                                  v.norm * sina) / v.norm ** 2

            Tmat[3, 3] = 1.0

        Tmat[np.where(np.abs(Tmat) <= np.finfo(float).eps)] = 0.0

        return Tmat


def axis_angle_from_rotation_matrix(rmatrix):
    """Compute the rotation axis and angle from a rotation matrix.

    Parameters
    ----------
    rmatrix : :class:`~numpy:numpy.ndarray`

    Returns
    -------
    axis : :class:`~sknano.core.math.Vector`
    angle : :class:`~numpy:numpy.float`

    .. todo::

       Fix code to compute the correct sense of the rotation or
       direction of rotation axis vector.

    """
    if not isinstance(rmatrix, np.ndarray):
        raise TypeError('Expected matrix to be a 3x3 numpy array.')
    from ._vector import Vector

    angle = np.arccos((np.trace(rmatrix) - 1) / 2)

    # compute random vector
    x = Vector(np.random.random(3)).unit_vector
    axis = (np.dot(rmatrix, x) +
            np.dot(rmatrix.T, x) +
            (1 - np.trace(rmatrix)) * x).unit_vector

    axis.rezero()

    return axis, angle
