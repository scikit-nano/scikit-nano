# -*- coding: utf-8 -*-
"""
=======================================================================
3D geometric regions (:mod:`sknano.core.geometric_regions._3D_regions`)
=======================================================================

.. currentmodule:: sknano.core.geometric_regions._3D_regions

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod

import numbers
import numpy as np

from sknano.core.math import Point, Vector, vector as vec
from ._base import GeometricRegion, GeometricTransformsMixin

__all__ = ['Geometric3DRegion', 'Parallelepiped', 'Cuboid', 'Cube',
           'Ellipsoid', 'Spheroid', 'Sphere', 'Cylinder']


class Geometric3DRegion(GeometricRegion, GeometricTransformsMixin,
                        metaclass=ABCMeta):
    """Abstract base class for representing 3D geometric regions."""

    @property
    @abstractmethod
    def volume(self):
        """Volume of 3D geometric region."""
        raise NotImplementedError

    @property
    def measure(self):
        """Measure of 3D geometric region."""
        return self.volume


class Parallelepiped(Geometric3DRegion):
    """Abstract representation of parallelepiped.

    .. versionadded:: 0.3.0

    Represents a parallelepiped with origin :math:`o` and directions
    :math:`u, v, w`.

    Parameters
    ----------
    o : array_like, optional
        parallelepiped origin
    u, v, w : array_like, optional
        parallelepiped direction vectors stemming from origin `o`.


    """
    def __init__(self, o=None, u=None, v=None, w=None):
        super().__init__()

        if o is None:
            o = [0, 0, 0]
        self.o = o

        if u is None:
            u = [1, 0, 0]
        self.u = u

        if v is None:
            v = [1, 1, 0]
        self.v = v

        if w is None:
            w = [0, 1, 1]
        self.w = w

        self.points.append(self.o)
        self.vectors.extend([self.u, self.v, self.w])

        self.fmtstr = "o={o!r}, u={u!r}, v={v!r}, w={w!r}"

    @property
    def o(self):
        return self._o

    @o.setter
    def o(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._o = Point(value)

    @property
    def u(self):
        return self._u

    @u.setter
    def u(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._u = Vector(value, p0=self.o)

    @property
    def v(self):
        return self._v

    @v.setter
    def v(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._v = Vector(value, p0=self.o)

    @property
    def w(self):
        return self._w

    @w.setter
    def w(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._w = Vector(value, p0=self.o)

    @property
    def volume(self):
        u = self.u
        v = self.v
        w = self.w
        return vec.scalar_triple_product(u, v, w)

    @property
    def centroid(self):
        ox, oy, oz = self.o
        ux, uy, uz = self.u
        vx, vy, vz = self.v
        wx, wy, wz = self.w

        xcom = 0.5 * (2 * ox + ux + vx + wx)
        ycom = 0.5 * (2 * oy + uy + vy + wy)
        zcom = 0.5 * (2 * oy + uy + vy + wz)

        return Point([xcom, ycom, zcom])

    def contains(self, point):
        """Check if point is contained within volume of `Parallelepiped`."""
        x, y, z = Point(point)

        ox, oy, oz = self.o
        ux, uy, uz = self.u
        vx, vy, vz = self.v
        wx, wy, wz = self.w

        q1 = (vz * (wx * (y - oy) + wy * (ox - x)) +
              wz * (vx * (oy - y) + vy * (x - ox)) +
              oz * (vy * wx - vx * wy) +
              z * (vx * wy - vy * wx)) / \
            (uz * (vx * wy - vy * wx) +
             uy * (vz * wx - vx * wz) +
             ux * (vy * wz - vz * wy))

        q2 = (uz * (wx * (y - oy) + wy * (ox - x)) +
              wz * (ux * (oy - y) + uy * (x - ox)) +
              oz * (uy * wx - ux * wy) +
              z * (ux * wy - uy * wx)) / \
            (uz * (vy * wx - vx * wy) +
             uy * (vx * wz - vz * wx) +
             ux * (vz * wy - vy * wz))

        q3 = (uz * (vx * (y - oy) + vy * (ox - x)) +
              vz * (ux * (oy - y) + uy * (x - ox)) +
              oz * (uy * vx - ux * vy) +
              z * (ux * vy - uy * vx)) / \
            (uz * (vx * wy - vy * wx) +
             uy * (vz * wx - vx * wz) +
             ux * (vy * wz - vz * wy))

        return q1 >= 0 and q1 <= 1 and q2 >= 0 and q2 <= 1 and \
            q3 >= 0 and q3 <= 1

    def todict(self):
        return dict(o=self.o, u=self.u, v=self.v, w=self.w)


class Cuboid(Geometric3DRegion):
    """Abstract representation of cuboid.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    pmin, pmax : sequence, optional
    xmin, ymin, zmin : float, optional
    xmax, ymax, zmax : float, optional

    """
    def __init__(self, pmin=None, pmax=None,
                 xmin=None, ymin=None, zmin=None,
                 xmax=None, ymax=None, zmax=None):
        super().__init__()

        if pmin is None:
            pmin = [xmin, ymin, zmin]
        self.pmin = pmin

        if pmax is None:
            pmax = [xmax, ymax, zmax]
        self.pmax = pmax

        self.points.extend([self.pmin, self.pmax])
        self.fmtstr = "pmin={pmin!r}, pmax={pmax!r}"

    @property
    def pmin(self):
        return self._pmin

    @pmin.setter
    def pmin(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._pmin = Point(value)

    @property
    def pmax(self):
        return self._pmax

    @pmax.setter
    def pmax(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._pmax = Point(value)

    @property
    def xmin(self):
        return self.pmin.x

    @xmin.setter
    def xmin(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.pmin.x = value

    @property
    def ymin(self):
        return self.pmin.y

    @ymin.setter
    def ymin(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.pmin.y = value

    @property
    def zmin(self):
        return self.pmin.z

    @zmin.setter
    def zmin(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.pmin.z = value

    @property
    def xmax(self):
        return self.pmax.x

    @xmax.setter
    def xmax(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.pmax.x = value

    @property
    def ymax(self):
        return self.pmax.y

    @ymax.setter
    def ymax(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.pmax.y = value

    @property
    def zmax(self):
        return self.pmax.z

    @zmax.setter
    def zmax(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self.pmax.z = value

    @property
    def centroid(self):
        xcom = 0.5 * (self.xmin + self.xmax)
        ycom = 0.5 * (self.ymin + self.ymax)
        zcom = 0.5 * (self.zmin + self.zmax)
        return Point([xcom, ycom, zcom])

    @property
    def a(self):
        return self.xmax - self.xmin

    @property
    def b(self):
        return self.ymax - self.ymin

    @property
    def c(self):
        return self.zmax - self.zmin

    @property
    def volume(self):
        a = self.a
        b = self.b
        c = self.c
        return a * b * c

    def contains(self, point):
        """Check if point is contained within volume of `Cuboid`."""
        x, y, z = Point(point)

        return (x >= self.xmin) and (x <= self.xmax) and \
            (y >= self.ymin) and (y <= self.ymax) and \
            (z >= self.zmin) and (z <= self.zmax)

    def todict(self):
        return dict(pmin=self.pmin, pmax=self.pmax)


class Cube(Geometric3DRegion):
    """Abstract data structure representing a cube.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    center : sequence, optional
        Either a 3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the :math:`(x,y,z)`
        coordinates of the cube center.
    a : float, optional
        length of edge

    """
    def __init__(self, center=None, a=1):
        super().__init__()

        if center is None:
            center = [0, 0, 0]
        self.center = center
        self.a = a
        self.points.append(self.center)

        self.fmtstr = "center={center!r}, a={a:.3f}"

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._center = Point(value)

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, value):
        self._a = value

    @property
    def centroid(self):
        self.center

    @property
    def volume(self):
        a = self.a
        return a**3

    def contains(self, point):
        x, y, z = Point(point)
        h, k, l = self.center
        a = self.a
        xmin = h - a / 2
        ymin = k - a / 2
        zmin = l - a / 2
        xmax = h + a / 2
        ymax = k + a / 2
        zmax = l + a / 2
        return (x >= xmin) and (x <= xmax) and \
            (y >= ymin) and (y <= ymax) and \
            (z >= zmin) and (z <= zmax)

    def todict(self):
        return dict(center=self.center, a=self.a)


class Ellipsoid(Geometric3DRegion):
    """Abstract data structure representing an ellipsoid.

    .. versionadded:: 0.3.0

    The general ellipsoid is a quadratic surface with is given in
    Cartesian coordinates by:

    .. math::

       \\frac{x^2}{a^2} + \\frac{y^2}{b^2} + \\frac{z^2}{c^2} = 1

    Parameters
    ----------
    center : sequence or :class:`~sknano.core.Point`
        Either a 3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the :math:`(x,y,z)`
        coordinates of the :class:`Ellipsoid` center.
    a, b, c : float, optional
        Semi-principal axes :math:`a, b, c` of axis-aligned
        :class:`Ellipsoid`

    """
    def __init__(self, center=None, a=None, b=None, c=None):
        super().__init__()

        if center is None:
            center = [0, 0, 0]
        self.center = center

        if a is None:
            a = 0.5
        if b is None:
            b = 0.5
        if c is None:
            c = 0.5
        self.a, self.b, self.c = a, b, c

        self.points.append(self.center)

        self.fmtstr = "center={center!r}, a={a:.3f}, b={b:.3f}, c={c:.3f}"

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._center = Point(value)

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, value):
        self._a = value

    @property
    def b(self):
        return self._b

    @b.setter
    def b(self, value):
        self._b = value

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, value):
        self._c = value

    @property
    def centroid(self):
        return self.center

    @property
    def volume(self):
        a = self.a
        b = self.b
        c = self.c
        return 4 / 3 * np.pi * a * b * c

    def contains(self, point):
        """Check if point is contained within volume of :class:`Ellipsoid`."""
        x, y, z = Point(point)
        h, k, l = self.center
        a, b, c = self.a, self.b, self.c

        return (x - h)**2 / a**2 + (y - k)**2 / b**2 + \
            (z - l)**2 / c**2 <= 1.0

    def todict(self):
        return dict(center=self.center, a=self.a, b=self.b, c=self.c)


class Spheroid(Geometric3DRegion):
    """Abstract data structure representing a spheroid.

    .. versionadded:: 0.3.0

    The general spheroid is a quadratic surface with is given in
    Cartesian coordinates by:

    .. math::

       \\frac{x^2 + y^2}{a^2} + \\frac{z^2}{c^2} = 1

    Parameters
    ----------
    center : sequence, optional
        Either a 3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the :math:`(x,y,z)`
        coordinates of the :class:`Spheroid` center.
    a, c : float, optional
        Semi-axes :math:`a, c` of axis-aligned
        :class:`Spheroid` with symmetry axis along :math:`z` axis.
    """
    def __init__(self, center=None, a=None, c=None):
        super().__init__()

        if center is None:
            center = [0, 0, 0]
        self.center = center

        if a is None:
            a = 0.5
        if c is None:
            c = 0.5

        self.a, self.c = a, c

        self.points.append(self.center)

        self.fmtstr = "center={center!r}, a={a:.3f}, c={c:.3f}"

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._center = Point(value)

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, value):
        self._a = value

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, value):
        self._c = value

    @property
    def centroid(self):
        return self.center

    @property
    def volume(self):
        a = self.a
        c = self.c
        return 4 / 3 * np.pi * a**2 * c

    def contains(self, point=None):
        """Check if point is contained within volume of :class:`Spheroid`."""
        x, y, z = Point(point)
        h, k, l = self.center
        a, c = self.a, self.c

        return ((x - h)**2 + (y - k)**2) / a**2 + (z - l)**2 / c**2 <= 1.0

    def todict(self):
        return dict(center=self.center, a=self.a, c=self.c)


class Sphere(Geometric3DRegion):
    """Abstract data structure representing a sphere.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    center : sequence, optional
        Either a 3-tuple of floats or an instance of the
        :class:`~sknano.core.Point` class specifying the :math:`(x,y,z)`
        coordinates of the :class:`Sphere` center.
    r : float, optional
        Sphere radius :math:`r`
    """
    def __init__(self, center=None, r=1):
        super().__init__()

        if center is None:
            center = [0, 0, 0]
        self.center = center

        self.r = r

        self.points.append(self.center)

        self.fmtstr = "center={center!r}, r={r:.3f}"

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, value):
        self._r = value

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or \
                len(value) != 3:
            raise TypeError('Expected a 3-element array_like object')
        self._center = Point(value)

    @property
    def centroid(self):
        return self.center

    @property
    def volume(self):
        r = self.r
        return 4 / 3 * np.pi * r**3

    def contains(self, point):
        x, y, z = point
        h, k, l = self.center
        r = self.r

        return (x - h)**2 + (y - k)**2 + (z - l)**2 <= r**2

    def todict(self):
        return dict(center=self.center, r=self.r)


class Cylinder(Geometric3DRegion):
    """`Geometric3DRegion` for a cylinder.

    .. versionadded:: 0.3.10

    Parameters
    ----------
    p1, p2 : sequence, optional
        3-tuples or :class:`~sknano.core.Point` class instances
        specifying the `Cylinder` axis from point :math:`p_1=(x_1,y_1,z_1)`
        to :math:`p_2=(x_2,y_2,z_2)`.
    r : float, optional
        `Cylinder` radius :math:`r`
    """
    def __init__(self, p1=None, p2=None, r=1):
        super().__init__()

        if p1 is None:
            p1 = [0, 0, -1]
        if p2 is None:
            p2 = [0, 0, 1]

        self.p1 = p1
        self.p2 = p2
        self.points.extend([self.p1, self.p2])
        self.r = r
        self.fmtstr = "p1={p1!r}, p2={p2!r}, r={r:.3f}"

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, value):
        self._r = value

    @property
    def p1(self):
        return self._p1

    @p1.setter
    def p1(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)):
            raise TypeError('Expected a 3-element array_like object')
        self._p1 = Point(value)

    @property
    def p2(self):
        return self._p2

    @p2.setter
    def p2(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)):
            raise TypeError('Expected a 3-element array_like object')
        self._p2 = Point(value)

    @property
    def axis(self):
        return Vector(p0=self.p1, p=self.p2)

    @property
    def centroid(self):
        x1, y1, z1 = self.p1
        x2, y2, z2 = self.p2

        xcom = 0.5 * (x1 + x2)
        ycom = 0.5 * (y1 + y2)
        zcom = 0.5 * (z1 + z2)
        return Point([xcom, ycom, zcom])

    @property
    def volume(self):
        return np.pi * self.r**2 * self.axis.length

    def contains(self, point):
        px, py, pz = Point(point)
        x1, y1, z1 = self.p1
        x2, y2, z2 = self.p2
        r = self.r

        q1 = ((px - x1) * (x2 - x1) +
              (py - y1) * (y2 - y1) +
              (pz - z1) * (z2 - z1)) / \
            ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

        return r > 0 and (not np.allclose(x1, x2) or
                          not np.allclose(y1, y2) or
                          not np.allclose(z1, z2)) and \
            q1 >= 0 and q1 <= 1 and (x1 - px + (x2 - x1) * q1)**2 + \
            (y1 - py + (y2 - y1) * q1)**2 + (z1 - pz + (z2 - z1) * q1)**2 <= r

    def todict(self):
        return dict(p1=self.p1, p2=self.p2, r=self.r)
