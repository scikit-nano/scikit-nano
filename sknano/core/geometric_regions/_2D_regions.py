# -*- coding: utf-8 -*-
"""
=======================================================================
2D geometric regions (:mod:`sknano.core.geometric_regions._2D_regions`)
=======================================================================

.. currentmodule:: sknano.core.geometric_regions._2D_regions

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod

import numpy as np

from sknano.core.math import Point, Vector
from ._base import GeometricRegion, GeometricTransformsMixin

__all__ = ['Geometric2DRegion', 'Parallelogram', 'Rectangle', 'Square',
           'Ellipse', 'Circle', 'Triangle']


class Geometric2DRegion(GeometricRegion, GeometricTransformsMixin,
                        metaclass=ABCMeta):
    """Abstract base class for representing 2D geometric regions."""

    @property
    @abstractmethod
    def area(self):
        """Area of 2D geometric region."""
        raise NotImplementedError

    @property
    def measure(self):
        """Alias for :attr:`~Geometric2DRegion.area`, which is the measure \
            of a 2D geometric region."""
        return self.area


class Parallelogram(Geometric2DRegion):
    """`Geometric2DRegion` for a parallelogram.

    .. versionadded:: 0.3.0

    Represents a parallelogram with origin :math:`o=(o_x, o_y)` and
    direction vectors :math:`\\mathbf{u}=(u_x, u_y)` and
    :math:`\\mathbf{v}=(v_x, v_y)`.

    Parameters
    ----------
    o : array_like, optional
        Parallelogram origin. If `None`, it defaults to `o=[0, 0]`.
    u, v : array_like, optional
        Parallelogram direction vectors stemming from origin `o`.
        If `None`, then the default values are `u=[1, 0]` and `v=[1, 1]`.

    Notes
    -----
    :class:`Parallelogram` represents the bounded region
    :math:`\\left \\{o+\\lambda_1\\mathbf{u}+\\lambda_2\\mathbf{v}\\in R^2
    |0\\le\\lambda_i\\le 1\\right\\}`, where :math:`\\mathbf{u}` and
    :math:`\\mathbf{v}` have to be linearly independent.

    Calling :class:`Paralleogram` with no parameters is equivalent to
    :class:`Parallelogram`\ `(o=[0, 0], u=[1, 0], v=[1, 1])`

    """
    def __init__(self, o=None, u=None, v=None):

        super().__init__()

        self._o = Point(nd=2)
        self._u = Vector(nd=2)
        self._v = Vector(nd=2)

        if o is None:
            o = [0, 0]
        self.o = o

        if u is None:
            u = [1, 0]
        self.u = u

        if v is None:
            v = [1, 1]
        self.v = v

        self.points.append(self.o)
        self.vectors.extend([self.u, self.v])
        self.fmtstr = "o={o!r}, u={u!r}, v={v!r}"

    @property
    def o(self):
        """2D point coordinates :math:`(o_x, o_y)` of origin.

        Returns
        -------
        :class:`Point`
            2D :class:`Point` coordinates :math:`(o_x, o_y)` of origin.
        """
        return self._o

    @o.setter
    def o(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self.vectors.translate(Vector(p0=self.o, p=Point(value)))
        self._o[:] = Point(value)

    @property
    def u(self):
        """2D direction vector :math:`\\mathbf{u}=(u_x, u_y)`, with origin \
            :attr:`~Paralleogram.o`"""
        return self._u

    @u.setter
    def u(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._u[:] = Vector(value, p0=self.o)

    @property
    def v(self):
        """2D direction vector :math:`\\mathbf{v}=(v_x, v_y)`, with origin \
            :attr:`~Paralleogram.o`"""
        return self._v

    @v.setter
    def v(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._v[:] = Vector(value, p0=self.o)

    @property
    def area(self):
        """:class:`Paralleogram` area, \
            :math:`A=|\\mathbf{u}\\times\\mathbf{v}|`.

        Computed as:

        .. math::

           A = |\\mathbf{u}\\times\\mathbf{v}|

        """
        u = self.u
        v = self.v
        return np.abs(np.cross(u, v))

    @property
    def centroid(self):
        """Paralleogram centroid, :math:`(c_x, c_y)`.

        Computed as the 2D point :math:`(c_x, c_y)` with coordinates:

        .. math::

           c_x = o_x + \\frac{u_x + v_x}{2}

           c_y = o_y + \\frac{u_y + v_y}{2}

        where :math:`(o_x, o_y)`, :math:`(u_x, u_y)`, and :math:`(v_x, v_y)`
        are the :math:`(x, y)` coordinates of the origin :math:`o`
        and :math:`(x, y)` components of the direction vectors
        :math:`\\mathbf{u}` and :math:`\\mathbf{v}`, respectively.

        Returns
        -------
        :class:`Point`
            2D :class:`Point` of centroid.

        """
        ox, oy = self.o
        ux, uy = self.u
        vx, vy = self.v

        cx = ox + (ux + vx) / 2
        cy = oy + (uy + vy) / 2

        return Point([cx, cy])

    def contains(self, point):
        """Test region membership of `point` in :class:`Parallelogram`.

        Parameters
        ----------
        point : array_like

        Returns
        -------
        :class:`~python:bool`
            `True` if `point` is within :class:`Paralleogram`,
            `False` otherwise

        Notes
        -----
        A `point` :math:`(p_x, p_y)` is within the bounded region of
        a parallelogram with origin :math:`(o_x, o_y)` and direction
        vectors :math:`\\mathbf{u}=(u_x, u_y)` and
        :math:`\\mathbf{v}=(v_x, v_y)` if the following is true:

        .. math::

           0\\le\\frac{(p_y - o_y) v_x + (o_x - p_x) v_y}{u_y v_x - u_x v_y}
           \\le 1 \\land
           0\\le\\frac{(p_y - o_y) u_x + (o_x - p_x) u_y}{u_x v_y - u_y v_x}
           \\le 1

        """
        px, py = Point(point)

        ox, oy = self.o
        ux, uy = self.u
        vx, vy = self.v

        q1 = ((py - oy) * vx + (ox - px) * vy) / (uy * vx - ux * vy)
        q2 = ((py - oy) * ux + (ox - px) * uy) / (ux * vy - uy * vx)

        return q1 >= 0 and q1 <= 1 and q2 >= 0 and q2 <= 1

    def todict(self):
        """Returns a :class:`~python:dict` of the :class:`Paralleogram` \
            constructor parameters."""
        return dict(o=self.o, u=self.u, v=self.v)


class Rectangle(Geometric2DRegion):
    """`Geometric2DRegion` for a rectangle.

    .. versionadded:: 0.3.0

    Represents an axis-aligned bounded region from
    :math:`p_{\\mathrm{min}}=(x_{\\mathrm{min}},y_{\\mathrm{min}})` to
    :math:`p_{\\mathrm{max}}=(x_{\\mathrm{max}},y_{\\mathrm{max}})`.

    Parameters
    ----------
    pmin, pmax : array_like, optional
        The minimum and maximum 2D point coordinates of the axis-aligned
        rectangle from `pmin=[xmin, ymin]` to `pmax=[xmax, ymax]`.
    xmin, ymin : float, optional
        The minimum :math:`(x, y)` point of the axis-aligned rectangle.
    xmax, ymax : float, optional
        The maximum :math:`(x, y)` point of the axis-aligned rectangle.

    Notes
    -----
    :class:`Rectangle` represents the region
    :math:`\\left\\{\\{x, y\\}|x_{\\mathrm{min}}\\le x\\le x_{\\mathrm{max}}
    \\land y_{\\mathrm{min}}\\le y\\le y_{\\mathrm{max}}\\right\\}`

    Calling :class:`Rectangle` with no parameters is equivalent to
    :class:`Rectangle`\ `(pmin=[0, 0], pmax=[1, 1])`.

    """
    def __init__(self, pmin=None, pmax=None, xmin=0, ymin=0, xmax=1, ymax=1):

        super().__init__()

        self._pmin = Point(nd=2)
        self._pmax = Point(nd=2)

        if pmin is None:
            pmin = [xmin, ymin]
        self.pmin = pmin

        if pmax is None:
            pmax = [xmax, ymax]
        self.pmax = pmax

        self.points.append([self.pmin, self.pmax])

        self.fmtstr = "pmin={pmin!r}, pmax={pmax!r}"

    @property
    def pmin(self):
        """2D :class:`Point` at \
            (:attr:`~Rectangle.xmin`, :attr:`~Rectangle.ymin`)."""
        return self._pmin

    @pmin.setter
    def pmin(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._pmin[:] = Point(value)

    @property
    def pmax(self):
        """2D :class:`Point` at \
            (:attr:`~Rectangle.xmax`, :attr:`~Rectangle.ymax`)."""
        return self._pmax

    @pmax.setter
    def pmax(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._pmax[:] = Point(value)

    @property
    def xmin(self):
        """:math:`x_{\\mathrm{min}}` coordinate."""
        return self.pmin.x

    @xmin.setter
    def xmin(self, value):
        self.pmin.x = float(value)

    @property
    def xmax(self):
        """:math:`x_{\\mathrm{max}}` coordinate."""
        return self.pmax.x

    @xmax.setter
    def xmax(self, value):
        self.pmax.x = float(value)

    @property
    def ymin(self):
        """:math:`y_{\\mathrm{min}}` coordinate."""
        return self.pmin.y

    @ymin.setter
    def ymin(self, value):
        self.pmin.y = float(value)

    @property
    def ymax(self):
        """:math:`y_{\\mathrm{max}}` coordinate."""
        return self.pmax.y

    @ymax.setter
    def ymax(self, value):
        self.pmax.y = float(value)

    @property
    def a(self):
        """Distance between :math:`x_{\\mathrm{max}}-x_{\\mathrm{min}}`."""
        return self.xmax - self.xmin

    @property
    def b(self):
        """Distance between :math:`y_{\\mathrm{max}}-y_{\\mathrm{min}}`."""
        return self.ymax - self.ymin

    @property
    def area(self):
        """:class:`Rectangle` area, :math:`A=ab`"""
        a = self.a
        b = self.b
        return a * b

    @property
    def centroid(self):
        """:class:`Rectangle` centroid, :math:`(c_x, c_y)`.

        Computed as the 2D :class:`Point` :math:`(c_x, c_y)` with
        coordinates:

        .. math::

           c_x = \\frac{x_{\\mathrm{min}}+x_{\\mathrm{max}}}{2}

           c_y = \\frac{y_{\\mathrm{min}}+y_{\\mathrm{max}}}{2}

        Returns
        -------
        :class:`Point`
            2D :class:`Point` of centroid.

        """
        cx = (self.xmax + self.xmin) / 2
        cy = (self.ymax + self.ymin) / 2
        return Point([cx, cy])

    def contains(self, point):
        """Test region membership of `point` in :class:`Rectangle`.

        Parameters
        ----------
        point : array_like

        Returns
        -------
        :class:`~python:bool`
            `True` if `point` is within :class:`Rectangle`,
            `False`, otherwise.

        Notes
        -----
        A point :math:`(p_x, p_y)` is within the bounded region of a
        rectangle with lower corner at
        :math:`p_{\\mathrm{min}}=
        (x_{\\mathrm{min}}, y_{\\mathrm{min}})` and
        upper corner at
        :math:`p_{\\mathrm{max}}=
        (x_{\\mathrm{max}}, y_{\\mathrm{max}})` if the
        following is true:

        .. math::

           x_{\mathrm{min}}\\le x\\le x_{\\mathrm{max}}\\land

           y_{\mathrm{min}}\\le y\\le y_{\\mathrm{max}}

        """
        px, py = Point(point)
        xmin = self.xmin
        xmax = self.xmax
        ymin = self.ymin
        ymax = self.ymax

        return (px >= xmin) and (px <= xmax) and (py >= ymin) and (py <= ymax)

    def todict(self):
        """Returns a :class:`~python:dict` of the :class:`Rectangle` \
            constructor parameters."""
        return dict(pmin=self.pmin, pmax=self.pmax)


class Square(Geometric2DRegion):
    """`Geometric2DRegion` for a square.

    .. versionadded:: 0.3.0

    Represents the axis-aligned bounded region with center
    :math:`(c_x, c_y)` and side length :math:`a`.

    Parameters
    ----------
    center : array_like, optional
        The center point coordinate :math:`(c_x, c_y)`
        of the axis-aligned square.
    a : float, optional
        The side length :math:`a` of the axis-aligned square.

    Notes
    -----
    :class:`Square` represents the region
    :math:`\\left\\{\\{c_i\\pm\\frac{a}{2}\\}|a>0\\forall
    i\\in\\{x,y\\}\\right\\}`.

    Calling :class:`Square` with no parameters is equivalent to
    :class:`Square`\ `(center=[0, 0], a=1)`.

    """
    def __init__(self, center=None, a=1):

        super().__init__()

        self._center = Point(nd=2)

        if center is None:
            center = [0, 0]
        self.center = center

        self.a = a

        self.points.append(self.center)
        self.fmtstr = "center={center!r}, a={a:.2f}"

    @property
    def center(self):
        """Center point :math:`(c_x, c_y)` of axis-aligned square."""
        return self._center

    @center.setter
    def center(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._center[:] = Point(value)

    @property
    def a(self):
        """Side length :math:`a` of the axis-aligned square."""
        return self._a

    @a.setter
    def a(self, value):
        self._a = float(value)

    @property
    def area(self):
        """:class:`Square` area, :math:`A=a^2`."""
        return self.a ** 2

    @property
    def centroid(self):
        """Alias for :attr:`~Square.center`."""
        return self.center

    def contains(self, point):
        """Test region membership of `point` in :class:`Square`.

        Parameters
        ----------
        point : array_like

        Returns
        -------
        :class:`~python:bool`
            `True` if `point` is within :class:`Square`,
            `False` otherwise.

        Notes
        -----
        A `point` :math:`(p_x, p_y)` is within the bounded region of a
        square with center :math:`(c_x, c_y)` and side length :math:`a`
        if the following is true:

        .. math::

           c_i - \\frac{a}{2}\\le p_i\\le c_i + \\frac{a}{2}\\forall
           i\\in \\{x, y\\}

        """
        px, py = Point(point)
        cx, cy = self.center
        a = self.a
        xmin = cx - a / 2
        xmax = cx + a / 2
        ymin = cy - a / 2
        ymax = cy + a / 2
        return (px >= xmin) and (px <= xmax) and (py >= ymin) and (py <= ymax)

    def todict(self):
        """Returns a :class:`~python:dict` of the :class:`Square` \
            constructor parameters."""
        return dict(center=self.center, a=self.a)


class Triangle(Geometric2DRegion):
    """`Geometric2DRegion` for a triangle.

    .. versionadded:: 0.3.10

    Represents the bounded region with corner points
    :math:`p_1=(x_1,y_1)`, :math:`p_2=(x_2,y_2)`, and
    :math:`p_3=(x_3,y_3)`.

    Parameters
    ----------
    p1, p2, p3 : array_like, optional
        2-tuples or :class:`~sknano.core.Point` class instances
        specifying the `Triangle` corner points :math:`p_1=(x_1,y_1)`,
        :math:`p_2=(x_2,y_2)`, and :math:`p_3=(x_3, y_3)`.

    Notes
    -----
    :class:`Triangle` represents a 2D geometric region consisting of all
    combinations of corner points :math:`p_i`,
    :math:`\\left\\{\\lambda_1 p_1+\\lambda_2 p_2 + \\lambda_3 p_3|
    \\lambda_i\\ge0\\land\\lambda_1+\\lambda_2+\\lambda_3=1\\right\\}`.

    Calling :class:`Triangle` with no parameters is equivalent to
    :class:`Triangle`\ `(p1=[0, 0], p2=[0, 1], p3=[1, 0])`.

    """
    def __init__(self, p1=None, p2=None, p3=None):
        super().__init__()

        self._p1 = Point(nd=2)
        self._p2 = Point(nd=2)
        self._p3 = Point(nd=2)

        if p1 is None:
            p1 = [0, 0]
        if p2 is None:
            p2 = [0, 1]
        if p3 is None:
            p3 = [1, 0]

        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.points.extend([self.p1, self.p2, self.p3])
        self.fmtstr = "p1={p1!r}, p2={p2!r}, p3={p3!r}"

    @property
    def p1(self):
        """Corner point :math:`p_1=(x_1, y_1)`."""
        return self._p1

    @p1.setter
    def p1(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._p1[:] = Point(value)

    @property
    def p2(self):
        """Corner point :math:`p_2=(x_2, y_2)`."""
        return self._p2

    @p2.setter
    def p2(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._p2[:] = Point(value)

    @property
    def p3(self):
        """Corner point :math:`p_3=(x_3, y_3)`."""
        return self._p3

    @p3.setter
    def p3(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._p3[:] = Point(value)

    @property
    def area(self):
        """:class:`Triangle` area.

        Computed as:

        .. math::

           A = \\frac{1}{2}|-x_2 y_1 + x_3 y_1 + x_1 y_2 - x_3 y_2 - x_1 y_3
           + x_2 y_3|

        """
        x1, y1 = self.p1
        x2, y2 = self.p2
        x3, y3 = self.p3

        return np.abs(-x2 * y1 + x3 * y1 + x1 * y2 - x3 * y2 - x1 * y3 +
                      x2 * y3) / 2

    @property
    def centroid(self):
        """:class:`Triangle` centroid, :math:`(c_x, c_y)`.

        Computed as 2D :class:`Point` :math:`(c_x, c_y)` with
        coordinates:

        .. math::

           c_x = \\frac{x_1 + x_2 + x_3}{3}

           c_y = \\frac{y_1 + y_2 + y_3}{3}

        Returns
        -------
        :class:`Point`
            2D :class:`Point` of centroid.

        """
        x1, y1 = self.p1
        x2, y2 = self.p2
        x3, y3 = self.p3

        cx = (x1 + x2 + x3) / 3
        cy = (y1 + y2 + y3) / 3
        return Point([cx, cy])

    def contains(self, point):
        """Test region membership of `point` in :class:`Triangle`.

        Parameters
        ----------
        point : array_like

        Returns
        -------
        :class:`~python:bool`
            `True` if `point` is within :class:`Triangle`,
            `False`, otherwise.

        Notes
        -----
        A point :math:`(p_x, p_y)` is within the bounded region of a
        triangle with corner points :math:`p_1=(x_1, y_1)`,
        :math:`p_2=(x_2, y_2)`, and :math:`p_3=(x_3, y_3)`, if the
        following is true:

        .. math::

           \\frac{(x_1 - x_3) p_y + (x_3 - p_x) y_1 + (p_x - x_1) y_3}{
           (y_1-y_2) x_3 + (y_2 - y_3) x_1 + (y_3 - y_1) x_2}\\ge 0\\land

           \\frac{(x_2 - x_1) p_y + (p_x - x_2) y_1 + (x_1 - p_x) y_2}{
           (y_1 - y_2) x_3 + (y_2 - y_3) x_1 + (y_3 - y_1) x_2}\\ge 0\\land

           \\frac{(x_2 - x_3) p_y + (x_3 - p_x) y_2 + (p_x - x_2) y_3}{
           (y_1 - y_2) x_3 + (y_2 - y_3) x_1 + (y_3 - y_1) x_2}\\le 0


        """
        px, py = Point(point)
        x1, y1 = self.p1
        x2, y2 = self.p2
        x3, y3 = self.p3

        q1 = ((x1 - x3) * py + (x3 - px) * y1 + (px - x1) * y3) / \
            ((y1 - y2) * x3 + (y2 - y3) * x1 + (y3 - y1) * x2)

        q2 = ((x2 - x1) * py + (px - x2) * y1 + (x1 - px) * y2) / \
            ((y1 - y2) * x3 + (y2 - y3) * x1 + (y3 - y1) * x2)

        q3 = ((x2 - x3) * py + (x3 - px) * y2 + (px - x2) * y3) / \
            ((y1 - y2) * x3 + (y2 - y3) * x1 + (y3 - y1) * x2)

        return q1 >= 0 and q2 >= 0 and q3 <= 0

    def todict(self):
        """Returns a :class:`~python:dict` of the :class:`Triangle` \
            constructor parameters."""
        return dict(p1=self.p1, p2=self.p2, p3=self.p3)


class Ellipse(Geometric2DRegion):
    """`Geometric2DRegion` for an ellipse.

    .. versionadded:: 0.3.0

    Represents an axis-aligned ellipse centered at :math:`(c_x, c_y)`
    with semi-axes lengths :math:`r_x, r_y`.

    Parameters
    ----------
    center : array_like, optional
        Center of axis-aligned ellipse with semi-axes lengths :math:`r_x, r_y`
    rx, ry : float
        Lengths of semi-axes :math:`r_x, r_y`

    Notes
    -----
    :class:`Ellipse` represents the axis-aligned ellipsoid:

    .. math::

       \\left\\{\\{x, y, z\\}\\in R^3|
       \\left(\\frac{x-c_x}{r_x}\\right)^2 +
       \\left(\\frac{y-c_y}{r_y}\\right)^2 +
       \\left(\\frac{z-c_z}{r_z}\\right)^2\\le 1\\right\\}

    Calling :class:`Ellipse` with no parameters is equivalent to
    :class:`Ellipse`\ `(center=[0, 0], rx=1, ry=1)`.

    """
    def __init__(self, center=None, rx=1, ry=1):
        super().__init__()

        self._center = Point(nd=2)

        if center is None:
            center = [0, 0]
        self.center = center
        self.rx = rx
        self.ry = ry

        self.points.append(self.center)
        self.fmtstr = "center={center!r}, rx={rx:.3f}, ry={ry:.3f}"

    @property
    def center(self):
        """:class:`Ellipse` center point :math:`(c_x, c_y)`."""
        return self._center

    @center.setter
    def center(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._center[:] = Point(value)

    @property
    def rx(self):
        """Length of semi-axis :math:`r_x`."""
        return self._rx

    @rx.setter
    def rx(self, value):
        self._rx = float(value)

    @property
    def ry(self):
        """Length of semi-axis :math:`r_y`."""
        return self._ry

    @ry.setter
    def ry(self, value):
        self._ry = float(value)

    @property
    def area(self):
        """:class:`Ellipse` area, :math:`A=\\pi r_x r_y`."""
        return np.pi * self.rx * self.ry

    @property
    def centroid(self):
        """Alias for :attr:`~Ellipse.center`."""
        return self.center

    def contains(self, point):
        """Test region membership of `point` in :class:`Ellipse`.

        Parameters
        ----------
        point : array_like

        Returns
        -------
        :class:`~python:bool`
            `True` if `point` is within :class:`Ellipse`,
            `False` otherwise

        Notes
        -----
        A `point` :math:`(p_x, p_y)` is within the bounded region of
        an ellipse with center :math:`(c_x, c_y)` and semi-axes lengths
        :math:`r_x, r_y` if the following is true:

        .. math::

           \\left(\\frac{p_x - c_x}{r_x}\\right)^2 +
           \\left(\\frac{p_y - c_y}{r_y}\\right)^2\\le 1

        """
        px, py = Point(point)
        cx, cy = self.center
        rx, ry = self.rx, self.ry

        q1 = (px - cx) ** 2 / rx ** 2
        q2 = (py - cy) ** 2 / ry ** 2

        return q1 + q2 <= 1

    def todict(self):
        """Returns a :class:`~python:dict` of the :class:`Ellipse` \
            constructor parameters."""
        return dict(center=self.center, rx=self.rx, ry=self.ry)


class Circle(Geometric2DRegion):
    """`Geometric2DRegion` for a circle.

    .. versionadded:: 0.3.0

    Represents the bounded region with center :math:`(h, k)` and
    radius :math:`r`.

    Parameters
    ----------
    center : array_like, optional
        Center point :math:`(h, k)` of circle.
    r : float, optional
        Radius :math:`r` of circle.

    Notes
    -----
    Calling :class:`Circle` with no parameters is equivalent to
    :class:`Circle`\ `(center=[0, 0], r=1.0)`.

    """
    def __init__(self, center=None, r=1.0):

        super().__init__()

        self._center = Point(nd=2)

        if center is None:
            center = [0, 0]
        self.center = center
        self.r = r

        self.points.append(self.center)
        self.fmtstr = "center={center!r}, r={r:.3f}"

    @property
    def center(self):
        """Center point :math:`(h, k)` of circle."""
        return self._center

    @center.setter
    def center(self, value):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise TypeError('Expected a 2-element array_like object')
        self._center[:] = Point(value)

    @property
    def r(self):
        """:class:`Circle` radius, :math:`r`."""
        return self._r

    @r.setter
    def r(self, value):
        self._r = float(value)

    @property
    def area(self):
        """:class:`Circle` area, :math:`A=\\pi r^2`."""
        return np.pi * self.r ** 2

    @property
    def centroid(self):
        """Alias for :attr:`~Circle.center`."""
        return self.center

    def contains(self, point):
        """Test region membership of `point` in :class:`Circle`.

        Parameters
        ----------
        point : array_like

        Returns
        -------
        :class:`~python:bool`
            `True` if `point` is within :class:`Circle`,
            `False` otherwise.

        Notes
        -----
        A `point` :math:`(p_x, p_y)` is within the bounded region
        of a circle with center :math:`(h, k)` and radius :math:`r`
        if the following is true:

        .. math::

           (p_x - h)^2 + (p_y - k)^2 \\le r^2

        """
        x, y = Point(point)
        h, k = self.center
        r = self.r

        return (x - h) ** 2 + (y - k) ** 2 <= r ** 2

    def todict(self):
        """Returns a :class:`~python:dict` of the :class:`Circle` \
            constructor parameters."""
        return dict(center=self.center, r=self.r)
