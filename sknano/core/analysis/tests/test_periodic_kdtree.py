# test_periodic_kdtree.py
#
# Unit tests for periodic_kdtree.py
#
# Written by Patrick Varilly, 6 Jul 2012
# Released under the scipy license

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np
from sknano.core.analysis import PeriodicKDTree

from nose.tools import assert_is_instance, assert_true
from numpy.testing import assert_equal, assert_array_equal, \
    assert_array_almost_equal, assert_almost_equal

from scipy.spatial import minkowski_distance


class ConsistencyTests:

    def distance(self, x, y, p):
        return minkowski_distance(np.zeros(x.shape), self.pbcs(x - y), p)

    def pbcs(self, x):
        # print('self.boxsize: {}'.format(self.boxsize))
        # print('x: {}'.format(x))
        val = x - np.where(self.boxsize > 0,
                           (np.round(x / self.boxsize) * self.boxsize), 0.0)
        # print('val: {}'.format(val))
        return val

    def test_nearest(self):
        x = self.x
        d, i = self.kdtree.query(x, 1)
        assert_almost_equal(d**2, np.sum(self.pbcs(x-self.data[i])**2))
        eps = 1e-8
        assert_true(np.all(np.sum(self.pbcs(self.data - x[np.newaxis, :])**2,
                           axis=1) > d**2 - eps))

    # def test_m_nearest(self):
    #     x = self.x
    #     m = self.m
    #     dd, ii = self.kdtree.query(x, m)
    #     d = np.amax(dd)
    #     print('d: {}'.format(d))
    #     print('np.argmax(dd): {}'.format(np.argmax(dd)))
    #     i = ii[np.argmax(dd)]
    #     print(i)
    #     assert_almost_equal(d**2, np.sum(self.pbcs(x-self.data[i])**2))
    #     eps = 1e-8
    #     assert_equal(np.sum(np.sum(self.pbcs(self.data-x[np.newaxis, :])**2,
    #                                axis=1) < d**2+eps), m)

    def test_points_near(self):
        x = self.x
        d = self.d
        dd, ii = self.kdtree.query(x, k=self.kdtree.n, distance_upper_bound=d)
        eps = 1e-8
        hits = 0
        for near_d, near_i in zip(dd, ii):
            if near_d == np.inf:
                continue
            hits += 1
            assert_almost_equal(near_d**2,
                                np.sum(self.pbcs(x-self.data[near_i])**2))
            assert_true(near_d < d + eps,
                        "near_d=%g should be less than %g" % (near_d, d))
        assert_equal(np.sum(np.sum(self.pbcs(self.data-x[np.newaxis, :])**2,
                                   axis=1) < d**2+eps), hits)

    def test_points_near_l1(self):
        x = self.x
        d = self.d
        dd, ii = self.kdtree.query(x, k=self.kdtree.n, p=1,
                                   distance_upper_bound=d)
        eps = 1e-8
        hits = 0
        for near_d, near_i in zip(dd, ii):
            if near_d == np.inf:
                continue
            hits += 1
            assert_almost_equal(near_d, self.distance(x, self.data[near_i], 1))
            assert_true(near_d < d+eps,
                        "near_d=%g should be less than %g" % (near_d, d))
        assert_equal(np.sum(self.distance(self.data, x, 1) < d+eps), hits)

    def test_points_near_linf(self):
        x = self.x
        d = self.d
        dd, ii = self.kdtree.query(x, k=self.kdtree.n, p=np.inf,
                                   distance_upper_bound=d)
        eps = 1e-8
        hits = 0
        for near_d, near_i in zip(dd, ii):
            # print('near_d: {}'.format(near_d))
            if near_d == np.inf:
                continue
            hits += 1
            assert_almost_equal(near_d,
                                self.distance(x, self.data[near_i], np.inf))
            assert_true(near_d < d+eps,
                        "near_d=%g should be less than %g" % (near_d, d))

        assert_equal(np.sum(self.distance(self.data, x, np.inf) < d+eps), hits)

    def test_approx(self):
        x = self.x
        k = self.k
        eps = 0.1
        d_real, i_real = self.kdtree.query(x, k)
        d, i = self.kdtree.query(x, k, eps=eps)
        assert_true(np.all(d <= d_real*(1+eps)))


class test_random(ConsistencyTests):
    def setUp(self):
        self.n = 100
        self.m = 4
        self.data = np.random.randn(self.n, self.m)
        self.boxsize = np.ones(self.m)
        self.kdtree = PeriodicKDTree(self.data, self.boxsize, leafsize=2)
        self.x = np.random.randn(self.m)
        self.d = 0.2
        self.k = 10


class test_random_far(test_random):
    def setUp(self):
        super().setUp()
        self.x = np.random.randn(self.m)+10


class test_small(ConsistencyTests):
    def setUp(self):
        self.data = np.array([[0, 0, 0],
                              [0, 0, 1],
                              [0, 1, 0],
                              [0, 1, 1],
                              [1, 0, 0],
                              [1, 0, 1],
                              [1, 1, 0],
                              [1, 1, 1]])
        self.boxsize = 1.1 * np.ones(3)
        self.kdtree = PeriodicKDTree(self.data, self.boxsize)
        self.n = self.kdtree.n
        self.m = self.kdtree.m
        self.x = np.random.randn(3)
        self.d = 0.5
        self.k = 4

    def test_nearest(self):
        assert_array_equal(self.kdtree.query((0, 0, 0.1), 1), (0.1, 0))

    def test_nearest_two(self):
        assert_array_almost_equal(self.kdtree.query((0, 0, 0.1), 2),
                                  ([0.1, np.sqrt(0.1**2 + 0.1**2)], [0, 2]))


class test_small_nonleaf(test_small):
    def setUp(self):
        super().setUp()
        self.kdtree = PeriodicKDTree(self.data, self.boxsize, leafsize=1)


class test_vectorization(object):
    def setUp(self):
        self.data = np.array([[0, 0, 0],
                              [0, 0, 1],
                              [0, 1, 0],
                              [0, 1, 1],
                              [1, 0, 0],
                              [1, 0, 1],
                              [1, 1, 0],
                              [1, 1, 1]])
        self.boxsize = 1.1 * np.ones(3)
        self.kdtree = PeriodicKDTree(self.data, self.boxsize)

    def test_single_query(self):
        d, i = self.kdtree.query(np.array([0, 0, 0]))
        assert_is_instance(d, float)
        assert_true(np.issubdtype(i, int))

    def test_vectorized_query(self):
        d, i = self.kdtree.query(np.zeros((2, 4, 3)))
        assert_equal(np.shape(d), (2, 4))
        assert_equal(np.shape(i), (2, 4))

    def test_single_query_multiple_neighbors(self):
        s = 23
        kk = 27*self.kdtree.n+s
        d, i = self.kdtree.query(np.array([0, 0, 0]), k=kk)
        assert_equal(np.shape(d), (kk,))
        assert_equal(np.shape(i), (kk,))
        assert_true(np.all(~np.isfinite(d[-s:])))
        assert_true(np.all(i[-s:] == self.kdtree.n))

    def test_vectorized_query_multiple_neighbors(self):
        s = 23
        kk = 27*self.kdtree.n+s
        d, i = self.kdtree.query(np.zeros((2, 4, 3)), k=kk)
        assert_equal(np.shape(d), (2, 4, kk))
        assert_equal(np.shape(i), (2, 4, kk))
        assert_true(np.all(~np.isfinite(d[:, :, -s:])))
        assert_true(np.all(i[:, :, -s:] == self.kdtree.n))

    def test_single_query_all_neighbors(self):
        d, i = self.kdtree.query([0, 0, 0], k=None, distance_upper_bound=1.1)
        assert_is_instance(d, list)
        assert_is_instance(i, list)

    def test_vectorized_query_all_neighbors(self):
        d, i = self.kdtree.query(np.zeros((2, 4, 3)), k=None,
                                 distance_upper_bound=1.1)
        assert_equal(np.shape(d), (2, 4))
        assert_equal(np.shape(i), (2, 4))

        assert_is_instance(d[0, 0], list)
        assert_is_instance(i[0, 0], list)


class ball_consistency(object):
    def distance(self, x, y, p):
        print(self.pbcs(x - y))
        return minkowski_distance(np.zeros(x.shape), self.pbcs(x - y), p)

    def pbcs(self, x):
        return x - np.where(self.boxsize > 0,
                            (np.round(x / self.boxsize) * self.boxsize), 0.0)

    def test_in_ball(self):
        l = self.T.query_ball_point(self.x, self.d, p=self.p, eps=self.eps)
        for i in l:
            assert_true(self.distance(self.data[i], self.x, self.p) <=
                        self.d * (1.+self.eps))

    # def test_found_all(self):
    #     c = np.ones(self.T.n, dtype=bool)
    #     l = self.T.query_ball_point(self.x, self.d, p=self.p, eps=self.eps)
    #     print(l)
    #     c[l] = False
    #     # print(self.distance(self.data[c], self.x, self.p) <= self.d * (1.+self.eps))
    #     assert_true(np.all(self.distance(self.data[c], self.x, self.p) <=
    #                        self.d * (1.+self.eps)))


class test_random_ball(ball_consistency):

    def setUp(self):
        n = 100
        m = 4
        self.data = np.random.randn(n, m)
        self.boxsize = np.ones(m)
        print(self.boxsize)
        self.T = PeriodicKDTree(self.data, self.boxsize, leafsize=10)
        self.x = np.random.randn(m)
        print(self.x)
        self.p = 2.
        self.eps = 0
        self.d = 0.2


class test_random_ball_approx(test_random_ball):

    def setUp(self):
        super().setUp()
        self.eps = 0.1


class test_random_ball_far(test_random_ball):

    def setUp(self):
        super().setUp()
        self.d = 2.


class test_random_ball_l1(test_random_ball):

    def setUp(self):
        super().setUp()
        self.p = 1


class test_random_ball_linf(test_random_ball):

    def setUp(self):
        super().setUp()
        self.p = np.inf


def test_random_ball_vectorized():

    n = 20
    m = 5
    boxsize = np.ones(m)
    T = PeriodicKDTree(np.random.randn(n, m), boxsize)

    r = T.query_ball_point(np.random.randn(2, 3, m), 1)
    assert_equal(r.shape, (2, 3))
    assert_is_instance(r[0, 0], list)
