# -*- coding: utf-8 -*-
"""
========================================================================
Stopwatch Timer class (:mod:`sknano.testing._timer`)
========================================================================

.. currentmodule:: sknano.testing._timer

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

__all__ = ['Timer']


import time


class Timer:
    """Stopwatch Timer class."""
    def __init__(self, func=time.perf_counter):
        self.elapsed = 0.0
        self._func = func
        self._start = None

    def start(self):
        if self._start is not None:
            raise RuntimeError('Already started')
        self._start = self._func()

    def stop(self):
        if self._start is None:
            raise RuntimeError('Not started')
        end = self._func()
        self.elapsed += end - self._start
        self._start = None

    def reset(self):
        self.elapsed = 0.0

    @property
    def running(self):
        return self._start is not None

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, *args):
        self.stop()
