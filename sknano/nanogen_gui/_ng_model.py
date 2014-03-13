# -*- coding: utf-8 -*-
"""
====================================================
NanoGen model (:mod:`sknano.nanogen_gui._ng_model`)
====================================================

.. currentmodule:: sknano.nanogen_gui._ng_model

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ..nanogen import NanotubeBundle

__all__ = ['NGModel']


class NGModel(NanotubeBundle):
    """:py:mod:`~sknano.nanogen_gui` MVC model class."""
    def __init__(self):

        self.observers = []
        super(NGModel, self).__init__(n=10, m=10, verbose=True)

    def init(self):
        self.notify_observers()

    def _compute_bundle_params(self):
        super(NGModel, self).compute_bundle_params()
        self.notify_observers()

    @property
    def n(self):
        """Chiral index :math:`n`"""
        return self.n

    @n.setter
    def n(self, value):
        self.n = value
        self._compute_bundle_params()

    @property
    def m(self):
        """Chiral index :math:`m`"""
        return self.m

    @m.setter
    def m(self, value):
        self.m = value
        self._compute_bundle_params()

    #@property
    #def bond(self):
    #    return self.bond

    #@bond.setter
    #def bond(self, value):
    #    self.bond = value
    #    self._compute_bundle_params()

    @property
    def Lz(self):
        return self.Lz

    @Lz.setter
    def Lz(self, value):
        self.Lz = value
        self._compute_bundle_params()

    @property
    def nz(self):
        return self.nz

    @nz.setter
    def nz(self, value):
        self.nz = value
        self._compute_bundle_params()

    def register_observer(self, observer):
        self.observers.append(observer)

    def remove_observer(self, observer):
        self.observers.remove(observer)

    def notify_observers(self):
        for observer in self.observers[:]:
            observer.update_app_view()
