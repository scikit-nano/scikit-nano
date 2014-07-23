# -*- coding: utf-8 -*-
"""
====================================================
NanoGen model (:mod:`sknano.nanogen_gui._ng_model`)
====================================================

.. currentmodule:: sknano.nanogen_gui._ng_model

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ..nanogen import NanotubeBundle, Graphene

__all__ = ['NGModel']


class NGModel(object):
    """:mod:`~sknano.nanogen_gui` MVC model class."""
    def __init__(self):
        self._observers = []
        self._swntbundle = NanotubeBundle(n=10, m=10, bond=1.42)
        self._graphene = Graphene(length=10, width=1, bond=1.42, edge='AC')

    @property
    def n(self):
        """Chiral index :math:`n`"""
        return self._swntbundle.n

    @n.setter
    def n(self, value):
        self._swntbundle.n = value
        self.notify_observers()

    @property
    def m(self):
        """Chiral index :math:`m`"""
        return self._swntbundle.m

    @m.setter
    def m(self, value):
        self._swntbundle.m = value
        self.notify_observers()

    @property
    def nanotube_bond(self):
        return self._swntbundle.bond

    @nanotube_bond.setter
    def nanotube_bond(self, value):
        self._swntbundle.bond = value
        self.notify_observers()

    @property
    def nanotube_element1(self):
        return self._swntbundle.element1

    @nanotube_element1.setter
    def nanotube_element1(self, value):
        self._swntbundle.element1 = value
        self.notify_observers()

    @property
    def nanotube_element2(self):
        return self._swntbundle.element2

    @nanotube_element2.setter
    def nanotube_element2(self, value):
        self._swntbundle.element2 = value
        self.notify_observers()

    @property
    def Lx(self):
        return self._swntbundle.Lx

    @property
    def Ly(self):
        return self._swntbundle.Ly

    @property
    def Lz(self):
        return self._swntbundle.Lz

    @Lz.setter
    def Lz(self, value):
        self._swntbundle.Lz = value
        self.notify_observers()

    @property
    def nx(self):
        return self._swntbundle.nx

    @nx.setter
    def nx(self, value):
        self._swntbundle.nx = value
        self.notify_observers()

    @property
    def ny(self):
        return self._swntbundle.ny

    @ny.setter
    def ny(self, value):
        self._swntbundle.ny = value
        self.notify_observers()

    @property
    def nz(self):
        return self._swntbundle.nz

    @nz.setter
    def nz(self, value):
        self._swntbundle.nz = value
        self.notify_observers()

    @property
    def length(self):
        return self._graphene.length

    @length.setter
    def length(self, value):
        self._graphene.length = value
        self.notify_observers()

    @property
    def width(self):
        return self._graphene.width

    @width.setter
    def width(self, value):
        self._graphene.width = value
        self.notify_observers()

    @property
    def edge(self):
        return self._graphene.edge

    @edge.setter
    def edge(self, value):
        self._graphene.edge = value

    @property
    def nlayers(self):
        return self._graphene.nlayers

    @nlayers.setter
    def nlayers(self, value):
        self._graphene.nlayers = value
        self.notify_observers()

    @property
    def graphene_bond(self):
        return self._graphene.bond

    @graphene_bond.setter
    def graphene_bond(self, value):
        self._graphene.bond = value
        self.notify_observers()

    @property
    def graphene_element1(self):
        return self._graphene.element1

    @graphene_element1.setter
    def graphene_element1(self, value):
        self._graphene.element1 = value
        self.notify_observers()

    @property
    def graphene_element2(self):
        return self._graphene.element2

    @graphene_element2.setter
    def graphene_element2(self, value):
        self._graphene.element2 = value
        self.notify_observers()

    def register_observer(self, observer):
        self._observers.append(observer)

    def remove_observer(self, observer):
        self._observers.remove(observer)

    def notify_observers(self):
        self._swntbundle.compute_bundle_params()
        self._graphene.compute_layer_params()
        for observer in self._observers[:]:
            observer.update_app_view()
