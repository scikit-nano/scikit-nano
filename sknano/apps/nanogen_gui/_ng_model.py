# -*- coding: utf-8 -*-
"""
==========================================================
NanoGen model (:mod:`sknano.apps.nanogen_gui._ng_model`)
==========================================================

.. currentmodule:: sknano.apps.nanogen_gui._ng_model

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from sknano.core.refdata import CCbond, dVDW
from sknano.structures import SWNT, SWNTBundle, Graphene

__all__ = ['NGModel']


class NGModel(object):
    """:mod:`~sknano.apps.nanogen_gui` MVC model class."""
    def __init__(self):
        self._observers = []
        self.swnt = SWNT(n=10, m=10, bond=CCbond)
        self.swntbundle = SWNTBundle(n=10, m=10, bond=CCbond)
        self.graphene = Graphene(length=10, width=1, bond=CCbond, edge='AC')

    @property
    def n(self):
        """Chiral index :math:`n`"""
        return self.swntbundle.n

    @n.setter
    def n(self, value):
        self.swntbundle.n = value
        self.notify_observers()

    @property
    def m(self):
        """Chiral index :math:`m`"""
        return self.swntbundle.m

    @m.setter
    def m(self, value):
        self.swntbundle.m = value
        self.notify_observers()

    @property
    def nanotube_bond(self):
        return self.swntbundle.bond

    @nanotube_bond.setter
    def nanotube_bond(self, value):
        self.swntbundle.bond = value
        self.notify_observers()

    @property
    def nanotube_element1(self):
        return self.swntbundle.element1

    @nanotube_element1.setter
    def nanotube_element1(self, value):
        self.swntbundle.element1 = value
        self.notify_observers()

    @property
    def nanotube_element2(self):
        return self.swntbundle.element2

    @nanotube_element2.setter
    def nanotube_element2(self, value):
        self.swntbundle.element2 = value
        self.notify_observers()

    @property
    def Lx(self):
        return self.swntbundle.Lx

    @property
    def Ly(self):
        return self.swntbundle.Ly

    @property
    def Lz(self):
        return self.swntbundle.Lz

    @Lz.setter
    def Lz(self, value):
        #self.swntbundle.Lz = value
        #self.notify_observers()
        self.swntbundle.nz = value / self.swntbundle.T

    @property
    def nx(self):
        return self.swntbundle.nx

    @nx.setter
    def nx(self, value):
        self.swntbundle.nx = value
        self.notify_observers()

    @property
    def ny(self):
        return self.swntbundle.ny

    @ny.setter
    def ny(self, value):
        self.swntbundle.ny = value
        self.notify_observers()

    @property
    def nz(self):
        return self.swntbundle.nz

    @nz.setter
    def nz(self, value):
        self.swntbundle.nz = value
        self.notify_observers()

    @property
    def length(self):
        return self.graphene.length

    @length.setter
    def length(self, value):
        self.graphene.length = value
        self.notify_observers()

    @property
    def width(self):
        return self.graphene.width

    @width.setter
    def width(self, value):
        self.graphene.width = value
        self.notify_observers()

    @property
    def edge(self):
        return self.graphene.edge

    @edge.setter
    def edge(self, value):
        self.graphene.edge = value

    @property
    def nlayers(self):
        return self.graphene.nlayers

    @nlayers.setter
    def nlayers(self, value):
        self.graphene.nlayers = value
        self.notify_observers()

    @property
    def graphene_bond(self):
        return self.graphene.bond

    @graphene_bond.setter
    def graphene_bond(self, value):
        self.graphene.bond = value
        self.notify_observers()

    @property
    def graphene_element1(self):
        return self.graphene.element1

    @graphene_element1.setter
    def graphene_element1(self, value):
        self.graphene.element1 = value
        self.notify_observers()

    @property
    def graphene_element2(self):
        return self.graphene.element2

    @graphene_element2.setter
    def graphene_element2(self, value):
        self.graphene.element2 = value
        self.notify_observers()

    def register_observer(self, observer):
        self._observers.append(observer)

    def remove_observer(self, observer):
        self._observers.remove(observer)

    def notify_observers(self):
        self.swntbundle.compute_bundle_params()
        #self.graphene.compute_layer_params()
        for observer in self._observers[:]:
            observer.update_app_view()
