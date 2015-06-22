# -*- coding: utf-8 -*-
"""
==========================================================
NanoGen model (:mod:`sknano.apps.nanogen_gui._ng_model`)
==========================================================

.. currentmodule:: sknano.apps.nanogen_gui._ng_model

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from sknano.core.refdata import aCC, dVDW
from sknano.structures import compute_Lx, compute_Ly, compute_Lz, \
    compute_T, SWNT, SWNTBundle, MWNT, MWNTBundle, Graphene, \
    UnrolledSWNT

__all__ = ['NGModel']


class NGModel:
    """:mod:`~sknano.apps.nanogen_gui` MVC model class."""
    def __init__(self):
        self._observers = []
        self._bond = aCC
        self._element1 = self._element2 = 'C'
        self._n = self._m = 10

        self.swnt = SWNT(self.Ch, basis=self.basis, bond=self.bond)
        self.swnt_bundle = SWNTBundle(**self.swnt.todict())

        self.mwnt = MWNT(Ch_list=None, Nwalls=3, min_wall_diameter=5,
                         max_wall_diameter=100, wall_spacing=dVDW,
                         basis=self.basis, bond=self.bond)
        self.mwnt_bundle = MWNTBundle(**self.mwnt.todict())

        self.graphene = Graphene(armchair_edge_length=1,
                                 zigzag_edge_length=1, basis=self.basis,
                                 bond=self.bond)
        self.unrolled_swnt = UnrolledSWNT(**self.swnt.todict())
        self.nlayers = 1
        self.layer_rotation_increment = 0.0

        self.notify_observers()

    @property
    def bond(self):
        return self._bond

    @bond.setter
    def bond(self, value):
        self._bond = value
        self.notify_observers()

    @property
    def element1(self):
        return self._element1

    @element1.setter
    def element1(self, value):
        self._element1 = value
        self.notify_observers()

    @property
    def element2(self):
        return self._element2

    @element2.setter
    def element2(self, value):
        self._element2 = value
        self.notify_observers()

    @property
    def Ch(self):
        """Chiral indices :math:`(n, m)`"""
        return self.n, self.m

    @property
    def basis(self):
        return [self.element1, self.element2]

    @property
    def n(self):
        """Chiral index :math:`n`"""
        return self._n

    @n.setter
    def n(self, value):
        self._n = value
        self.swnt.n = self.swnt_bundle.n = self.unrolled_swnt.n = self.n
        self.notify_observers()

    @property
    def m(self):
        """Chiral index :math:`m`"""
        return self._m

    @m.setter
    def m(self, value):
        self._m = value
        self.swnt.m = self.swnt_bundle.m = self.unrolled_swnt.m = self.m
        self.notify_observers()

    @property
    def Lx(self):
        return compute_Lx(self.Ch, nx=self.nx, bond=self.bond, gutter=dVDW)

    @property
    def Ly(self):
        return compute_Ly(self.Ch, ny=self.ny, bond=self.bond, gutter=dVDW)

    @property
    def Lz(self):
        return compute_Lz(self.Ch, nz=self.nz, bond=self.bond, gutter=dVDW)

    @Lz.setter
    def Lz(self, value):
        #self.swntbundle.Lz = value
        #self.notify_observers()
        # self._nz = \
        #     10 * value / compute_T(self.Ch, bond=self.bond, length=True)
        self.swnt.nz = self.swnt_bundle.nz = self.unrolled_swnt.nz = \
            10 * value / compute_T(self.Ch, bond=self.bond, length=True)

    @property
    def nx(self):
        return self._nx

    @nx.setter
    def nx(self, value):
        self._nx = value
        self.swnt_bundle.nx = self.mwnt_bundle.nx = self.unrolled_swnt.nx = \
            self.nx
        self.notify_observers()

    @property
    def ny(self):
        return self._ny

    @ny.setter
    def ny(self, value):
        self._ny = value
        self.swnt_bundle.ny = self.mwnt_bundle.ny = self.ny
        self.notify_observers()

    @property
    def nz(self):
        return self._nz

    @nz.setter
    def nz(self, value):
        self._nz = value
        self.swnt.nz = self.swnt_bundle.nz = self.unrolled_swnt.nz = self.nz
        self.notify_observers()

    @property
    def armchair_edge_length(self):
        return self.graphene.armchair_edge_length

    @armchair_edge_length.setter
    def armchair_edge_length(self, value):
        self.graphene.armchair_edge_length = value
        self.notify_observers()

    @property
    def zigzag_edge_length(self):
        return self.graphene.zigzag_edge_length

    @zigzag_edge_length.setter
    def zigzag_edge_length(self, value):
        self.graphene.zigzag_edge_length = value
        self.notify_observers()

    @property
    def nlayers(self):
        return self._nlayers

    @nlayers.setter
    def nlayers(self, value):
        self._nlayers = value
        self.graphene.nlayers = self.unrolled_swnt.nlayers = value
        self.notify_observers()

    def register_observer(self, observer):
        self._observers.append(observer)

    def remove_observer(self, observer):
        self._observers.remove(observer)

    def notify_observers(self):
        for observer in self._observers[:]:
            observer.update_app_view()
