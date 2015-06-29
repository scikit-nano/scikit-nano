# -*- coding: utf-8 -*-
"""
========================================================================
NanoGen model classes (:mod:`sknano.apps.nanogen_gui._nanogen_models`)
========================================================================

.. currentmodule:: sknano.apps.nanogen_gui._nanogen_models

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from sknano.core.refdata import aCC, dVDW
from sknano.structures import compute_Lx, compute_Ly, compute_Lz, \
    compute_Ch, compute_T, SWNT, SWNTBundle, MWNT, MWNTBundle, Graphene, \
    UnrolledSWNT

__all__ = ['NanoGenModel', 'SWNTModel', 'MWNTModel',
           'GrapheneModel', 'FullereneModel',
           'BulkStructureModel', 'ObserverModelMixin']


class ObserverModelMixin:
    def register_observer(self, observer):
        self._observers.append(observer)

    def remove_observer(self, observer):
        self._observers.remove(observer)

    def notify_observers(self):
        for observer in self._observers[:]:
            observer.update_app_view()


class NanoGenModel(ObserverModelMixin):
    """:mod:`~sknano.apps.nanogen_gui` MVC model class."""
    def __init__(self):
        self._observers = []
        self.notify_observers()


class GeneratorModelBase(ObserverModelMixin):

    def __init__(self):
        self._observers = []
        self._bond = aCC
        self._element1 = self._element2 = 'C'

    @property
    def basis(self):
        return [self.element1, self.element2]

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


class SWNTModelMixin:
    @property
    def Ch(self):
        """Chiral indices :math:`(n, m)`"""
        return self.n, self.m

    @property
    def Lz(self):
        return compute_Lz(self.Ch, nz=self.nz, bond=self.bond, gutter=dVDW)

    @Lz.setter
    def Lz(self, value):
        # self._nz = \
        #     10 * value / compute_T(self.Ch, bond=self.bond, length=True)
        # self.structure.nz = self.nz
        # if self.structure.fix_Lz:
        nz = 10 * value / compute_T(self.Ch, bond=self.bond, length=True)
        if not self.structure.fix_Lz:
            nz = int(nz)
        self.structure.nz = nz
        self.notify_observers()

    @property
    def nz(self):
        # return self._nz
        return self.structure.nz

    @nz.setter
    def nz(self, value):
        # self._nz = value
        # self.structure.nz = self.nz
        self.structure.nz = value
        self.notify_observers()

    @property
    def fix_Lz(self):
        return self.structure.fix_Lz

    @fix_Lz.setter
    def fix_Lz(self, value):
        self.structure.fix_Lz = value


class MWNTModelMixin:
    @property
    def Lz(self):
        return self.structure.Lz

    @Lz.setter
    def Lz(self, value):
        # self._Lz = value
        # self.structure.Lz = self.Lz
        self.structure.Lz = value
        self.notify_observers()


class BundleModelMixin:
    @property
    def Lx(self):
        return compute_Lx(self.Ch, nx=self.nx, bond=self.bond, gutter=dVDW)

    @property
    def Ly(self):
        return compute_Ly(self.Ch, ny=self.ny, bond=self.bond, gutter=dVDW)

    @property
    def nx(self):
        # return self._nx
        return self.structure.nx

    @nx.setter
    def nx(self, value):
        # self._nx = value
        # self.structure.nx = self.nx
        self.structure.nx = value
        self.notify_observers()

    @property
    def ny(self):
        # return self._ny
        return self.structure.ny

    @ny.setter
    def ny(self, value):
        # self._ny = value
        # self.structure.ny = self.ny
        self.structure.ny = value
        self.notify_observers()


class SWNTModel(GeneratorModelBase, SWNTModelMixin, BundleModelMixin):
    def __init__(self):
        super().__init__()
        # self._n = self._m = 10
        self.structure = SWNTBundle((10, 10), basis=self.basis, bond=self.bond,
                                    nx=1, ny=1, nz=1, bundle_packing='hcp')
        self.notify_observers()

    @property
    def n(self):
        """Chiral index :math:`n`"""
        # return self._n
        return self.structure.n

    @n.setter
    def n(self, value):
        # self._n = value
        # self.structure.n = self.n
        self.structure.n = value
        self.notify_observers()

    @property
    def m(self):
        """Chiral index :math:`m`"""
        # return self._m
        return self.structure.m

    @m.setter
    def m(self, value):
        # self._m = value
        # self.structure.m = self.m
        self.structure.m = value
        self.notify_observers()


class MWNTModel(GeneratorModelBase, MWNTModelMixin, BundleModelMixin):
    def __init__(self):
        super().__init__()
        self.structure = \
            MWNTBundle(Ch_list=None, Nwalls=3, min_wall_diameter=5,
                       max_wall_diameter=100, wall_spacing=dVDW,
                       basis=self.basis, bond=self.bond, nx=1, ny=1, Lz=1,
                       bundle_packing='hcp')

        self.notify_observers()

    @property
    def Ch_list(self):
        # return self._Ch_list
        return self.structure.Ch_list

    @Ch_list.setter
    def Ch_list(self, value):
        # self._Ch_list = value
        # self.structure.Ch_list = self.Ch_list
        self.structure.Ch_list = value
        self.notify_observers()


class UnrolledSWNTModelMixin(SWNTModelMixin):
    @property
    def Lx(self):
        # return self.nx * compute_Ch(self.Ch, bond=self.bond) / 10
        return self.structure.Lx

    @Lx.setter
    def Lx(self, value):
        # self._nx = 10 * value / compute_Ch(self.Ch, bond=self.bond)
        # self.structure.nx = self.nx
        self.structure.nx = 10 * value / compute_Ch(self.Ch, bond=self.bond)
        self.notify_observers()

    @property
    def nx(self):
        # return self._nx
        return self.structure.nx

    @nx.setter
    def nx(self, value):
        # self._nx = value
        # self.structure.nx = self.nx
        self.structure.nx = value
        self.notify_observers()


class GrapheneModel(GeneratorModelBase, UnrolledSWNTModelMixin):

    def __init__(self):
        super().__init__()

        self.graphene = Graphene(armchair_edge_length=1,
                                 zigzag_edge_length=1, basis=self.basis,
                                 bond=self.bond)
        # self._n = self._m = 10
        self.structure = UnrolledSWNT((10, 10), basis=self.basis,
                                      bond=self.bond, nx=1, nz=1)
        self.nlayers = 1
        self.layer_rotation_increment = 0.0
        self.notify_observers()

    @property
    def armchair_edge_length(self):
        return self._armchair_edge_length

    @armchair_edge_length.setter
    def armchair_edge_length(self, value):
        self._armchair_edge_length = value
        self.notify_observers()

    @property
    def zigzag_edge_length(self):
        return self._zigzag_edge_length

    @zigzag_edge_length.setter
    def zigzag_edge_length(self, value):
        self._zigzag_edge_length = value
        self.notify_observers()

    @property
    def nlayers(self):
        return self._nlayers

    @nlayers.setter
    def nlayers(self, value):
        self._nlayers = value
        # self.graphene.nlayers = self.structure.nlayers = self.nlayers
        self.notify_observers()

    @property
    def layer_rotation_increment(self):
        return self._layer_rotation_increment

    @layer_rotation_increment.setter
    def layer_rotation_increment(self, value):
        self._layer_rotation_increment = value
        self.notify_observers()


class FullereneModel(ObserverModelMixin):
    pass


class BulkStructureModel(ObserverModelMixin):
    pass
