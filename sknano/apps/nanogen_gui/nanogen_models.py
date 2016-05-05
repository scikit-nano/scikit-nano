# -*- coding: utf-8 -*-
"""
========================================================================
NanoGen model classes (:mod:`sknano.apps.nanogen_gui.nanogen_models`)
========================================================================

.. currentmodule:: sknano.apps.nanogen_gui.nanogen_models

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from sknano.core.refdata import aCC, element_data
from sknano.core.structures import compute_L, compute_T, \
    SWNT, MWNT, Graphene, UnrolledSWNT, Fullerene

_r_CC_vdw = element_data['C']['VanDerWaalsRadius']

__all__ = ['NanoGenModel', 'SWNTModel', 'MWNTModel',
           'GrapheneModel', 'FullereneModel',
           'CrystalStructureModel', 'ObserverModelMixin']


class ObserverModelMixin:
    """Observer model mixin for registering/removing/notifying observers."""
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
    def n(self):
        """Chiral index :math:`n`"""
        return self.structure.n

    @n.setter
    def n(self, value):
        self.structure.n = value
        self.notify_observers()

    @property
    def m(self):
        """Chiral index :math:`m`"""
        return self.structure.m

    @m.setter
    def m(self, value):
        self.structure.m = value
        self.notify_observers()

    @property
    def L(self):
        return compute_L(self.Ch, n3=self.n3, bond=self.bond)

    @L.setter
    def L(self, value):
        n3 = value / compute_T(self.Ch, bond=self.bond, length=True)
        if not self.structure.fix_L:
            n3 = int(n3)
        self.structure.n3 = n3
        self.notify_observers()

    @property
    def n3(self):
        return self.structure.n3

    @n3.setter
    def n3(self, value):
        self.structure.n3 = value
        self.notify_observers()

    @property
    def fix_L(self):
        return self.structure.fix_L

    @fix_L.setter
    def fix_L(self, value):
        self.structure.fix_L = value
        self.notify_observers()


class MWNTModelMixin:
    @property
    def L(self):
        return self.structure.L

    @L.setter
    def L(self, value):
        self.structure.L = value
        self.notify_observers()


class BundleModelMixin:
    @property
    def l1(self):
        # return compute_l1(self.Ch, n1=self.n1, bond=self.bond,
        #                   gutter=_r_CC_vdw)
        return self.structure.lattice.a1.length

    @property
    def l2(self):
        # return compute_l2(self.Ch, n2=self.n2, bond=self.bond,
        #                   gutter=_r_CC_vdw)
        return self.structure.lattice.a2.length

    @property
    def n1(self):
        return self.structure.n1

    @n1.setter
    def n1(self, value):
        self.structure.n1 = value
        self.notify_observers()

    @property
    def n2(self):
        return self.structure.n2

    @n2.setter
    def n2(self, value):
        self.structure.n2 = value
        self.notify_observers()


class SWNTModel(GeneratorModelBase, SWNTModelMixin, BundleModelMixin):
    def __init__(self):
        super().__init__()
        self.structure = SWNT((10, 10), basis=self.basis, bond=self.bond,
                              n1=1, n2=1, n3=1, bundle_packing='hcp')
        self.notify_observers()


class MWNTModel(GeneratorModelBase, MWNTModelMixin, BundleModelMixin):
    def __init__(self):
        super().__init__()
        self.structure = \
            MWNT(Ch_list=[(10, 10), (20, 20), (30, 30)], min_wall_diameter=5,
                 max_wall_diameter=50, wall_spacing=2 * _r_CC_vdw,
                 basis=self.basis, bond=self.bond, n1=1, n2=1, L=5,
                 bundle_packing='hcp')

        self.notify_observers()

    @property
    def Ch(self):
        return self.Ch_list[-1]

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

    @property
    def chiral_types(self):
        return self.structure.chiral_types

    @chiral_types.setter
    def chiral_types(self, value):
        self.structure.chiral_types = value
        self.notify_observers()

    @property
    def Nwalls(self):
        return self.structure.Nwalls

    @Nwalls.setter
    def Nwalls(self, value):
        self.structure.Nwalls = value
        self.notify_observers()

    @property
    def min_wall_diameter(self):
        return self.structure.min_wall_diameter

    @min_wall_diameter.setter
    def min_wall_diameter(self, value):
        self.structure.min_wall_diameter = value
        self.notify_observers()

    @property
    def max_wall_diameter(self):
        return self.structure.max_wall_diameter

    @max_wall_diameter.setter
    def max_wall_diameter(self, value):
        self.structure.max_wall_diameter = value
        self.notify_observers()

    @property
    def wall_spacing(self):
        return self.structure.wall_spacing

    @wall_spacing.setter
    def wall_spacing(self, value):
        self.structure.wall_spacing = value
        self.notify_observers()

    # def generate_Ch_list(self):
    #     self.structure = \
    #         MWNT(Nwalls=self.structure.Nwalls,
    #                    min_wall_diameter=self.structure.min_wall_diameter,
    #                    max_wall_diameter=self.structure.max_wall_diameter,
    #                    wall_spacing=self.structure.wall_spacing)


class UnrolledSWNTModelMixin(SWNTModelMixin):

    @property
    def n1(self):
        return self.structure.n1

    @n1.setter
    def n1(self, value):
        self.structure.n1 = value
        self.notify_observers()


class GrapheneModel(GeneratorModelBase, UnrolledSWNTModelMixin):

    def __init__(self):
        super().__init__()

        self.conventional_cell_graphene = \
            Graphene.from_conventional_cell(armchair_edge_length=1,
                                            zigzag_edge_length=1,
                                            basis=self.basis,
                                            bond=self.bond)
        self.primitive_cell_graphene = \
            Graphene.from_primitive_cell(edge_length=1,
                                         basis=self.basis,
                                         bond=self.bond)

        self.structure = UnrolledSWNT((10, 10), basis=self.basis,
                                      bond=self.bond, n1=1, n3=1)
        self.nlayers = 1
        self.layer_rotation_increment = 0.0
        self.notify_observers()

    @property
    def armchair_edge_length(self):
        return self.conventional_cell_graphene.armchair_edge_length

    @armchair_edge_length.setter
    def armchair_edge_length(self, value):
        self.conventional_cell_graphene.armchair_edge_length = value
        self.notify_observers()

    @property
    def zigzag_edge_length(self):
        return self.conventional_cell_graphene.zigzag_edge_length

    @zigzag_edge_length.setter
    def zigzag_edge_length(self, value):
        self.conventional_cell_graphene.zigzag_edge_length = value
        self.notify_observers()

    @property
    def edge_length(self):
        return self.primitive_cell_graphene.edge_length

    @edge_length.setter
    def edge_length(self, value):
        self.primitive_cell_graphene.edge_length = value
        self.notify_observers()

    @property
    def nlayers(self):
        return self._nlayers

    @nlayers.setter
    def nlayers(self, value):
        self._nlayers = value
        self.notify_observers()

    @property
    def layer_rotation_increment(self):
        return self._layer_rotation_increment

    @layer_rotation_increment.setter
    def layer_rotation_increment(self, value):
        self._layer_rotation_increment = value
        self.notify_observers()


class FullereneModel(ObserverModelMixin):
    def __init__(self):
        self._observers = []
        super().__init__()
        self.structure = Fullerene()
        self.notify_observers()

    @property
    def fullerenes(self):
        return self.structure.fullerenes


class CrystalStructureModel(ObserverModelMixin):
    pass
