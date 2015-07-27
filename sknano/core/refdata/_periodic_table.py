# -*- coding: utf-8 -*-
"""
==========================================================================
Periodic Table of Elements (:mod:`sknano.core.refdata._periodic_table`)
==========================================================================

.. currentmodule:: sknano.core.refdata._periodic_table

.. autosummary::
   :toctree: generated/

"""
from __future__ import division, absolute_import, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from collections import OrderedDict

from ._element_data import element_data
# from ._nist import atomic_masses_nist_list as _atomic_masses

__all__ = ['atomic_masses', 'atomic_mass_symbol_map',
           'atomic_numbers', 'atomic_number_symbol_map',
           'element_symbols', 'element_names']

element_symbols = \
    ['H', 'He',  # Period 1
     'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',  # Period 2
     'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',  # Period 3
     'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',  # Period 4
     'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',  # Period 4
     'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',  # Period 5
     'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',  # Period 5
     'Cs', 'Ba',  # Period 6
     'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',  # Lanthanides
     'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',  # Lanthanides
     'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',  # Period 6
     'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',  # Period 6
     'Fr', 'Ra',  # Period 7
     'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',  # Actinides
     'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',  # Actinides
     'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',  # Period 7
     'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo']  # Period 7

element_names = [
    'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron',
    'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium',
    'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur',
    'Chlorine', 'Argon', 'Potassium', 'Calcium', 'Scandium',
    'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron',
    'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium', 'Germanium',
    'Arsenic', 'Selenium', 'Bromine', 'Krypton', 'Rubidium',
    'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum',
    'Technetium', 'Ruthenium', 'Rhodium', 'Palladium', 'Silver',
    'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium',
    'Iodine', 'Xenon', 'Caesium', 'Barium', 'Lanthanum',
    'Cerium', 'Praseodymium', 'Neodymium', 'Promethium',
    'Samarium', 'Europium', 'Gadolinium', 'Terbium',
    'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium',
    'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium',
    'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury',
    'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine',
    'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium',
    'Protactinium', 'Uranium', 'Neptunium', 'Plutonium',
    'Americium', 'Curium', 'Berkelium', 'Californium',
    'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium',
    'Lawrencium', 'Rutherfordium', 'Dubnium', 'Seaborgium',
    'Bohrium', 'Hasium', 'Meitnerium', 'Darmstadtium',
    'Roentgenium', 'Copernicium', 'Ununtrium', 'Flerovium',
    'Ununpentium', 'Livermorium', 'Ununseptium', 'Ununoctium']

atomic_numbers = OrderedDict()
[atomic_numbers.update({symbol: element_data[symbol]['AtomicNumber']})
 for symbol in element_symbols]

atomic_number_symbol_map = {v: k for k, v in list(atomic_numbers.items())}

atomic_masses = OrderedDict()
# [atomic_masses.update({symbol: float(_atomic_masses[i])})
#  if _atomic_masses[i] is not None else
#  atomic_masses.update({symbol: None}) for i, symbol in
#  enumerate(element_symbols)]
[atomic_masses.update({symbol: element_data[symbol]['AtomicMass']})
 for symbol in element_symbols]

atomic_mass_symbol_map = {v: k for k, v in list(atomic_masses.items()) if
                          v is not None}


class PeriodicTable:
    """Class for creating abstract object representing Periodic Table."""

    def __init__(self):
        pass
