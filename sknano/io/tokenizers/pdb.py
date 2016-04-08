# -*- coding: utf-8 -*-
"""
=======================================================
PDB tokenizer class (:mod:`sknano.io.tokenizers.pdb`)
=======================================================

.. currentmodule:: sknano.io.tokenizers.pdb

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import OrderedDict

from pyparsing import Combine, Empty, Forward, Group, Keyword, Literal, \
    OneOrMore, Optional, Word, White, Suppress, Regex, ZeroOrMore, \
    alphas, alphanums, oneOf, nums, matchOnlyAtCol, restOfLine

from sknano.core import integer, real, asint
from sknano.core.crystallography import Crystal3DLattice
from sknano.core.refdata import element_symbols
from sknano.core.atoms import StructureAtom as Atom

__all__ = ['PDBTokenizer']


records = OrderedDict()

records['HEADER'] = header = {}
records['OBSLTE'] = obslte = {}
records['TITLE'] = title = {}
title['format'] = ''

records['SPLIT'] = split = {}
records['CAVEAT'] = caveat = {}
records['COMPND'] = compnd = {}
records['SOURCE'] = source = {}
records['KEYWDS'] = keywds = {}
records['EXPDTA'] = expdta = {}
records['NUMMDL'] = nummdl = {}
records['MDLTYP'] = mdltyp = {}
records['AUTHOR'] = author = {}
records['REVDAT'] = revdat = {}
records['SPRSDE'] = sprsde = {}
records['JRNL'] = jrnl = {}
records['REMARK 0'] = remark_0 = {}
records['REMARK 1'] = remark_1 = {}
records['REMARK 2'] = remark_2 = {}
records['REMARK 3'] = remark_3 = {}
records['REMARK N'] = remark_N = {}
records['DBREF'] = dbref = {}
records['DBREF1'] = dbref1 = {}
records['DBREF2'] = dbref2 = {}
records['SEQADV'] = seqadv = {}
records['SEQRES'] = seqres = {}
records['MODRES'] = modres = {}
records['HET'] = het = {}
records['HETNAM'] = hetnam = {}
records['HETSYN'] = hetsyn = {}
records['FORMUL'] = formul = {}
records['HELIX'] = helix = {}
records['SHEET'] = sheet = {}
records['SSBOND'] = ssbond = {}
records['LINK'] = link = {}
records['CISPEP'] = cispep = {}
records['SITE'] = site = {}
records['CRYST1'] = cryst1 = {}
records['ORIGX1'] = origx1 = {}
records['ORIGX2'] = origx2 = {}
records['ORIGX3'] = origx3 = {}
records['SCALE1'] = scale1 = {}
records['SCALE2'] = scale2 = {}
records['SCALE3'] = scale3 = {}
records['MTRIX1'] = mtrix1 = {}
records['MTRIX2'] = mtrix2 = {}
records['MTRIX3'] = mtrix3 = {}
records['MODEL'] = model = {}
records['ATOM'] = atom = {}

atom['format'] = \
    '{:6}{:>5} {:>4}{:1}{:>3} {:1}{:>4}{:1}' + \
    '{:3}'.format('') + \
    '{:>8.3f}{:>8.3f}{:>8.3f}' + \
    '{:>6.2f}{:>6.2f}' + \
    '{:10}'.format('') + \
    '{:2}{:2}'

records['ANISOU'] = {}
records['TER'] = {}
records['HETATM'] = {}
records['ENDMDL'] = {}
records['CONECT'] = conect = {}
records['MASTER'] = master = {}
records['END'] = None

DASH = Literal('-') | Literal('-')
SPACE = Literal(' ')
BACKSLASH = Literal('\\')
COMMA = Literal(',')
COLON = Literal(':')
SEMICOLON = Literal(';')
NULL = Literal('NULL')

pdbcharset = alphanums + " `-=[]'./~!@#$%^&*()_+{}|\"<>?"

# classification = Word(alphanums, alphanums + '/', max=40)

special_chars = ',:;'

# def string_field(exact=0):
#     return Word(pdbcharset, exact=exact)


def string_field_with_special_chars(chars, exact=0):
    return Word(pdbcharset + chars, exact=exact)

string_field = Word(pdbcharset)

# csv_list = Group(delimitedList(string_field, delim=','))
csv_list = Group(string_field + ZeroOrMore(Suppress(COMMA) + string_field))
slist = Group(string_field + ZeroOrMore(Suppress(SEMICOLON) + string_field))

lstring = string_field_with_special_chars(',:;') + \
    ZeroOrMore(OneOrMore(White(' ')) + string_field)

# slist = Group(delimitedList(string_field, delim=';'))
# lstring = Group(delimitedList(string_field, delim=' '))


def padexpr(expr, padwidth):
    return expr + nspaces(padwidth)


def nspaces(n):
    return Suppress(White(' ', exact=n))


def character_field(width=None, ci=None, cf=None, padwidth=None):
    if width is None:
        width = (cf - ci) + 1
    if padwidth is None:
        padwidth = 1

    field = nspaces(width)
    padding = nspaces(padwidth)

    if width == 1:
        field |= Word(alphanums + ' ', exact=width)
    else:
        n = 0
        while width > 0:
            width -= 1
            n += 1
            field |= (nspaces(width) + Word(alphanums + ' ', exact=n))
    return field + padding


def integer_field(width=None, ci=None, cf=None, padwidth=None):
    if width is None:
        width = (cf - ci) + 1
    if padwidth is None:
        padwidth = 1

    field = nspaces(width)
    padding = nspaces(padwidth)
    n = 0
    while width > 0:
        width -= 1
        n += 1
        field |= (nspaces(width) + Word(nums, exact=n).setParseAction(asint))
    if padwidth == 0:
        return field
    else:
        return field + padding


def lstring_field(width=None, ci=None, cf=None, padwidth=None):
    if width is None:
        width = (cf - ci) + 1
    if padwidth is None:
        padwidth = 1

    field = nspaces(width)
    padding = nspaces(padwidth)
    n = 0
    while width > 0:
        width -= 1
        n += 1
        field |= (nspaces(width) +
                  string_field_with_special_chars(',:;', exact=n))
    if padwidth == 0:
        return field
    else:
        return field + padding


def continuation_field(width):
    field = nspaces(width)
    n = 0
    while width > 0:
        width -= 1
        n += 1
        field |= (nspaces(width) + Word(nums, exact=n).addParseAction(asint))
    return field


def Record(name, width=None, ci=None, cf=None, padwidth=None):
    if width is None:
        if ci is not None and cf is not None:
            width = (cf - ci) + 1
        else:
            width = len(name)
    if padwidth is None:
        padwidth = 0
    if padwidth != 0:
        return Keyword(name) + nspaces(padwidth)
    else:
        return Keyword(name)


def new_atom(s, l, t):
    element = getattr(t, 'element', getattr(t, 'name', 'X'))
    params = dict(element=element, serial=t.serial, x=t.x, y=t.y, z=t.z)
    return Atom(**params)


def new_xtal_lattice(s, l, t):
    params = dict(a=t.a, b=t.b, c=t.c,
                  alpha=t.alpha, beta=t.beta, gamma=t.gamma)
    return Crystal3DLattice(**params)


day = Word(nums, exact=2)
month = Word(alphas, exact=3)
year = Word(nums, exact=2)
date_field = Combine(day + DASH + month + DASH + year)

id_code = Word(nums, alphanums, exact=4)
token = Word(alphas, alphanums + '_-() ') + Literal(':') + SPACE
token_value = Word(alphanums, alphanums + '_-() ')
speclist = Combine(OneOrMore(token + token_value) + Optional(SEMICOLON))

header_expr = \
    Record("HEADER", padwidth=4) + string_field + \
    Optional(date_field + id_code)
title_expr = \
    Record("TITLE ", padwidth=2) + continuation_field(2) + string_field
compnd_expr = Record("COMPND", padwidth=1) + continuation_field(3) + speclist
source_expr = Record("SOURCE", padwidth=1) + continuation_field(3) + speclist
keywds_expr = Record("KEYWDS", padwidth=2) + continuation_field(2) + csv_list
expdta_expr = Record("EXPDTA", padwidth=2) + continuation_field(2) + slist
author_expr = Record("AUTHOR", padwidth=2) + continuation_field(2) + csv_list
revdat_expr = Record("REVDAT", padwidth=1) + Word(nums, min=1, max=3) + \
    (continuation_field(2) + date_field | date_field) + \
    id_code + Word(nums, exact=1) + ZeroOrMore(lstring)

#     Optional(date_field("modDate") + Suppress(SPACE) + id_code("modId")) + \
#     Suppress(OneOrMore(SPACE)) + Word(nums, max=1)("modType") + \
#     Suppress(OneOrMore(SPACE)) + \
#     Optional(delimitedList(Word(alphas), delim=' ')) + StringEnd()

jrnl_records = Forward()
jrnl_expr = Record("JRNL", width=6, padwidth=6) + jrnl_records


jrnl_auth_expr = Keyword('AUTH') + continuation_field(2) + csv_list
jrnl_titl_expr = Keyword('TITL') + continuation_field(2) + lstring
jrnl_edit_expr = Keyword('EDIT') + continuation_field(2) + lstring
jrnl_ref_expr = Keyword('REF') + lstring
jrnl_publ_expr = Keyword('PUBL') + continuation_field(2) + lstring
jrnl_refn_expr = Keyword('REFN') + (Keyword("ISSN") | Keyword("ESSN")) + \
    lstring
jrnl_pmid_expr = Keyword('PMID') + integer
jrnl_doi_expr = Keyword('DOI') + lstring

jrnl_records << \
    (jrnl_auth_expr | jrnl_titl_expr | jrnl_edit_expr | jrnl_ref_expr |
     jrnl_publ_expr | jrnl_refn_expr | jrnl_pmid_expr | jrnl_doi_expr)

remark_expr = Record("REMARK", padwidth=1) + integer("remarkNum") + \
    (nspaces(2) + jrnl_records | nspaces(1) + lstring | Empty())

cryst_expr = Record("CRYST1") + \
    real("a") + real("b") + real("c") + \
    real("alpha") + real("beta") + real("gamma") + \
    Group(Literal('P') + integer + Optional(integer + integer)) + integer
cryst_expr.setParseAction(new_xtal_lattice)

origx_expr = Regex(r'ORIGX[123]') + \
    real("On1") + real("On2") + real("On3") + real("Tn")

scale_expr = Regex(r'SCALE[123]') + \
    real("Sn1") + real("Sn2") + real("Sn3") + real("Un")

master_expr = Record("MASTER", padwidth=4) + integer("numRemark") + \
    Suppress(Literal("0")) + integer("numHet") + integer("numHelix") + \
    integer("numSheet") + integer("numTurn") + integer("numSite") + \
    integer("numXform") + integer("numCoord") + integer("numTer") + \
    integer("numConect") + integer("numSeq")
obslte_expr = Keyword("OBSLTE") + restOfLine
split_expr = Keyword("SPLIT ") + restOfLine
caveat_expr = Keyword("CAVEAT") + restOfLine
nummdl_expr = Keyword("NUMMDL") + restOfLine
mdltyp_expr = Keyword("MDLTYP") + restOfLine
sprsde_expr = Keyword("SPRSDE") + restOfLine
dbref_expr = Record("DBREF ") + padexpr(id_code, padwidth=1) + \
    character_field(width=1, padwidth=1) + \
    integer_field(width=4, padwidth=0) + \
    character_field(width=1, padwidth=1) + \
    integer_field(width=4, padwidth=0) + \
    character_field(width=1, padwidth=1) + restOfLine
    # lstring_field(ci=27, cf=32, padwidth=1) + restOfLine

seqadv_expr = Record("SEQADV") + id_code + restOfLine
seqres_expr = Record("SEQRES") + restOfLine
modres_expr = Record("MODRES") + restOfLine
het_expr = Record("HET", width=6) + restOfLine
hetnam_expr = Record("HETNAM") + restOfLine
hetsyn_expr = Record("HETSYN") + restOfLine
formul_expr = Record("FORMUL") + restOfLine
helix_expr = Record("HELIX ") + restOfLine
sheet_expr = Record("SHEET ") + restOfLine
ssbond_expr = Record("SSBOND") + restOfLine
link_expr = Record("LINK  ") + restOfLine
cispep_expr = Record("CISPEP") + restOfLine
site_expr = Record("SITE  ") + restOfLine
mtrix_expr = Regex(r'MTRIX[123]') + integer("serial") + \
    real("Mn1") + real("Mn2") + real("Mn3") + real("Vn") + (integer | Empty())

model_expr = Record("MODEL ", width=6, padwidth=4) + integer("serial")
endmdl_expr = Record("ENDMDL")

anisou_expr = Record("ANISOU") + integer("serial") + \
    Word(alphanums, min=1, max=4)("name") + \
    (Word(alphanums, exact=1).addParseAction(matchOnlyAtCol(17))("altLoc") |
     nspaces(1)) + \
    Word(alphanums, min=1, max=3)("resName") + \
    (Word(alphanums, exact=1).addParseAction(matchOnlyAtCol(22))("chainID") |
     nspaces(1)) + \
    Word(nums, min=1, max=4).setParseAction(asint)("resSeq") + \
    (Word(alphas, exact=1).addParseAction(matchOnlyAtCol(27))("iCode") |
     nspaces(1)) + \
    integer("U11") + integer("U22") + integer("U33") + \
    integer("U12") + integer("U13") + integer("U23") + \
    Optional(oneOf(' '.join(element_symbols))("element")) + \
    Optional(Combine(integer + oneOf('+ -'))("charge"))

ter_expr = Record("TER", width=6) + integer("serial") + \
    Word(alphanums, min=1, max=3)("resName") + \
    Word(alphanums, exact=1).addParseAction(matchOnlyAtCol(22))("chainID") + \
    Word(nums, min=1, max=4).setParseAction(asint)("resSeq") + \
    Optional(Word(alphas, exact=1)
             .addParseAction(matchOnlyAtCol(27))("iCode") | nspaces(1))

conect_expr = Record("CONECT") + integer("serial") + OneOrMore(integer)
end_expr = Record("END", width=6)

atom_expr = Record("ATOM", width=6) + integer("serial") + \
    Word(alphanums, min=1, max=4)("name") + \
    (Word(alphanums, exact=1).addParseAction(matchOnlyAtCol(17))("altLoc") |
     nspaces(1)) + \
    Word(alphanums, min=1, max=3)("resName") + \
    (Word(alphanums, exact=1).addParseAction(matchOnlyAtCol(22))("chainID") |
     nspaces(1)) + \
    Word(nums, min=1, max=4).setParseAction(asint)("resSeq") + \
    (Word(alphas, exact=1).addParseAction(matchOnlyAtCol(27))("iCode") |
     nspaces(1)) + \
    real("x") + real("y") + real("z") + \
    real("occupancy") + real("tempFactor") + \
    Optional(oneOf(' '.join(element_symbols))("element")) + \
    Optional(Combine(integer + oneOf('+ -'))("charge"))
atom_expr.setParseAction(new_atom)

hetatm_expr = Record("HETATM") + integer("serial") + \
    Word(alphanums, min=1, max=4)("name") + \
    (Word(alphanums, exact=1).addParseAction(matchOnlyAtCol(17))("altLoc") |
     nspaces(1)) + \
    (Word(alphanums, min=1, max=3)("resName") | nspaces(3)) + \
    (Word(alphanums, exact=1).addParseAction(matchOnlyAtCol(22))("chainID") |
     nspaces(1)) + \
    (Word(nums, min=1, max=4).setParseAction(asint)("resSeq") |
     nspaces(4)) + \
    (Word(alphas, exact=1).addParseAction(matchOnlyAtCol(27))("iCode") |
     nspaces(1)) + \
    real("x") + real("y") + real("z") + \
    real("occupancy") + real("tempFactor") + \
    Optional(oneOf(' '.join(element_symbols))("element")) + \
    Optional(Combine(integer + oneOf('+ -'))("charge"))
hetatm_expr.setParseAction(new_atom)

mandatory_records = ['HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS',
                     'EXPDTA', 'AUTHOR', 'REVDAT', 'REMARK 2',
                     'REMARK 3', 'SEQRES', 'CRYST1', 'ORIGX1',
                     'ORIGX2', 'ORIGX3', 'SCALE1', 'SCALE2', 'SCALE3',
                     'MASTER', 'END']
optional_records = ['OBSLTE', 'SPLIT', 'CAVEAT', 'NUMMDL', 'MDLTYP',
                    'SPRSDE', 'JRNL', 'REMARK 0', 'REMARK 1',
                    'REMARK N', 'DBREF', 'DBREF1', 'DBREF2',
                    'SEQADV', 'MODRES', 'HET', 'HETNAM', 'HETSYN',
                    'FORMUL', 'HELIX', 'SHEET', 'SSBOND', 'LINK',
                    'CISPEP', 'SITE', 'MTRIX1', 'MTRIX2', 'MTRIX3',
                    'MODEL', 'ATOM', 'ANISOU', 'TER', 'HETATM',
                    'ENDMDL', 'CONECT']

record_expr = Forward()
record_expr << (header_expr | title_expr | compnd_expr | source_expr |
                keywds_expr | expdta_expr | author_expr | revdat_expr |
                remark_expr | seqres_expr | cryst_expr | origx_expr |
                scale_expr | master_expr | end_expr | obslte_expr |
                split_expr | caveat_expr | nummdl_expr | mdltyp_expr |
                sprsde_expr | jrnl_expr | dbref_expr | seqadv_expr |
                modres_expr | het_expr | hetnam_expr | hetsyn_expr |
                formul_expr | helix_expr | sheet_expr | ssbond_expr |
                link_expr | cispep_expr | site_expr | mtrix_expr |
                model_expr | atom_expr | anisou_expr | ter_expr |
                hetatm_expr | endmdl_expr | conect_expr)


class PDBTokenizer:
    """PDB Record tokenizer."""

    def __init__(self):
        self.parse_results = OrderedDict()
        self.raw_fields = OrderedDict()

    def _update_records(self, record, fields, parse_result):
        self.raw_fields[record] = fields
        self.parse_results[record] = parse_result

    def _parse_header(self, fields):
        # print('fields: {}'.format(fields))
        result = header_expr.parseString(fields)
        self._update_records('HEADER', fields, result)
        return result

        # classification = fields[10:50]
        # dep_date = fields[50:59]
        # id_code = fields[63:66]
        # print(classification)
        # print(dep_date)
        # print(id_code)

    def _parse_title(self, fields):
        result = title_expr.parseString(fields)
        self._update_records('TITLE', fields, result)
        return result

    def _parse_compnd(self, fields):
        pass

    def _parse_source(self, fields):
        pass

    def _parse_keywds(self, fields):
        pass

    def _parse_expdta(self, fields):
        pass

    def _parse_author(self, fields):
        pass

    def _parse_revdat(self, fields):
        pass

    def _parse_jrnl(self, fields):
        pass

    def _parse_remark(self, fields):
        pass

    def _parse_dbref(self, fields):
        pass

    def _parse_seqadv(self, fields):
        pass

    def _parse_seqres(self, fields):
        pass

    def _parse_formul(self, fields):
        pass

    def _parse_helix(self, fields):
        pass

    def _parse_sheet(self, fields):
        pass

    def _parse_ssbond(self, fields):
        pass

    def _parse_cryst1(self, fields):
        pass

    def _parse_origx1(self, fields):
        pass

    def _parse_origx2(self, fields):
        pass

    def _parse_origx3(self, fields):
        pass

    def _parse_scale1(self, fields):
        pass

    def _parse_scale2(self, fields):
        pass

    def _parse_scale3(self, fields):
        pass

    def _parse_atom(self, fields):
        result = atom_expr.parseString(fields, parseAll=True)[0]
        self._update_records('ATOM', fields, result)
        return result

    def _parse_ter(self, fields):
        pass

    def _parse_hetatm(self, fields):
        result = hetatm_expr.parseString(fields, parseAll=True)[0]
        self._update_records('HETATM', fields, result)
        return result

    def _parse_conect(self, fields):
        pass

    def _parse_master(self, fields):
        pass

    def _parse_end(self, fields):
        pass
