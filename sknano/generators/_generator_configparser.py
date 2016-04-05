from __future__ import absolute_import, division, print_function, \
    unicode_literals

from collections import OrderedDict
from configparser import ConfigParser
import importlib

from sknano.core import BaseClass, call_signature

__all__ = ['GeneratorConfigParser']


class GeneratorConfigParser(BaseClass):
    """Class for reading/writing ini config files."""
    call_signature = call_signature.copy()

    def __init__(self, cfgfile=None, structure=None, **kwargs):
        super().__init__(**kwargs)
        self.cfgfile = cfgfile
        self.structure = structure

        self.fmtstr = "{cfgfile!r}"

        self.parser = ConfigParser()
        self.config = OrderedDict()
        if cfgfile is not None:
            self._parse_config()

        if structure is not None:
            self._update_config()

    def _parse_config(self):
        parser = self.parser
        parser.read(self.cfgfile)
        generator_module = 'sknano.generators'

        for generator_class in parser.sections():
            # generator_module = parser[section]['generator_module']
            parameters = parser[generator_class]['parameters']
            # fname = '{}({})'.format(generator_class[:-len('Generator')],
            #                         parameters)
            # self.fnames.append(fname.replace(' ', ''))
            # self.config.update({section: {'generator_class': generator_class,
            #                               'parameters': parameters}})

            call_sig = \
                self.call_signature.parseString(parameters, parseAll=True)[0]
            try:
                args, kwargs = call_sig
            except ValueError:
                args, kwargs = tuple(), call_sig[0]

            generator = getattr(importlib.import_module(generator_module),
                                generator_class)(*args, **kwargs)
            print(generator)

    def _update_config(self):
        if self.structure is not None:
            parameter_dict = self.structure.todict()
            for param, value in parameter_dict.items():
                print('{} = {}'.format(param, value))

    # def write(self, fp, space_around_delimiters=True):
    #     """Write an .ini-format representation of the configuration state.

    #     If `space_around_delimiters' is True (the default), delimiters
    #     between keys and values are surrounded by spaces.
    #     """
    #     if self.structure is not None:
    #     for k, v in
    #     # if space_around_delimiters:
    #     #     d = " {} ".format(self._delimiters[0])
        # else:
        #     d = self._delimiters[0]
        # if self._defaults:
        #     self._write_section(fp, self.default_section,
        #                         self._defaults.items(), d)
        # for section in self._sections:
        #     self._write_section(fp, section,
        #                         self._sections[section].items(), d)

    # def _write_section(self, fp, section_name, section_items, delimiter):
    #     """Write a single section to the specified `fp'."""
    #     fp.write("[{}]\n".format(section_name))
    #     for key, value in section_items:
    #         value = self._interpolation.before_write(self, section_name, key,
    #                                                  value)
    #         if value is not None or not self._allow_no_value:
    #             value = delimiter + str(value).replace('\n', '\n\t')
    #         else:
    #             value = ""
    #         fp.write("{}{}\n".format(key, value))
    #     fp.write("\n")

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(cfgfile=self.cfgfile, structure=self.structure)
