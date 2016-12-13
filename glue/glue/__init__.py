"""
The glue package has been renamed to ligo.glue. Provide import hook so that
importing `glue.abc.def` maps to `ligo.glue.abc.def` and prints a deprecation
warning.
"""
from __future__ import absolute_import
from importlib import import_module as _import_module
from warnings import warn as _warn
import sys as _sys


def _complain(name):
    _msg = '{0}: deprecated module, please import ligo.{0} instead.'
    _warn(_msg.format(name), DeprecationWarning, 3)


class _DeprecatedGlueImportHook(object):

    def find_module(self, fullname, path=None):
        toplevel, _, _ = fullname.partition('.')
        if toplevel == 'glue':
            return self
        else:
            return None

    def load_module(self, fullname):
        _complain(fullname)
        realname = 'ligo.' + fullname
        try:
            module = _sys.modules[realname]
        except KeyError:
            module = _sys.modules[realname] = _import_module(realname)
        return module


_complain('glue')
_sys.meta_path.append(_DeprecatedGlueImportHook())
from ligo.glue import *
