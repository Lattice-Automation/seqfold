# -*- coding: utf-8 -*-
from importlib.metadata import version, PackageNotFoundError

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = version(dist_name)
except PackageNotFoundError:
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError

from .fold import fold, dg, dg_cache, fold, Struct, dot_bracket
from .tm import tm, tm_cache, gc_cache
from .types import Cache
