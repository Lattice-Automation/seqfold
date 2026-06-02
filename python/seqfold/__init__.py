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

from ._core import (
    Struct,
    dg,
    dg_cache,
    dot_bracket,
    fold,
    gc_cache,
    tm,
    tm_cache,
)

__all__ = [
    "Struct",
    "dg",
    "dg_cache",
    "dot_bracket",
    "fold",
    "gc_cache",
    "tm",
    "tm_cache",
]
