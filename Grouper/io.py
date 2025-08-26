"""Module for working with external libraries."""

import importlib
import inspect
import os
import sys
import textwrap
from unittest import SkipTest


class DelayImportError(ImportError, SkipTest):
    """Error to allow better import handling."""

    pass


MESSAGES = dict()
MESSAGES["mbuild"] = """
The code at {filename}:{line_number} requires the "mbuild" package

mbuild can be installed with conda using:

# conda install -c conda-forge mbuild
"""

MESSAGES["torch"] = """
The code at {filename}:{line_number} requires the "torch" package

torch can be installed using:

# conda install torchvision -c pytorch
"""

MESSAGES["rdkit"] = """
The code at {filename}:{line_number} requires the "rdkit" package

rdkit can be installed with conda using:

# conda install -c conda-forge rdkit

or from source following instructions at:

https://www.rdkit.org/docs/Install.html#installation-from-source
"""


def import_(module):
    """Import a module and issue a nice message to stderr if it isn't installed.

    Parameters
    ----------
    module : str
        The module you'd like to import, as a string

    Returns
    -------
    module : {module, object}
        The module object

    Examples
    --------
    >>> # the following two lines are equivalent. the difference is that the
    >>> # second will check for an ImportError and print you a very nice
    >>> # user-facing message about what's wrong (where you can install the
    >>> # module from, etc) if the import fails
    >>> import tables
    >>> tables = import_('tables')

    """
    try:
        return importlib.import_module(module)
    except ImportError:
        try:
            message = MESSAGES[module]
        except KeyError:
            message = (
                "The code at {filename}:{line_number} requires the " f"{module} package"
            )

        (
            frame,
            filename,
            line_number,
            function_name,
            lines,
            index,
        ) = inspect.getouterframes(inspect.currentframe())[1]

        m = message.format(filename=os.path.basename(filename), line_number=line_number)
        m = textwrap.dedent(m)

        bar = (
            "\033[91m"
            + "#" * max(len(line) for line in m.split(os.linesep))
            + "\033[0m"
        )

        print("", file=sys.stderr)
        print(bar, file=sys.stderr)
        print(m, file=sys.stderr)
        print(bar, file=sys.stderr)
        raise DelayImportError(m)


try:
    import mbuild

    has_mbuild = True
    del mbuild
except ImportError:
    has_mbuild = False

try:
    import torch

    has_torch = True
    del torch
except ImportError:
    has_torch = False

try:
    import rdkit

    has_rdkit = True
    del rdkit
except ImportError:
    has_rdkit = False
