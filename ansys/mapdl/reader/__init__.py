import os

import appdirs

# Per contract with Sphinx-Gallery, this method must be available at top level
from ansys.mapdl.reader.misc import _configure_pyvista
try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata

# Setup data directory
USER_DATA_PATH = appdirs.user_data_dir(appname="ansys_mapdl_reader", appauthor="Ansys")

EXAMPLES_PATH = os.path.join(USER_DATA_PATH, "examples")

# set pyvista defaults
_configure_pyvista()

__version__ = importlib_metadata.version(__name__.replace(".", "-"))
"""ansys-mapdl-reader version."""
