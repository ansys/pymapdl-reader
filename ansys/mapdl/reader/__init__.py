import os

import appdirs

from ansys.mapdl.reader._version import __version__
from ansys.mapdl.reader.misc import _configure_pyvista

from . import _archive

# Setup data directory
USER_DATA_PATH = appdirs.user_data_dir(appname="ansys_mapdl_reader", appauthor="Ansys")

EXAMPLES_PATH = os.path.join(USER_DATA_PATH, "examples")

# set pyvista defaults
_configure_pyvista()
