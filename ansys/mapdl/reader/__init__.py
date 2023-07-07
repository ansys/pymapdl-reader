import os

import appdirs

# Per contract with Sphinx-Gallery, this method must be available at top level
try:
    # for pyvista >= 0.40
    from pyvista.plotting.utilities import _get_sg_image_scraper
except ImportError:
    from pyvista.utilities.sphinx_gallery import _get_sg_image_scraper


from ansys.mapdl.reader import examples
from ansys.mapdl.reader._version import __version__
from ansys.mapdl.reader.archive import (
    Archive,
    save_as_archive,
    write_cmblock,
    write_nblock,
)
from ansys.mapdl.reader.cell_quality import quality
from ansys.mapdl.reader.common import read_binary
from ansys.mapdl.reader.misc import Report, _configure_pyvista

from . import _archive

# Setup data directory
USER_DATA_PATH = appdirs.user_data_dir(appname="ansys_mapdl_reader", appauthor="Ansys")

EXAMPLES_PATH = os.path.join(USER_DATA_PATH, "examples")

# set pyvista defaults
_configure_pyvista()
