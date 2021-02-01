import appdirs
import os

# Per contract with Sphinx-Gallery, this method must be available at top level
from pyvista.utilities.sphinx_gallery import _get_sg_image_scraper

from ansys.mapdl.reader._version import __version__
from ansys.mapdl.reader.archive import (Archive, write_cmblock, write_nblock,
                                        save_as_archive)
from ansys.mapdl.reader.cell_quality import quality
from ansys.mapdl.reader.common import read_binary
from ansys.mapdl.reader.misc import Report, _configure_pyvista
from ansys.mapdl.reader import examples
from . import _archive

# Setup data directory
try:
    USER_DATA_PATH = appdirs.user_data_dir('ansys_mapdl_reader')
    if not os.path.exists(USER_DATA_PATH):  # pragma: no cover
        os.makedirs(USER_DATA_PATH)

    EXAMPLES_PATH = os.path.join(USER_DATA_PATH, 'examples')
    if not os.path.exists(EXAMPLES_PATH):  # pragma: no cover
        os.makedirs(EXAMPLES_PATH)

except:  # pragma: no cover
    pass

# set pyvista defaults
_configure_pyvista()
