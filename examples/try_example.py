import pyvista
import numpy as np

from ansys.mapdl import reader as pymapdl_reader
from ansys.mapdl.reader import examples

# Download an example shaft modal analysis result file
shaft = examples.download_shaft_modal()
