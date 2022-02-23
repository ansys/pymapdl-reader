Installation
------------
Installation is simply::

    pip install ansys-mapdl-reader

You can also visit `pymapdl-reader <https://github.com/pyansys/pymapdl-reader>`_
to download the source or releases from GitHub.

If you have any installation (or other) issues, please open an issue
at `pymapdl-reader Issues <https://github.com/pyansys/pymapdl-reader/issues>`_.

The ``ansys-mapdl-reader`` supports 64-Bit Windows, Mac OS, and Linux
for Python 3.7 through Pythou 3.9.


Python 3.10 Extra Instructions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PyMAPDL-Reader requires the `VTK library <https://vtk.org/>`_ which, at the
moment, is not available for Python 3.10 in `their official channel
<https://pypi.org/project/vtk/>`_.

If you wish to install PyMAPDL-Reader in Python 3.10, you can still do it by
using the unofficial VTK wheel from PyVista using ``--find-links``. This tells ``pip`` to look for vtk at `wheels.pyvista.org <https://wheels.pyvista.org/>`_. Use this with::

    pip install ansys-mapdl-reader --find-links https://wheels.pyvista.org/

Please visit `Unofficial VTK Wheels for Python 3.10
<https://github.com/pyvista/pyvista/discussions/2064>`_ for further details.
