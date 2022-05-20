======================================================
ansys-mapdl-reader-support - Support library for pyMAPDL Reader
======================================================

This is a support library for the legacy pyMAPDL Reader library.  It
only contains Cython libraries used by the same.


Installation
------------
Installation through pip::

   pip install ansys-mapdl-reader-support

You can also visit `pymapdl-reader <https://github.com/pyansys/pymapdl-reader>`_
to download the source or releases from GitHub.

Developing on Windows
---------------------

This package is designed to be developed on Linux, and if you need to develop on Windows
you will need to install your own C++ compiler. We recommend:

 1. Install Visual C++
        a. See `here <https://wiki.python.org/moin/WindowsCompilers>`_ for a list of which Python versions correspond to which Visual C++ version
        b. Only Python <= 3.8 appears to be supported at the moment.
 2. Install the development version of pymapdl-reader-support to your Python environment
        a. Navigate to the pyMAPDL Reader project's top level (4 levels up from this file)
        b. run ``pip install -e ansys/mapdl/reader/cython``

To get the package up and running.
