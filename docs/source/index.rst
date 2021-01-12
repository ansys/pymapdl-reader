PyMAPDL Legacy Binary and Archive Reader
========================================
This is the legacy module for reading in binary and ASCII files
generated from MAPDL.

This Python module allows you to:
 - Extract data directly from binary ANSYS v14.5+ files and to display
   or animate them.  Specifically, the following formats:

    - ``*.rst``: Result file from structural analysis
    - ``*.rth``: Result file from a thermal analysis
    - ``*.emat``: Stores data related to element matrices
    - ``*.full``: Stores the full stiffness-mass matrix
    - Reading nodes and elements from MAPDL ASCII block archive
      ``*.cdb`` and ``*.dat`` files

This module will be maintained provided that it provides unique
support for reading files from MAPDL and may be subject to
depreciation when ANSYS provides better support for a variety of file
formats.

In addition to this module, you are encouraged to checkout the new
Data Processing Framework (DPF) modules at `DPF-Core
<https://github.com/pyansys/DPF-Core>`_ and `DPF-Post
<https://github.com/pyansys/DPF-Post>`_ as they provide a modern
interface to ANSYS result files using a client/server interface using
the same software used within ANSYS Workbench, but via a Python
client.

Please see the :ref:`ref_example_gallery` for several demos using
``ansys-mapdl-reader``.

Installation
------------
Installation through pip::

    pip install ansys-mapdl-reader

You can also visit `pymapdl-reader <https://github.com/pyansys/pymapdl-reader>`_
to download the source or releases from GitHub.

If you have any installation (or other) issues, please open an issue
at `pymapdl-reader Issues <https://github.com/pyansys/pymapdl-reader/issues>`_.


Quick Examples
--------------

Direct Access to Binary Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Here's a quick example code block to show how easy it is to load and
plots results directly from an ANSYS result file using ``pymapdl-reader``:

.. code:: python

    >>> import ansys.mapdl.reader as pymapdl_reader
    >>> result = pymapdl_reader.read_binary('rotor.rst')
    >>> nnum, stress = result.principal_nodal_stress(0)
    >>> print(stress[:3])
    array([[-1.2874937e-06,  1.2874934e-06,  5.6843419e-14,  0.0000000e+00,
             8.1756007e-06, -8.1756007e-06],
           [-1.1674185e-04, -1.1674478e-04, -3.0856981e-04, -1.7892545e-06,
            -2.5823609e-05,  2.5835518e-05],
           [-5.7354209e-05, -5.5398770e-05, -1.4944717e-04, -1.0580692e-06,
            -1.7659733e-05, -3.5462126e-06]], dtype=float32)

    Element stress for the first result

    >>> estress, elem, enode = result.element_stress(0)
    >>> estress[0]  # element stress for element 0
    array([[ 1.0236604e+04, -9.2875127e+03, -4.0922625e+04, -2.3697146e+03,
            -1.9239732e+04,  3.0364934e+03]
           [ 5.9612605e+04,  2.6905924e+01, -3.6161423e+03,  6.6281304e+03,
             3.1407712e+02,  2.3195926e+04]
	   [ 3.8178301e+04,  1.7534495e+03, -2.5156013e+02, -6.4841372e+03,
            -5.0892783e+03,  5.2503605e+00]
	   [ 4.9787645e+04,  8.7987168e+03, -2.1928742e+04, -7.3025332e+03,
             1.1294199e+04,  4.3000205e+03]])

    >>> elem[0]  # corresponding ansys element number for element 0
        32423

    >>> enode[0]  # corresponding nodes belonging to for element 0
        array([ 9012,  7614,  9009, 10920], dtype=int32)

    >>> result.plot_nodal_solution(0)

.. figure:: ./images/rotor.jpg
    :width: 500pt

--------

Contents
========

.. toctree::
   :maxdepth: 1
   :caption: ANSYS File Support

   archive
   loading_results
   examples
   loading_km
   loading_emat

.. toctree::
   :maxdepth: 1
   :caption: Example Gallery

   examples/index

.. toctree::
   :maxdepth: 1
   :caption: Miscellaneous
   :hidden:

   quality

--------

License and Acknowledgments
---------------------------
The ``ansys-mapdl-reader`` module is licensed under the MIT license.
