======================================================
PyMAPDL Reader - Legacy Binary and Archive File Reader
======================================================
This is the legacy module for reading in binary and ASCII files
generated from MAPDL.

This Python module allows you to extract data directly from binary
ANSYS v14.5+ files and to display or animate them rapidly using a
straightforward API coupled with C libraries based on header files
provided by ANSYS.

The ``ansys-mapdl-reader`` module supports the following formats:
  - ``*.rst`` - Structural analysis result file
  - ``*.rth`` - Thermal analysis result file 
  - ``*.emat`` - Element matrice data file
  - ``*.full`` - Full stiffness-mass matrix file
  - ``*.cdb`` or ``*.dat`` - MAPDL ASCII block archive and
    Mechanical Workbench input files

Please see the :ref:`ref_example_gallery` for several demos using
``ansys-mapdl-reader``.

.. note::

   This module will likely change or be depreciated in the future.

   You are encouraged to use the new Data Processing Framework (DPF)
   modules at `DPF-Core <https://github.com/pyansys/DPF-Core>`_ and
   `DPF-Post <https://github.com/pyansys/DPF-Post>`_ as they provide a
   modern interface to ANSYS result files using a client/server
   interface using the same software used within ANSYS Workbench, but
   via a Python client.


Brief Demo: Direct Access to Binary Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Here's a quick example code block to show how easy it is to load and
plots results directly from an ANSYS result file using
``ansys-mapdl-reader``:

.. code:: python

    >>> from ansys.mapdl import reader as pymapdl_reader
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
    >>> estress[0]
    array([[ 1.0236604e+04, -9.2875127e+03, -4.0922625e+04, -2.3697146e+03,
            -1.9239732e+04,  3.0364934e+03]
           [ 5.9612605e+04,  2.6905924e+01, -3.6161423e+03,  6.6281304e+03,
             3.1407712e+02,  2.3195926e+04]
	   [ 3.8178301e+04,  1.7534495e+03, -2.5156013e+02, -6.4841372e+03,
            -5.0892783e+03,  5.2503605e+00]
	   [ 4.9787645e+04,  8.7987168e+03, -2.1928742e+04, -7.3025332e+03,
             1.1294199e+04,  4.3000205e+03]])

    Get the corresponding ansys element number for the first element.

    >>> elem[0]
        32423

    Get the nodes belonging to for the first element.

    >>> enode[0]
        array([ 9012,  7614,  9009, 10920], dtype=int32)

    Plot the nodal displacement of the first result.

    >>> result.plot_nodal_displacement(0)

.. figure:: ./images/rotor.jpg
    :width: 500pt


.. toctree::
   :maxdepth: 1
   :hidden:

   getting_started
   user_guide/index
   api/index
   examples/index
   contributing


License and Acknowledgments
---------------------------
The ``ansys-mapdl-reader`` module is licensed under the MIT license.
