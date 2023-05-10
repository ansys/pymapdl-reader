======================================================
PyMAPDL Reader - Legacy Binary and Archive File Reader
======================================================
|pyansys| |pypi| |PyPIact| |GH-CI| |codecov| |MIT| |black| |pre-commit|

.. |pyansys| image:: https://img.shields.io/badge/Py-Ansys-ffc107.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAIAAACQkWg2AAABDklEQVQ4jWNgoDfg5mD8vE7q/3bpVyskbW0sMRUwofHD7Dh5OBkZGBgW7/3W2tZpa2tLQEOyOzeEsfumlK2tbVpaGj4N6jIs1lpsDAwMJ278sveMY2BgCA0NFRISwqkhyQ1q/Nyd3zg4OBgYGNjZ2ePi4rB5loGBhZnhxTLJ/9ulv26Q4uVk1NXV/f///////69du4Zdg78lx//t0v+3S88rFISInD59GqIH2esIJ8G9O2/XVwhjzpw5EAam1xkkBJn/bJX+v1365hxxuCAfH9+3b9/+////48cPuNehNsS7cDEzMTAwMMzb+Q2u4dOnT2vWrMHu9ZtzxP9vl/69RVpCkBlZ3N7enoDXBwEAAA+YYitOilMVAAAAAElFTkSuQmCC
   :target: https://docs.pyansys.com/
   :alt: PyAnsys

.. |pypi| image:: https://img.shields.io/pypi/v/ansys-mapdl-reader.svg?logo=python&logoColor=white
   :target: https://pypi.org/project/ansys-mapdl-reader/

.. |PyPIact| image:: https://img.shields.io/pypi/dm/ansys-mapdl-reader.svg?label=PyPI%20downloads
   :target: https://pypi.org/project/ansys-mapdl-reader/

.. |codecov| image:: https://codecov.io/gh/pyansys/pymapdl-reader/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/pyansys/pymapdl-reader

.. |GH-CI| image:: https://github.com/pyansys/pymapdl-reader/actions/workflows/testing-and-deployment.yml/badge.svg
   :target: https://github.com/pyansys/pymapdl-reader/actions/workflows/testing-and-deployment.yml

.. |MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg?style=flat
  :target: https://github.com/psf/black
  :alt: black

.. |pre-commit| image:: https://results.pre-commit.ci/badge/github/pyansys/pymapdl-reader/main.svg
   :target: https://results.pre-commit.ci/latest/github/pyansys/pymapdl-reader/main
   :alt: pre-commit.ci status

This is the legacy module for reading in binary and ASCII files
generated from MAPDL.

This Python module allows you to extract data directly from binary
ANSYS v14.5+ files and to display or animate them rapidly using a
straightforward API coupled with C libraries based on header files
provided by ANSYS.

The ``ansys-mapdl-reader`` module supports the following formats:

- ``*.rst`` - Structural analysis result file
- ``*.rth`` - Thermal analysis result file 
- ``*.emat`` - Element matrix data file
- ``*.full`` - Full stiffness-mass matrix file
- ``*.cdb`` or ``*.dat`` - MAPDL ASCII block archive and
  Mechanical Workbench input files

Please see the `PyMAPDL-Reader Documentation
<https://readerdocs.pyansys.com>`_ for the full documentation.

.. note::

   This module may be depreciated in the future.

   You are encouraged to use the new Data Processing Framework (DPF)
   modules at `PyDPF-Core <https://github.com/pyansys/pydpf-core>`_ and
   `PyDPF-Post <https://github.com/pyansys/pydpf-post>`_ as they provide a
   modern interface to Ansys result files using a client/server
   interface using the same software used within Ansys Mechanical, but
   via a Python client.


Installation
------------
Installation through pip::

   pip install ansys-mapdl-reader

You can also visit `pymapdl-reader <https://github.com/pyansys/pymapdl-reader>`_
to download the source or releases from GitHub.


Examples
--------

Loading and Plotting a MAPDL Archive File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ANSYS archive files containing solid elements (both legacy and
modern), can be loaded using Archive and then converted to a vtk
object.

.. code:: python

    from ansys.mapdl import reader as pymapdl_reader
    from ansys.mapdl.reader import examples
    
    # Sample *.cdb
    filename = examples.hexarchivefile
    
    # Read ansys archive file
    archive = pymapdl_reader.Archive(filename)
    
    # Print raw data from cdb
    for key in archive.raw:
       print("%s : %s" % (key, archive.raw[key]))
    
    # Create a vtk unstructured grid from the raw data and plot it
    grid = archive.parse_vtk(force_linear=True)
    grid.plot(color='w', show_edges=True)
    
    # write this as a vtk xml file 
    grid.save('hex.vtu')

    # or as a vtk binary
    grid.save('hex.vtk')


.. figure:: https://github.com/pyansys/pymapdl-reader/blob/main/doc/source/images/hexbeam_small.png
   :alt: Hexahedral beam

You can then load this vtk file using ``pyvista`` or another program that uses VTK.
    
.. code:: python

    # Load this from vtk
    import pyvista as pv
    grid = pv.UnstructuredGrid('hex.vtu')
    grid.plot()


Loading the Result File
~~~~~~~~~~~~~~~~~~~~~~~
This example reads in binary results from a modal analysis of a beam
from ANSYS.

.. code:: python

    # Load the reader from pyansys
    from ansys.mapdl import reader as pymapdl_reader
    from ansys.mapdl.reader import examples
    
    # Sample result file
    rstfile = examples.rstfile
    
    # Create result object by loading the result file
    result = pymapdl_reader.read_binary(rstfile)
    
    # Beam natural frequencies
    freqs = result.time_values

.. code:: python

    >>> print(freq)
    [ 7366.49503969  7366.49503969 11504.89523664 17285.70459456
      17285.70459457 20137.19299035]
    
Get the 1st bending mode shape.  Results are ordered based on the
sorted node numbering.  Note that results are zero indexed

.. code:: python

    >>> nnum, disp = result.nodal_solution(0)
    >>> print(disp)
    [[ 2.89623914e+01 -2.82480489e+01 -3.09226692e-01]
     [ 2.89489249e+01 -2.82342416e+01  2.47536161e+01]
     [ 2.89177130e+01 -2.82745126e+01  6.05151053e+00]
     [ 2.88715048e+01 -2.82764960e+01  1.22913304e+01]
     [ 2.89221536e+01 -2.82479511e+01  1.84965333e+01]
     [ 2.89623914e+01 -2.82480489e+01  3.09226692e-01]
     ...


Plotting Nodal Results
~~~~~~~~~~~~~~~~~~~~~~
As the geometry of the model is contained within the result file, you
can plot the result without having to load any additional geometry.
Below, displacement for the first mode of the modal analysis beam is
plotted using ``VTK``.

.. code:: python
    
    # Plot the displacement of Mode 0 in the x direction
    result.plot_nodal_solution(0, 'x', label='Displacement')

.. figure:: https://github.com/pyansys/pymapdl-reader/blob/main/doc/source/images/hexbeam_disp_small.png


Results can be plotted non-interactively and screenshots saved by
setting up the camera and saving the result.  This can help with the
visualization and post-processing of a batch result.

First, get the camera position from an interactive plot:

.. code:: python

    >>> cpos = result.plot_nodal_solution(0)
    >>> print(cpos)
    [(5.2722879880979345, 4.308737919176047, 10.467694436036483),
     (0.5, 0.5, 2.5),
     (-0.2565529433509593, 0.9227952809887077, -0.28745339908049733)]

Then generate the plot:

.. code:: python

    result.plot_nodal_solution(0, 'x', label='Displacement', cpos=cpos,
                               screenshot='hexbeam_disp.png',
                               window_size=[800, 600], interactive=False)

Stress can be plotted as well using the below code.  The nodal stress
is computed in the same manner that ANSYS uses by to determine the
stress at each node by averaging the stress evaluated at that node for
all attached elements.  For now, only component stresses can be
displayed.

.. code:: python
    
    # Display node averaged stress in x direction for result 6
    result.plot_nodal_stress(5, 'Sx')

.. figure:: https://github.com/pyansys/pymapdl-reader/blob/main/doc/source/images/beam_stress_small.png


Nodal stress can also be generated non-interactively with:

.. code:: python

    result.plot_nodal_stress(5, 'Sx', cpos=cpos, screenshot=beam_stress.png,
                           window_size=[800, 600], interactive=False)


Animating a Modal Solution
~~~~~~~~~~~~~~~~~~~~~~~~~~
Mode shapes from a modal analysis can be animated using ``animate_nodal_solution``:

.. code:: python

    result.animate_nodal_solution(0)


.. figure:: https://github.com/pyansys/pymapdl-reader/blob/main/doc/source/images/beam_mode_shape_small.gif
   :alt: Modal shape animation

If you wish to save the animation to a file, specify the
movie_filename and animate it with:

.. code:: python

    result.animate_nodal_solution(0, movie_filename='/tmp/movie.mp4', cpos=cpos)


Reading a Full File
~~~~~~~~~~~~~~~~~~~
This example reads in the mass and stiffness matrices associated with
the above example.

.. code:: python

    # Load the reader from pyansys
    from ansys.mapdl import reader as pymapdl_reader
    from scipy import sparse
    
    # load the full file
    fobj = pymapdl_reader.FullReader('file.full')
    dofref, k, m = fobj.load_km()  # returns upper triangle only

    # make k, m full, symmetric matrices
    k += sparse.triu(k, 1).T
    m += sparse.triu(m, 1).T

If you have ``scipy`` installed, you can solve the eigensystem for its
natural frequencies and mode shapes.

.. code:: python

    from scipy.sparse import linalg

    # condition the k matrix
    # to avoid getting the "Factor is exactly singular" error
    k += sparse.diags(np.random.random(k.shape[0])/1E20, shape=k.shape)

    # Solve
    w, v = linalg.eigsh(k, k=20, M=m, sigma=10000)

    # System natural frequencies
    f = np.real(w)**0.5/(2*np.pi)
    
    print('First four natural frequencies')
    for i in range(4):
        print '{:.3f} Hz'.format(f[i])
    
.. code::

    First four natural frequencies
    1283.200 Hz
    1283.200 Hz
    5781.975 Hz
    6919.399 Hz

Developing on Windows
---------------------

This package is designed to be developed on Linux, and if you need to develop on Windows
you will need to install your own C++ compiler. We recommend:

1. Install Visual C++
       a. See `here <https://wiki.python.org/moin/WindowsCompilers>`_ for a list of which Python versions correspond to which Visual C++ version
2. Install the development version of pymapdl-reader to your Python environment
       a. Navigate to the project's top level (the same directory as this README)
       b. run ``pip install -e .``


License and Acknowledgments
---------------------------
The ``ansys-mapdl-reader`` library is licensed under the MIT license.
