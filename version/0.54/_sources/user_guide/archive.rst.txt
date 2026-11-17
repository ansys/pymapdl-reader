Reading and Writing Mapdl Archive Files
=======================================

Reading ANSYS Archives
----------------------
MAPDL archive ``*.cdb`` and ``*.dat`` files containing elements (both
legacy and modern) can be loaded using Archive and then converted to a
``vtk`` object:

.. code:: python

    from ansys.mapdl import reader as pymapdl_reader
    from ansys.mapdl.reader import examples

    # Read a sample archive file
    archive = pymapdl_reader.Archive(examples.hexarchivefile)

    # Print various raw data from cdb
    print(archive.nnum, archive.nodes)

    # access a vtk unstructured grid from the raw data and plot it
    grid = archive.grid
    archive.plot(color='w', show_edges=True)


You can also optionally read in any stored parameters within the
archive file by enabling the ``read_parameters`` parameter.

.. code:: python

    from ansys.mapdl import reader as pymapdl_reader
    archive = pymapdl_reader.Archive('mesh.cdb', read_parameters=True)

    # parameters are stored as a dictionary
    archive.parameters

See the `Archive` class documentation below for more details on the
class methods and properties.


Writing ANSYS Archives
----------------------
Unstructured grids generated using VTK can be converted to ANSYS APDL
archive files and loaded into any version of ANSYS using
``pymapdl_reader.save_as_archive``.  The following example using the
built-in archive file demonstrates this capability.

.. code:: python

    import pyvista as pv
    from pyvista import examples
    from ansys.mapdl import reader as pymapdl_reader

    # load in a vtk unstructured grid
    grid = pv.UnstructuredGrid(examples.hexbeamfile)
    script_filename = '/tmp/grid.cdb'
    pymapdl_reader.save_as_archive(script_filename, grid)

    # optionally read in archive in ANSYS and generate cell shape
    # quality report
    from ansys.mapdl.core import launch_mapdl
    mapdl = launch_mapdl()
    mapdl.cdread('db', script_filename)
    mapdl.prep7()
    mapdl.shpp('SUMM')

Resulting ANSYS quality report:

.. code::

    ------------------------------------------------------------------------------
               <<<<<<          SHAPE TESTING SUMMARY           >>>>>>
               <<<<<<        FOR ALL SELECTED ELEMENTS         >>>>>>
    ------------------------------------------------------------------------------
                       --------------------------------------
                       |  Element count        40 SOLID185  |
                       --------------------------------------
   
     Test                Number tested  Warning count  Error count    Warn+Err %
     ----                -------------  -------------  -----------    ----------
     Aspect Ratio                 40              0             0         0.00 %
     Parallel Deviation           40              0             0         0.00 %
     Maximum Angle                40              0             0         0.00 %
     Jacobian Ratio               40              0             0         0.00 %
     Warping Factor               40              0             0         0.00 %
   
     Any                          40              0             0         0.00 %
    ------------------------------------------------------------------------------


Converting a MAPDL Archive File to VTK for Paraview
---------------------------------------------------
MAPDL archive files containing solid elements (both legacy and modern)
can be loaded using Archive and then converted to a VTK object.

.. code:: python

    from ansys.mapdl import reader as pymapdl_reader
    from ansys.mapdl.reader import examples

    # Sample *.cdb
    filename = examples.hexarchivefile

    # Read ansys archive file
    archive = pymapdl_reader.Archive(filename)

    # Print overview of data read from cdb
    print(archive)

    # Create a vtk unstructured grid from the raw data and plot it
    grid = archive.parse_vtk(force_linear=True)
    grid.plot(color='w', show_edges=True)

    # save this as a vtk xml file 
    grid.save('hex.vtu')

.. image:: ../images/hexbeam.png


You can then load this vtk file using ``pyvista`` or another program that uses VTK.
    
.. code:: python

    # Load this from vtk
    import pyvista as pv
    grid = pv.read('hex.vtk')
    grid.plot()


Supported Elements
~~~~~~~~~~~~~~~~~~
At the moment, only solid elements are supported by the
``save_as_archive`` function, to include:

 - ``vtk.VTK_TETRA``
 - ``vtk.VTK_QUADRATIC_TETRA``
 - ``vtk.VTK_PYRAMID``
 - ``vtk.VTK_QUADRATIC_PYRAMID``
 - ``vtk.VTK_WEDGE``
 - ``vtk.VTK_QUADRATIC_WEDGE``
 - ``vtk.VTK_HEXAHEDRON``
 - ``vtk.VTK_QUADRATIC_HEXAHEDRON``

Linear element types will be written as SOLID185, quadratic elements
will be written as SOLID186, except for quadratic tetrahedrals, which
will be written as SOLID187.


Archive Class
-------------
.. autoclass:: ansys.mapdl.reader.archive.Archive
    :members:
    :inherited-members:
