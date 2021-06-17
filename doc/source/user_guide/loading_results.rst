Reading MAPDL Result Files
==========================
The `ansys-mapdl-reader` module supports the following result types from MAPDL:

- ``".rfl"``
- ``".rmg"``
- ``".rst"``
- ``".rth"``

The MAPDL result file is a FORTRAN formatted binary file containing
the results written from a MAPDL analysis.  The results, at a
minimum, contain the geometry of the model analyzed along with the
nodal and element results.  Depending on the analysis, these results
could be anything from modal displacements to nodal temperatures.
This includes (and is not limited to):

    - Nodal DOF results from a static analysis or modal analysis.
    - Nodal DOF results from a cyclic static or modal analysis.
    - Nodal averaged component stresses (i.e. x, y, z, xy, xz, yz)
    - Nodal principal stresses (i.e. S1, S2, S3, SEQV, SINT)
    - Nodal elastic, plastic, and thermal stress
    - Nodal time history
    - Nodal boundary conditions and force
    - Nodal temperatures
    - Nodal thermal strain
    - Various element results (see ``element_solution_data``)

This module will likely change or depreciated in the future, and you
are encouraged to checkout the new Data Processing Framework (DPF)
modules at `DPF-Core <https://github.com/pyansys/DPF-Core>`_ and
`DPF-Post <https://github.com/pyansys/DPF-Post>`_ as they provide a
modern interface to ANSYS result files using a client/server interface
using the same software used within ANSYS Workbench, but via a Python
client.


Loading the Result File
-----------------------
As the MAPDL result files are binary files, the entire file does not
need to be loaded into memory in order to retrieve results.  This
module accesses the results through a python object `result` which you
can initialize with:

.. code:: python

    from ansys.mapdl import reader as pymapdl_reader
    result = pymapdl_reader.read_binary('file.rst')
    
Upon initialization the ``Result`` object contains several
properties to include the time values from the analysis, node
numbering, element numbering, etc.

The ``ansys-mapdl-reader`` module can determine the correct result
type by reading the header of the file, which means that if it is an
MAPDL binary file, ``ansys-mapdl-reader`` can probably read it (at
least to some degree.  For example, a thermal result file can be read
with

.. code:: python

    rth = pymapdl_reader.read_binary('file.rth')


Result Properties
-----------------
The properties of the ``Result`` can be quickly shown by printing the
result file with:

.. code:: python

    >>> result = pymapdl_reader.read_binary('file.rst')
    >>> print(result)
    PyMAPDL Result file object
    Units       : User Defined
    Version     : 20.1
    Cyclic      : False
    Result Sets : 1
    Nodes       : 321
    Elements    : 40


    Available Results:
    EMS : Miscellaneous summable items (normally includes face pressures)
    ENF : Nodal forces
    ENS : Nodal stresses
    ENG : Element energies and volume
    EEL : Nodal elastic strains
    ETH : Nodal thermal strains (includes swelling strains)
    EUL : Element euler angles
    EPT : Nodal temperatures
    NSL : Nodal displacements
    RF  : Nodal reaction forces


To obtain the time or frequency values of an analysis use:
    
.. code:: python

    >>> result.time_values
    array([1.])


Individual results can be obtained with one of the many methods
available to the result object.  For example, the nodal displacement
for the first result can be accessed with:

.. code:: python

    >>> nnum, disp = rst.nodal_displacement(0)
    >>> nnum
    array([  1,   2,   3, ..., 318, 319, 320, 321], dtype=int32)

    >>> disp
    array([[-2.03146520e-09, -3.92491045e-03,  5.00047448e-05],
           [ 1.44630651e-09,  1.17747356e-02, -1.49992672e-04],
           [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
           ...
           [-7.14982194e-03,  3.12495002e-03,  5.74992265e-04],
           [-7.04982329e-03,  2.44996706e-03,  5.74992939e-04],
           [-6.94982520e-03,  1.77498362e-03,  5.74992891e-04]])


The sorted node and element numbering of a result can be obtained with:

.. code:: python

    >>> rst.geometry.nnum
    array([  1,   2,   3, ..., 318, 319, 320, 321], dtype=int32)

    >>> result.geometry.enum
    array([ 1,  3,  2,  4,  5,  7,  6,  8,  9, 11, 10, 12, 13, 15, 14, 16, 17,
           19, 18, 20, 21, 23, 22, 24, 25, 27, 26, 28, 29, 31, 30, 32, 33, 35,
           34, 36, 37, 39, 38, 40], dtype=int32)

Mesh
----
The mesh of the result can be found by querying the ``mesh`` property
of a result, which returns a ``ansys.mapdl.reader.mesh.Mesh`` class.

.. code:: python

    >>> from ansys.mapdl import reader as pymapdl_reader
    >>> from ansys.mapdl.reader import examples
    >>> rst = pymapdl_reader.read_binary(examples.rstfile)
    >>> print(rst.mesh)

.. code::

    ANSYS Mesh
      Number of Nodes:              321
      Number of Elements:           40
      Number of Element Types:      1
      Number of Node Components:    0
      Number of Element Components: 0


Which contains the following attributes:

.. autoclass:: ansys.mapdl.reader.mesh.Mesh
    :members:


Coordinate Systems
~~~~~~~~~~~~~~~~~~
Non-default coordinate systems are always saved to a MAPDL result
file.  The coordinate system is zero indexed and individual coordinate
systems can be accessed with:

.. code:: python

    >>> coord_idx = 12
    >>> result.geometry['coord systems'][coord_idx]
    {'transformation matrix': array([[ 0.0, -1.0,  0.0],
                                     [ 0.0,  0.0, -1.0],
                                     [ 1.0,  0.0,  0.0]]),
     'origin': array([0., 0., 0.]),
     'PAR1': 1.0,
     'PAR2': 1.0,
     'euler angles': array([ -0., -90.,  90.]),
     'theta singularity': 0.0,
     'phi singularity': 0.0,
     'type': 1,
     'reference num': 12}

A 4x4 transformation matrix can be constructed by concatenating the
transformation matrix and the origin into one array.  For example:

.. code:: python

    >>> cs = result.geometry['coord systems'][coord_idx]
    >>> trans = cs['transformation matrix']
    >>> origin = cs['origin']
    >>> bottom = np.zeros(4)
    >>> bottom[3] = 1
    >>> tmat = np.hstack((trans, origin.reshape(-1 ,1)))
    >>> tmat = np.vstack((tmat, bottom))

See ``parse_coordinate_system`` for more details regarding the
contents of the coordinate systems stored in the result file.


Accessing Solution Results
--------------------------
You can obtain detailed information using ``solution_info`` for each result:

.. code:: python

    # return a dictionary of solution info for the first result
    info = result.solution_info(0)

    for key in info:
        print(key, info[key])

This yields::

    timfrq 1.0
    lfacto 1.0
    lfactn 1.0
    cptime 50.9189941460218
    tref 0.0
    tunif 0.0
    tbulk 82.0
    volbase 0.0
    tstep 0.0
    __unused 0.0
    accel_x 0.0
    accel_y 0.0
    accel_z 0.0
    omega_v_x 0.0
    omega_v_y 0.0
    omega_v_z 100
    omega_a_x 0.0
    omega_a_y 0.0
    omega_a_z 0.0
    omegacg_v_x 0.0
    omegacg_v_y 0.0
    omegacg_v_z 0.0
    omegacg_a_x 0.0
    omegacg_a_y 0.0
    omegacg_a_z 0.0
    cgcent 0.0
    fatjack 0.0
    dval1 0.0
    pCnvVal 0.0


The DOF solution for an analysis for each node in the analysis can be
obtained using the code block below.  These results correspond to the
node numbers in the result file.  This array is sized by the number of
nodes by the number of degrees of freedom.

.. code:: python    

    # Return an array of results (nnod x dof)
    nnum, disp = result.nodal_solution(0) # uses 0 based indexing 
    
    # where nnum is the node numbers corresponding to the displacement results

    # The same results can be plotted using 
    result.plot_nodal_solution(0, 'x', label='Displacement') # x displacement

    # normalized displacement can be plotted by excluding the direction string
    result.plot_nodal_solution(0, label='Normalized')

Stress can be obtained as well using the below code.  The nodal stress
is computed in the same manner as MAPDL by averaging the stress
evaluated at that node for all attached elements.

.. code:: python
    
    # obtain the component node averaged stress for the first result
    # organized with one [Sx, Sy Sz, Sxy, Syz, Sxz] entry for each node
    nnum, stress = result.nodal_stress(0) # results in a np array (nnod x 6)

    # Display node averaged stress in x direction for result 6
    result.plot_nodal_stress(5, 'Sx')

    # Compute principal nodal stresses and plot SEQV for result 1
    nnum, pstress = result.principal_nodal_stress(0)
    result.plot_principal_nodal_stress(0, 'SEQV')

Element stress can be obtained using the following segment of code.
Ensure that the element results are expanded for a modal analysis
within ANSYS with::

    /SOLU
    MXPAND, ALL, , , YES

This block of code shows how you can access the non-averaged stresses
for the first result from a modal analysis.

.. code:: python
    
    from ansys.mapdl import reader as pymapdl_reader
    result = pymapdl_reader.read_binary('file.rst')
    estress, elem, enode = result.element_stress(0)

    
These stresses can be verified using MAPDL using:

.. code:: python

    >>> estress[0]
    [[ 1.0236604e+04 -9.2875127e+03 -4.0922625e+04 -2.3697146e+03
      -1.9239732e+04  3.0364934e+03]
     [ 5.9612605e+04  2.6905924e+01 -3.6161423e+03  6.6281304e+03
       3.1407712e+02  2.3195926e+04]
     [ 3.8178301e+04  1.7534495e+03 -2.5156013e+02 -6.4841372e+03
      -5.0892783e+03  5.2503605e+00]
     [ 4.9787645e+04  8.7987168e+03 -2.1928742e+04 -7.3025332e+03
       1.1294199e+04  4.3000205e+03]]

    >>> elem[0]
        32423

    >>> enode[0]
        array([ 9012,  7614,  9009, 10920], dtype=int32)

Which are identical to the results from MAPDL:

.. code::

  POST1:
  ESEL, S, ELEM, , 32423
  PRESOL, S

  ***** POST1 ELEMENT NODAL STRESS LISTING *****                                
 
  LOAD STEP=     1  SUBSTEP=     1                                             
   FREQ=    47.852      LOAD CASE=   0                                         
 
  THE FOLLOWING X,Y,Z VALUES ARE IN GLOBAL COORDINATES                         
 
  ELEMENT=   32423        SOLID187
    NODE    SX          SY          SZ          SXY         SYZ         SXZ     
    9012   10237.     -9287.5     -40923.     -2369.7     -19240.      3036.5    
    7614   59613.      26.906     -3616.1      6628.1      314.08      23196.    
    9009   38178.      1753.4     -251.56     -6484.1     -5089.3      5.2504    
   10920   49788.      8798.7     -21929.     -7302.5      11294.      4300.0    


Loading a Results from a Modal Analysis Result File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example reads in binary results from a modal analysis of a beam
from ANSYS.  This section of code does not rely on ``VTK`` and can be
used with only ``numpy`` installed.

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

    >>> print(freqs)
    [ 7366.49503969  7366.49503969 11504.89523664 17285.70459456
      17285.70459457 20137.19299035]
    
Get the 1st bending mode shape.  Results are ordered based on the
sorted node numbering.  Note that results are zero indexed.

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


Accessing Element Solution Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Individual element results for the entire solution can be accessed
using the ``element_solution_data`` method.  For example, to get the
volume of each element:

.. code:: python

    import numpy as np
    from ansys.mapdl import reader as pymapdl_reader

    rst = pymapdl_reader.read_binary('./file.rst')
    enum, edata = rst.element_solution_data(0, datatype='ENG')

    # output as a list, but can be viewed as an array since
    # the results for each element are the same size
    edata = np.asarray(edata)
    volume = edata[:, 0]


Animiating a Modal Solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Solutions from a modal analysis can be animated using
``animate_nodal_solution``.  For example:

.. code:: python

    from ansys.mapdl.reader import examples
    from ansys.mapdl import reader as pymapdl_reader

    result = pymapdl_reader.read_binary(examples.rstfile)
    result.animate_nodal_solution(3)


Plotting Nodal Results
~~~~~~~~~~~~~~~~~~~~~~
As the geometry of the model is contained within the result file, you
can plot the result without having to load any additional geometry.
Below, displacement for the first mode of the modal analysis beam is
plotted using ``VTK``.

Here, we plot the displacement of Mode 0 in the x direction:

.. code:: python
    
    result.plot_nodal_solution(0, 'x', label='Displacement')

.. image:: ../images/hexbeam_disp.png


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

.. image:: ../images/beam_stress.png

Nodal stress can also be generated non-interactively with:

.. code:: python

    result.plot_nodal_stress(5, 'Sx', cpos=cpos, screenshot=beam_stress.png,
                             window_size=[800, 600], interactive=False)

Animating a Modal Solution
~~~~~~~~~~~~~~~~~~~~~~~~~~
Mode shapes from a modal analysis can be animated using
``animate_nodal_solution``:

.. code:: python

    result.animate_nodal_solution(0)

If you wish to save the animation to a file, specify the
movie_filename and animate it with:

.. code:: python

    result.animate_nodal_solution(0, movie_filename='movie.mp4', cpos=cpos)

.. image:: ../images/beam_mode_shape.gif


Results from a Cyclic Analysis
------------------------------
The ``ansys-mapdl-reader`` module can load and display the results of
a cyclic analysis:

.. code:: python

    from ansys.mapdl import reader as pymapdl_reader

    # load the result file    
    result = pymapdl_reader.read_binary('rotor.rst')
    
You can reference the load step table and harmonic index tables by
printing the result header dictionary keys ``'ls_table'`` and
``'hindex'``:

.. code:: python

    >>> print(result.resultheader['ls_table'])
    # load step, sub step, cumulative index
    array([[ 1,  1,  1], 
           [ 1,  2,  2],
           [ 1,  3,  3],
           [ 1,  4,  4],
           [ 1,  5,  5],
           [ 2,  1,  6],

    >>> print(result.resultheader['hindex'])
    array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4,
           4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7], dtype=int32)

Where each harmonic index entry corresponds a cumulative index.  For
example, result number 11 is the first mode for the 2nd harmonic
index:

.. code:: python

    >>> result.resultheader['ls_table'][10] # Result 11 (using zero based indexing)
    array([ 3,  1, 11], dtype=int32)
    
    >>> result.resultheader['hindex'][10]
    2

Alternatively, the result number can be obtained by using:

.. code:: python

    >>> mode = 1
    >>> harmonic_index = 2
    >>> result.harmonic_index_to_cumulative(mode, harmonic_index)
    24

Using this indexing method, repeated modes are indexed by the same
mode index.  To access the other repeated mode, use a negative
harmonic index.  Should a result not exist, ``ansys-mapdl-reader``
will return which modes are available:

.. code:: python

    >>> mode = 1
    >>> harmonic_index = 20
    >>> result.harmonic_index_to_cumulative(mode, harmonic_index)
    Exception: Invalid mode for harmonic index 1
    Available modes: [0 1 2 3 4 5 6 7 8 9]

Results from a cyclic analysis require additional post processing to
be interperted correctly.  Mode shapes are stored within the result
file as unprocessed parts of the real and imaginary parts of a modal
solution.  ``ansys-mapdl-reader`` combines these values into a single
complex array and then returns the real result of that array.

.. code:: python

    >>> nnum, ms = result.nodal_solution(10) # mode shape of result 11
    >>> print(ms[:3])
    [[ 44.700, 45.953, 38.717]
     [ 42.339, 48.516, 52.475]
     [ 36.000, 33.121, 39.044]]

Sometimes it is necessary to determine the maximum displacement of a
mode.  To do so, return the complex solution with:

.. code:: python

    nnum, ms = result.nodal_solution(0, as_complex=True)
    norm = np.abs((ms*ms).sum(1)**0.5)
    idx = np.nanargmax(norm)
    ang = np.angle(ms[idx, 0])

    # rotate the solution by the angle of the maximum nodal response
    ms *= np.cos(ang) - 1j*np.sin(ang)

    # get only the real response
    ms = np.real(ms)
    
See ``help(result.nodal_solution)`` for more details.

The real displacement of the sector is always the real component of
the mode shape ``ms``, and this can be varied by multiplying the mode
shape by a complex value for a given phase.

The results of a single sector can be displayed as well using the
``plot_nodal_solution``

.. code:: python

    rnum = result.harmonic_index_to_cumulative(0, 2)
    result.plot_nodal_solution(rnum, label='Displacement', expand=False)
    
.. image:: ../images/rotor.jpg

The phase of the result can be changed by modifying the ``phase``
option.  See ``help(result.plot_nodal_solution)`` for details on its
implementation.


Exporting to ParaView
---------------------
ParaView is a visualization application that can be used for rapid
generation of plots and graphs using VTK through a GUI.
``ansys-mapdl-reader`` can translate the MAPDL result files to
ParaView compatible files containing the geometry and nodal results
from the analysis:

.. code:: python

    from ansys.mapdl import reader as pymapdl_reader
    from ansys.mapdl.reader import examples

    # load example beam result file
    result = pymapdl_reader.read_binary(examples.rstfile)
    
    # save as a binary vtk xml file
    result.save_as_vtk('beam.vtu')

The vtk xml file can now be loaded using ParaView.  This screenshot
shows the nodal displacement of the first result from the result file
plotted within `ParaView <https://www.paraview.org/>`_.  Within the
vtk file are two point arrays (``NodalResult`` and ``nodal_stress``)
for each result in the result file.  The nodal result values will
depend on the analysis type, while nodal stress will always be the
node average stress in the Sx, Sy Sz, Sxy, Syz, and Sxz directions.

.. image:: ../images/paraview.jpg


Result Object Methods
---------------------
.. autoclass:: ansys.mapdl.reader.rst.Result
    :members:

.. autoclass:: ansys.mapdl.reader.cyclic_reader.CyclicResult
    :members:
