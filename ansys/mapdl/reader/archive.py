"""Module to read ANSYS ASCII block formatted CDB files"""
import io
import os
import logging
from functools import wraps

import numpy as np
from vtk import (VTK_TETRA, VTK_QUADRATIC_TETRA, VTK_PYRAMID,
                 VTK_QUADRATIC_PYRAMID, VTK_WEDGE, VTK_QUADRATIC_WEDGE,
                 VTK_HEXAHEDRON, VTK_QUADRATIC_HEXAHEDRON, VTK_TRIANGLE, VTK_QUAD,
                 VTK_QUADRATIC_TRIANGLE, VTK_QUADRATIC_QUAD)
import vtk
import pyvista as pv

from ansys.mapdl.reader import _reader, _archive
from ansys.mapdl.reader.misc import vtk_cell_info, chunks
from ansys.mapdl.reader.mesh import Mesh
from ansys.mapdl.reader.cell_quality import quality

VTK9 = vtk.vtkVersion().GetVTKMajorVersion() >= 9

log = logging.getLogger(__name__)
log.setLevel('CRITICAL')


class Archive(Mesh):
    """Read a blocked ANSYS archive file or input file.

    Reads a blocked CDB file and optionally parses it to a vtk grid.
    This can be used to read in files written from MAPDL using the
    ``CDWRITE`` command or input files (``'.dat'``) files written from
    ANSYS Workbench.

    Write the archive file using ``CDWRITE, DB, archive.cdb``

    Parameters
    ----------
    filename : string
        Filename of block formatted cdb file

    read_parameters : bool, optional
        Optionally read parameters from the archive file.  Default
        ``False``.

    parse_vtk : bool, optional
        When ``True``, parse the raw data into to VTK format.

    force_linear : bool, optional
        This parser creates quadratic elements if available.  Set
        this to True to always create linear elements.  Defaults
        to False.

    allowable_types : list, optional
        Allowable element types.  Defaults to all valid element
        types in ``ansys.mapdl.reader.elements.valid_types``

        See ``help(ansys.mapdl.reader.elements)`` for available element types.

    null_unallowed : bool, optional
        Elements types not matching element types will be stored
        as empty (null) elements.  Useful for debug or tracking
        element numbers.  Default False.

    verbose : bool, optional
        Print out each step when reading the archive file.  Used for
        debug purposes and defaults to ``False``.

    name : str, optional
        Internally used parameter used to have a custom ``__repr__``.

    Examples
    --------
    >>> from ansys.mapdl import reader as pymapdl_reader
    >>> from ansys.mapdl.reader import examples
    >>> hex_beam = pymapdl_reader.Archive(examples.hexarchivefile)
    >>> print(hex_beam)
    ANSYS Archive File HexBeam.cdb
      Number of Nodes:              40
      Number of Elements:           321
      Number of Element Types:      1
      Number of Node Components:    2
      Number of Element Components: 2

    Print the node array

    >>> hex_beam.nodes
    array([[0.  , 0.  , 0.  ],
           [1.  , 0.  , 0.  ],
           [0.25, 0.  , 0.  ],
           ...,
           [0.75, 0.5 , 3.5 ],
           [0.75, 0.5 , 4.  ],
           [0.75, 0.5 , 4.5 ]])

    Read an ANSYS workbench input file

    >>> my_archive = pymapdl_reader.Archive('C:\\Users\\user\\stuff.dat')

    Notes
    -----
    This class only reads EBLOCK records with SOLID records.  For
    example, the record ``EBLOCK,19,SOLID,,3588`` will be read, but
    ``EBLOCK,10,,,3588`` will not be read.  Generally, MAPDL will only
    write SOLID records and Mechanical Workbench may write SOLID
    records.  These additional records will be ignored.
    """

    def __init__(self, filename, read_parameters=False,
                 parse_vtk=True, force_linear=False,
                 allowable_types=None, null_unallowed=False,
                 verbose=False, name=''):
        """Initializes an instance of the archive class."""
        self._read_parameters = read_parameters
        self._filename = filename
        self._name = name
        self._raw = _reader.read(filename, read_parameters=read_parameters,
                                 debug=verbose)
        super().__init__(self._raw['nnum'],
                         self._raw['nodes'],
                         self._raw['elem'],
                         self._raw['elem_off'],
                         self._raw['ekey'],
                         node_comps=self._raw['node_comps'],
                         elem_comps=self._raw['elem_comps'],
                         rdat=self._raw['rdat'],
                         rnum=self._raw['rnum'],
                         keyopt=self._raw['keyopt'])

        self._allowable_types = allowable_types
        self._force_linear = force_linear
        self._null_unallowed = null_unallowed

        if parse_vtk:
            self._grid = self._parse_vtk(allowable_types,
                                         force_linear, null_unallowed)

    @property
    def raw(self):  # pragma: no cover
        raise AttributeError('The `raw` attribute has been depreciated.  Access'
                             ' the values directly from the archive object.\n\n'
                             '    Instead of:\n'
                             '    archive.raw["nodes"]\n'
                             '    \n'
                             '    Use\n'
                             '    print(archive.nodes)')

    @property
    def parameters(self):
        """Parameters stored in the archive file

        Examples
        --------
        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> from ansys.mapdl.reader import examples
        >>> archive = pymapdl_reader.Archive(examples.hexarchivefile,
                                             read_parameters=True)
        >>> archive.parameters
        {}
        """
        if not self._read_parameters:
            raise AttributeError('No parameters read.  Read the archive again '
                                 ' with ``read_parameters=True``')
        return self._raw['parameters']

    def __repr__(self):
        if self._name:
            txt = 'MAPDL %s\n' % self._name
        else:
            basename = os.path.basename(self._filename)
            txt = 'ANSYS Archive File %s\n' % basename

        txt += '  Number of Nodes:              %d\n' % len(self.nnum)
        txt += '  Number of Elements:           %d\n' % len(self.enum)
        txt += '  Number of Element Types:      %d\n' % len(self.ekey)
        txt += '  Number of Node Components:    %d\n' % len(self.node_components)
        txt += '  Number of Element Components: %d\n' % len(self.element_components)
        return txt

    @property
    def grid(self):
        """``vtk.UnstructuredGrid`` of the archive file.

        Examples
        --------
        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> from ansys.mapdl.reader import examples
        >>> archive = pymapdl_reader.Archive(examples.hexarchivefile)
        >>> archive.grid
        UnstructuredGrid (0x7ffa237f08a0)
          N Cells:      40
          N Points:     321
          X Bounds:     0.000e+00, 1.000e+00
          Y Bounds:     0.000e+00, 1.000e+00
          Z Bounds:     0.000e+00, 5.000e+00
          N Arrays:     13
        """
        if self._grid is None:  # parse the grid using the cached parameters
            self._grid = self._parse_vtk(self._allowable_types,
                                         self._force_linear,
                                         self._null_unallowed)
        return self._grid

    @property
    def quality(self):
        """Minimum scaled jacobian cell quality.

        Negative values indicate invalid cells while positive values
        indicate valid cells.  Varies between -1 and 1.

        Examples
        --------
        >>> from ansys.mapdl import reader as pymapdl_reader
        >>> from ansys.mapdl.reader import examples
        >>> archive = pymapdl_reader.Archive(examples.hexarchivefile)
        >>> archive.quality
        array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
               1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
               1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])
        """
        if self._grid is None:  # pragma: no cover
            raise AttributeError('Archive must be parsed as a vtk grid.\n'
                                 'Set `parse_vtk=True`')
        return quality(self._grid)

    @wraps(pv.plot)
    def plot(self, *args, **kwargs):
        """Plot the mesh"""
        if self._grid is None:  # pragma: no cover
            raise AttributeError('Archive must be parsed as a vtk grid.\n'
                                 'Set `parse_vtk=True`')
        kwargs.setdefault('color', 'w')
        kwargs.setdefault('show_edges', True)
        self.grid.plot(*args, **kwargs)


def save_as_archive(filename, grid, mtype_start=1, etype_start=1,
                    real_constant_start=1, mode='w', nblock=True,
                    enum_start=1, nnum_start=1,
                    include_etype_header=True, reset_etype=False,
                    allow_missing=True, include_surface_elements=True,
                    include_solid_elements=True, include_components=True):
    """Writes FEM as an ANSYS APDL archive file.  This function
    supports the following element types:

        - ``vtk.VTK_TETRA``
        - ``vtk.VTK_QUADRATIC_TETRA``
        - ``vtk.VTK_PYRAMID``
        - ``vtk.VTK_QUADRATIC_PYRAMID``
        - ``vtk.VTK_WEDGE``
        - ``vtk.VTK_QUADRATIC_WEDGE``
        - ``vtk.VTK_HEXAHEDRON``
        - ``vtk.VTK_QUADRATIC_HEXAHEDRON``
        - ``vtk.VTK_TRIANGLE``
        - ``vtk.VTK_QUAD``

    Will automatically renumber nodes and elements if the FEM does not
    contain ANSYS node or element numbers.  Node numbers are stored as
    a point array "ansys_node_num", and cell numbers are stored as
    cell array "ansys_elem_num".

    Parameters
    ----------
    filename : str
       Filename to write archive file.

    grid : vtk.UnstructuredGrid
        VTK UnstructuredGrid to convert to an APDL archive file.
        PolyData will automatically be converted to an unstructured
        mesh.

    mtype_start : int, optional
        Material number to assign to elements.  Can be set manually by
        adding the cell array "mtype" to the unstructured grid.

    etype_start : int, optional
        Starting element type number.  Can be manually set by adding
        the cell array "ansys_etype" to the unstructured grid.

    real_constant_start : int, optional
        Starting real constant to assign to unset cells.  Can be
        manually set by adding the cell array "ansys_real_constant" to
        the unstructured grid.

    mode : str, optional
        File mode.  See ``help(open)``

    nblock : bool, optional
        Write node block when writing archive file.

    enum_start : int, optional
        Starting element number to assign to unset cells.  Can be
        manually set by adding the cell array "ansys_elem_num" to the
        unstructured grid.

    nnum_start : int, optional
        Starting element number to assign to unset points.  Can be
        manually set by adding the point array "ansys_node_num" to the
        unstructured grid.

    include_etype_header : bool, optional
        For each element type, includes element type command
        (e.g. "ET, 1, 186") in the archive file.

    reset_etype : bool, optional
        Resets element type.  Element types will automatically be
        determined by the shape of the element (i.e. quadradic
        tetrahedrals will be saved as SOLID187, linear hexahedrals as
        SOLID185).  Default True.

    include_surface_elements : bool, optional
        Includes surface elements when writing the archive file and
        saves them as SHELL181.

    include_solid_elements : bool, optional
        Includes solid elements when writing the archive file and
        saves them as SOLID185, SOLID186, or SOLID187.

    include_components : bool, optional
        Writes note components to file.  Node components must be
        stored within the unstructured grid as uint8 or bool arrays.

    Examples
    --------
    Write a VTK of pyvista UnstructuredGrid to archive.cdb

    >>> from ansys.mapdl import reader as pymapdl_reader
    >>> from pyvista import examples
    >>> grid = examples.load_hexbeam()
    >>> pymapdl_reader.save_as_archive('archive.cdb', grid)

    """

    if isinstance(grid, pv.PolyData):
        grid = grid.cast_to_unstructured_grid()

    if not isinstance(grid, vtk.vtkUnstructuredGrid):
        raise TypeError('``grid`` argument must be an UnstructuredGrid')

    allowable = []
    if include_solid_elements:
        allowable.extend([VTK_TETRA,
                          VTK_QUADRATIC_TETRA,
                          VTK_PYRAMID,
                          VTK_QUADRATIC_PYRAMID,
                          VTK_WEDGE,
                          VTK_QUADRATIC_WEDGE,
                          VTK_HEXAHEDRON,
                          VTK_QUADRATIC_HEXAHEDRON])

    if include_surface_elements:
        allowable.extend([VTK_TRIANGLE,
                          VTK_QUAD])
                          # VTK_QUADRATIC_TRIANGLE,
                          # VTK_QUADRATIC_QUAD

    # extract allowable cell types
    mask = np.in1d(grid.celltypes, allowable)
    grid = grid.extract_cells(mask)

    header = '/PREP7\n'

    # node numbers
    if 'ansys_node_num' in grid.point_arrays:
        nodenum = grid.point_arrays['ansys_node_num']
    else:
        log.info('No ANSYS node numbers set in input. Adding default range')
        nodenum = np.arange(1, grid.number_of_points + 1, dtype=np.int32)

    if np.any(nodenum == -1):
        if not allow_missing:
            raise Exception('Missing node numbers.  Exiting due "allow_missing=False"')
        start_num = nodenum.max() + 1
        if nnum_start > start_num:
            start_num = nnum_start
        nadd = np.sum(nodenum == -1)
        end_num = start_num + nadd
        log.info('FEM missing some node numbers.  Adding node numbering '
                 'from %d to %d', start_num, end_num)
        nodenum[nodenum == -1] = np.arange(start_num, end_num, dtype=np.int32)

    # element block
    ncells = grid.number_of_cells
    if 'ansys_elem_num' in grid.cell_arrays:
        enum = grid.cell_arrays['ansys_elem_num']
    else:
        if not allow_missing:
            raise Exception('Missing node numbers.  Exiting due "allow_missing=False"')
        log.info('No ANSYS element numbers set in input.  ' +
                 'Adding default range starting from %d' % enum_start)
        enum = np.arange(1, ncells + 1, dtype=np.int32)

    if np.any(enum == -1):
        if not allow_missing:
            raise Exception('-1 encountered in "ansys_elem_num".\n'
                            'Exiting due "allow_missing=False"')

        start_num = enum.max() + 1
        if enum_start > start_num:
            start_num = enum_start
        nadd = np.sum(enum == -1)
        end_num = start_num + nadd
        log.info('FEM missing some cell numbers.  Adding numbering ' +
                 'from %d to %d' % (start_num, end_num))
        enum[enum == -1] = np.arange(start_num, end_num, dtype=np.int32)

    # material type
    if 'ansys_material_type' in grid.cell_arrays:
        mtype = grid.cell_arrays['ansys_material_type']
    else:
        log.info('No ANSYS element numbers set in input.  ' +
                 'Adding default range starting from %d', mtype_start)
        mtype = np.arange(1, ncells + 1, dtype=np.int32)

    if np.any(mtype == -1):
        log.info('FEM missing some material type numbers.  Adding...')
        mtype[mtype == -1] = mtype_start

    # real constant
    if 'ansys_real_constant' in grid.cell_arrays:
        rcon = grid.cell_arrays['ansys_real_constant']
    else:
        log.info('No ANSYS element numbers set in input.  ' +
                 'Adding default range starting from %d', real_constant_start)
        rcon = np.arange(1, ncells + 1, dtype=np.int32)

    if np.any(rcon == -1):
        log.info('FEM missing some material type numbers.  Adding...')
        rcon[rcon == -1] = real_constant_start

    # element type
    invalid = False
    if 'ansys_etype' in grid.cell_arrays and not reset_etype:
        missing = False
        typenum = grid.cell_arrays['ansys_elem_type_num']
        etype = grid.cell_arrays['ansys_etype']
        if np.any(etype == -1):
            log.warning('Some elements are missing element type numbers.')
            invalid = True

        if include_etype_header and not invalid:
            _, ind = np.unique(etype, return_index=True)
            for idx in ind:
                header += 'ET, %d, %d\n' % (etype[idx], typenum[idx])
    else:
        missing = True

    # check if valid
    if not missing:
        mask = grid.celltypes < 20
        if np.any(grid.cell_arrays['ansys_elem_type_num'][mask] == 186):
            invalid = True
            log.warning('Invalid ANSYS element types.')

    if invalid or missing:
        if not allow_missing:
            raise Exception('Invalid or missing data in "ansys_elem_type_num"' +
                            ' or "ansys_etype".  Exiting due "allow_missing=False"')

        log.info('No ANSYS element type or invalid data input.  ' +
                 'Adding default range starting from %d' % etype_start)

        etype = np.empty(grid.number_of_cells, np.int32)
        etype_185 = etype_start + 2
        etype[grid.celltypes == VTK_TETRA] = etype_185
        etype[grid.celltypes == VTK_HEXAHEDRON] = etype_185
        etype[grid.celltypes == VTK_WEDGE] = etype_185
        etype[grid.celltypes == VTK_PYRAMID] = etype_185

        etype_186 = etype_start
        etype[grid.celltypes == VTK_QUADRATIC_HEXAHEDRON] = etype_186
        etype[grid.celltypes == VTK_QUADRATIC_WEDGE] = etype_186
        etype[grid.celltypes == VTK_QUADRATIC_PYRAMID] = etype_186

        etype_187 = etype_start + 1
        etype[grid.celltypes == VTK_QUADRATIC_TETRA] = etype_187

        # Surface elements
        etype_181 = etype_start + 3
        etype[grid.celltypes == VTK_TRIANGLE] = etype_181
        etype[grid.celltypes == VTK_QUAD] = etype_181

        typenum = np.empty_like(etype)
        typenum[etype == etype_185] = 185
        typenum[etype == etype_186] = 186
        typenum[etype == etype_187] = 187
        typenum[etype == etype_181] = 181

        header += 'ET, %d, 185\n' % etype_185
        header += 'ET, %d, 186\n' % etype_186
        header += 'ET, %d, 187\n' % etype_187
        header += 'ET, %d, 181\n' % etype_181

    # number of nodes written per element
    elem_nnodes = np.empty(etype.size, np.int32)
    elem_nnodes[typenum == 181] = 4
    elem_nnodes[typenum == 185] = 8
    elem_nnodes[typenum == 186] = 20
    elem_nnodes[typenum == 187] = 10

    # write the EBLOCK
    with open(str(filename), mode) as f:
        f.write(header)

    if not isinstance(filename, str):
        filename = str(filename)
    write_nblock(filename, nodenum, grid.points, mode='a')

    # write remainder of eblock
    cells, offset = vtk_cell_info(grid, shift_offset=False)
    _archive.py_write_eblock(filename,
                             enum,
                             etype,
                             mtype,
                             rcon,
                             elem_nnodes,
                             cells,
                             offset,
                             grid.celltypes,
                             typenum,
                             nodenum,
                             VTK9,
                             mode='a')

    if include_components:
        with open(filename, 'a') as fid:

            # write node components
            for node_key in grid.point_arrays:
                arr = grid.point_arrays[node_key]
                if arr.dtype in [np.uint8, np.bool]:
                    items = nodenum[arr.view(np.bool)]
                    write_cmblock(fid, items, node_key, 'NODE')

            # write element components
            for node_key in grid.cell_arrays:
                arr = grid.cell_arrays[node_key]
                if arr.dtype in [np.uint8, np.bool]:
                    items = enum[arr.view(np.bool)]
                    write_cmblock(fid, items, node_key, 'ELEMENT')


def write_nblock(filename, node_id, pos, angles=None, mode='w'):
    """Writes nodes and node angles to file.

    Parameters
    ----------
    filename : str or file handle
        Filename to write node block to.

    node_id : np.ndarray
        ANSYS node numbers.

    pos : np.ndarray
        Node coordinates.

    angles : np.ndarray, optional
        Writes the node angles for each node when included.

    mode : str, optional
        Write mode.  Default ``'w'``.
    """
    assert pos.ndim == 2 and pos.shape[1] == 3, 'Invalid position array'
    if angles is not None:
        assert angles.ndim == 2 and angles.shape[1] == 3, 'Invalid angle array'

    if node_id.dtype != np.int32:
        node_id = node_id.astype(np.int32)

    # node array must be sorted
    # note, this is sort check is most suited for pre-sorted arrays
    # see https://stackoverflow.com/questions/3755136/
    if not np.all(node_id[:-1] <= node_id[1:]):
        sidx = np.argsort(node_id)
        node_id = node_id[sidx]
        pos = pos[sidx]

    if angles is not None:
        if pos.dtype == np.float32:
            if angles.dtype != pos.dtype:
                angles = angles.astype(pos.dtype)
            _archive.py_write_nblock_float(filename,
                                           node_id,
                                           node_id[-1],
                                           pos,
                                           angles,
                                           mode)
        else:
            _archive.py_write_nblock(filename,
                                     node_id,
                                     node_id[-1],
                                     pos,
                                     angles,
                                     mode)

    else:
        if pos.dtype == np.float32:
            _archive.py_write_nblock_float(filename,
                                           node_id,
                                           node_id[-1],
                                           pos,
                                           np.empty((0, 0), dtype=np.float32),
                                           mode)
        else:
            _archive.py_write_nblock(filename,
                                     node_id,
                                     node_id[-1],
                                     pos,
                                     np.empty((0, 0), dtype=np.float64),
                                     mode)


def write_cmblock(filename, items, comp_name, comp_type,
                  digit_width=10, mode='w'):
    """Writes a component block (CMBLOCK) to a file.

    Parameters
    ----------
    filename : str or file handle
        File to write CMBLOCK component to.

    items : list or np.ndarray
        Element or node numbers to write.

    comp_name : str
        Name of the component

    comp_type : str
        Component type to write.  Should be either 'ELEMENT' or 'NODE'.

    digit_width : int, optional
        Default 10

    mode : str, optional
        Write mode.  Default ``'w'``.
    """
    comp_name = comp_name.upper()
    comp_type = comp_type.upper()
    if comp_type.upper() not in ['ELEMENT', 'NODE']:
        raise ValueError("`comp_type` must be either 'ELEMENT' or 'NODE'")

    opened_file = False
    if isinstance(filename, io.TextIOBase):
        fid = filename
    else:
        fid = open(filename, mode)
        opened_file = True

    if not isinstance(items, np.ndarray):
        items = np.array(items, dtype=np.int32)
    elif items.dtype != np.int32:
        items = items.astype(np.int32)

    # All this python writing could be a bottleneck for non-contiguous CMBLOCKs.
    # consider cythonizing this in the future
    cmblock_items = _archive.cmblock_items_from_array(items)
    nitems = len(cmblock_items)
    print(f'CMBLOCK,{comp_name},{comp_type},{nitems:8d}', file=fid)
    print(f'(8i{digit_width})', file=fid)
    digit_formatter = f'%{digit_width}d'

    # use np savetxt here as it's faster than looping through and
    # writing each line.
    # nearest multiple of 8
    up_to = len(cmblock_items) % 8
    np.savetxt(fid, cmblock_items[:-up_to].reshape(-1, 8), digit_formatter*8)

    # write the final line
    chunk = cmblock_items[-up_to:]
    print(''.join([digit_formatter] * len(chunk)) % tuple(chunk), file=fid)
    print('', file=fid)

    if opened_file:
        fid.close()
