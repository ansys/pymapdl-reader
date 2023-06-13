import os
import pathlib

import numpy as np
import pytest
import pyvista as pv
from pyvista import CellType
from pyvista import examples as pyvista_examples

from ansys.mapdl import reader as pymapdl_reader
from ansys.mapdl.reader import _archive, archive, examples

LINEAR_CELL_TYPES = [
    CellType.TETRA,
    CellType.PYRAMID,
    CellType.WEDGE,
    CellType.HEXAHEDRON,
]
QUADRATIC_CELL_TYPES = [
    CellType.QUADRATIC_TETRA,
    CellType.QUADRATIC_PYRAMID,
    CellType.QUADRATIC_WEDGE,
    CellType.QUADRATIC_HEXAHEDRON,
]

TEST_PATH = os.path.dirname(os.path.abspath(__file__))
TESTFILES_PATH = os.path.join(TEST_PATH, "test_data")
TESTFILES_PATH_PATHLIB = pathlib.Path(TESTFILES_PATH)
DAT_FILE = os.path.join(TESTFILES_PATH, "Panel_Transient.dat")


def proto_cmblock(array):
    """prototype cmblock code"""
    items = np.zeros_like(array)
    items[0] = array[0]

    c = 1
    in_list = False
    for i in range(array.size - 1):
        # check if part of a range
        if array[i + 1] - array[i] == 1:
            in_list = True
        elif array[i + 1] - array[i] > 1:
            if in_list:
                items[c] = -array[i]
                c += 1
                items[c] = array[i + 1]
                c += 1
            else:
                items[c] = array[i + 1]
                c += 1
            in_list = False

    # check if we've ended on a list
    # catch if last item is part of a list
    if items[c - 1] != abs(array[-1]):
        items[c] = -array[i + 1]
        c += 1

    return items[:c]


@pytest.fixture()
def pathlib_archive():
    filename = TESTFILES_PATH_PATHLIB / "ErnoRadiation.cdb"
    return pymapdl_reader.Archive(filename)


@pytest.fixture()
def hex_archive():
    return pymapdl_reader.Archive(examples.hexarchivefile)


@pytest.fixture(scope="module")
def all_solid_cells_archive():
    return pymapdl_reader.Archive(os.path.join(TESTFILES_PATH, "all_solid_cells.cdb"))


@pytest.fixture(scope="module")
def all_solid_cells_archive_linear():
    return pymapdl_reader.Archive(
        os.path.join(TESTFILES_PATH, "all_solid_cells.cdb"), force_linear=True
    )


@pytest.mark.parametrize(
    "array",
    (
        np.arange(1, 10, dtype=np.int32),
        np.array([1, 5, 10, 20, 40, 80], dtype=np.int32),
        np.array([1, 2, 3, 10, 20, 40, 51, 52, 53], dtype=np.int32),
        np.array([1, 2, 3, 10, 20, 40], dtype=np.int32),
        np.array([10, 20, 40, 50, 51, 52], dtype=np.int32),
    ),
)
def test_cython_cmblock(array):
    """Simply verify it's identical to the prototype python code"""
    assert np.allclose(proto_cmblock(array), _archive.cmblock_items_from_array(array))


def test_load_dat():
    arch = pymapdl_reader.Archive(DAT_FILE, read_parameters=True)
    assert arch.n_node == 1263  # through inspection of the dat file
    assert arch.n_elem == 160  # through inspection of the dat file
    assert "Panelflattern" in arch.parameters["_wb_userfiles_dir"]


def test_repr(hex_archive):
    assert "%s" % hex_archive.n_node in str(hex_archive)
    assert "%s" % hex_archive.n_elem in str(hex_archive)


def test_read_mesh200():
    archive = pymapdl_reader.Archive(os.path.join(TESTFILES_PATH, "mesh200.cdb"))
    assert archive.grid.n_cells == 1000


def test_archive_init(hex_archive):
    assert isinstance(hex_archive._raw, dict)
    assert isinstance(hex_archive.grid, pv.UnstructuredGrid)


def test_parse_vtk(hex_archive):
    grid = hex_archive.grid
    assert grid.points.size
    assert grid.cells.size
    assert "ansys_node_num" in grid.point_data
    assert np.all(hex_archive.quality > 0)

    with pytest.raises(TypeError):
        hex_archive._parse_vtk(allowable_types=-1)

    with pytest.raises(TypeError):
        hex_archive._parse_vtk(allowable_types=3.0)


def test_invalid_archive(tmpdir, hex_archive):
    nblock_filename = str(tmpdir.mkdir("tmpdir").join("nblock.cdb"))
    pymapdl_reader.write_nblock(nblock_filename, hex_archive.nnum, hex_archive.nodes)

    archive = pymapdl_reader.Archive(nblock_filename)
    assert archive.grid is None


def test_write_angle(tmpdir, hex_archive):
    nblock_filename = str(tmpdir.mkdir("tmpdir").join("nblock.cdb"))
    pymapdl_reader.write_nblock(
        nblock_filename, hex_archive.nnum, hex_archive.nodes, hex_archive.node_angles
    )

    archive = pymapdl_reader.Archive(nblock_filename, parse_vtk=False)
    assert np.allclose(archive.nodes, hex_archive.nodes)


def test_missing_midside():
    allowable_types = [45, 95, 185, 186, 92, 187]
    archive_file = os.path.join(TESTFILES_PATH, "mixed_missing_midside.cdb")
    archive = pymapdl_reader.Archive(archive_file, allowable_types=allowable_types)

    assert (archive.quality > 0.0).all()
    assert not np.any(archive.grid.celltypes == CellType.TETRA)


def test_missing_midside_write(tmpdir):
    allowable_types = [45, 95, 185, 186, 92, 187]
    archive_file = os.path.join(TESTFILES_PATH, "mixed_missing_midside.cdb")
    archive = pymapdl_reader.Archive(archive_file, allowable_types=allowable_types)

    filename = str(tmpdir.join("tmp.cdb"))
    with pytest.raises(RuntimeError, match="Unsupported element types"):
        pymapdl_reader.save_as_archive(filename, archive.grid, exclude_missing=True)

    pymapdl_reader.save_as_archive(
        filename, archive.grid, exclude_missing=True, reset_etype=True
    )
    archive_new = pymapdl_reader.Archive(filename)


def test_writehex(tmpdir, hex_archive):
    filename = str(tmpdir.mkdir("tmpdir").join("tmp.cdb"))
    pymapdl_reader.save_as_archive(filename, hex_archive.grid)
    archive_new = pymapdl_reader.Archive(filename)
    assert np.allclose(hex_archive.grid.points, archive_new.grid.points)
    assert np.allclose(
        hex_archive.grid.cell_connectivity,
        archive_new.grid.cell_connectivity,
    )

    for node_component in hex_archive.node_components:
        assert np.allclose(
            hex_archive.node_components[node_component],
            archive_new.node_components[node_component],
        )

    for element_component in hex_archive.element_components:
        assert np.allclose(
            hex_archive.element_components[element_component],
            archive_new.element_components[element_component],
        )


def test_write_voxel(tmpdir):
    filename = str(tmpdir.join("tmp.cdb"))
    grid = pv.UniformGrid(dimensions=(10, 10, 10))
    pymapdl_reader.save_as_archive(filename, grid)

    archive = pymapdl_reader.Archive(filename)
    assert np.allclose(archive.grid.points, grid.points)
    assert np.allclose(archive.grid.point_data["ansys_node_num"], range(1, 1001))
    assert archive.grid.n_cells, grid.n_cells


def test_writesector(tmpdir):
    archive = pymapdl_reader.Archive(examples.sector_archive_file)
    filename = str(tmpdir.mkdir("tmpdir").join("tmp.cdb"))
    pymapdl_reader.save_as_archive(filename, archive.grid)
    archive_new = pymapdl_reader.Archive(filename)

    assert np.allclose(archive.grid.points, archive_new.grid.points)
    assert np.allclose(archive.grid.cells, archive_new.grid.cells)


def test_writehex_missing_elem_num(tmpdir, hex_archive):
    grid = hex_archive.grid
    grid.cell_data["ansys_elem_num"][:10] = -1
    grid.cell_data["ansys_etype"] = np.ones(grid.number_of_cells) * -1
    grid.cell_data["ansys_elem_type_num"] = np.ones(grid.number_of_cells) * -1

    filename = str(tmpdir.mkdir("tmpdir").join("tmp.cdb"))
    pymapdl_reader.save_as_archive(filename, grid)
    archive_new = pymapdl_reader.Archive(filename)

    assert np.allclose(hex_archive.grid.points, archive_new.grid.points)
    assert np.allclose(hex_archive.grid.cells, archive_new.grid.cells)


def test_writehex_missing_node_num(tmpdir, hex_archive):
    hex_archive.grid.point_data["ansys_node_num"][:-1] = -1

    filename = str(tmpdir.mkdir("tmpdir").join("tmp.cdb"))
    pymapdl_reader.save_as_archive(filename, hex_archive.grid)
    archive_new = pymapdl_reader.Archive(filename)
    assert np.allclose(hex_archive.grid.points.shape, archive_new.grid.points.shape)
    assert np.allclose(hex_archive.grid.cells.size, archive_new.grid.cells.size)


def test_write_non_ansys_grid(tmpdir):
    grid = pv.UnstructuredGrid(pyvista_examples.hexbeamfile)
    del grid.point_data["sample_point_scalars"]
    del grid.cell_data["sample_cell_scalars"]
    archive_file = str(tmpdir.mkdir("tmpdir").join("tmp.cdb"))
    pymapdl_reader.save_as_archive(archive_file, grid)


def test_read_complex_archive(all_solid_cells_archive):
    nblock_expected = np.array(
        [
            [3.7826539829200e00, 1.2788958692644e00, -1.0220880953640e00],
            [3.7987359490873e00, 1.2312085780780e00, -1.0001885444969e00],
            [3.8138798206653e00, 1.1833200772896e00, -9.7805743587145e-01],
            [3.7751258193793e00, 1.2956563072306e00, -9.9775569295981e-01],
            [3.7675976558386e00, 1.3124167451968e00, -9.7342329055565e-01],
            [3.8071756567432e00, 1.2018089624856e00, -9.5159140433025e-01],
            [3.8004714928212e00, 1.2202978476816e00, -9.2512537278904e-01],
            [3.7840345743299e00, 1.2663572964392e00, -9.4927433167235e-01],
            [3.8682501483615e00, 1.4211343558710e00, -9.2956245308371e-01],
            [3.8656154427804e00, 1.4283573726940e00, -9.3544082975315e-01],
            [3.8629807371994e00, 1.4355803895169e00, -9.4131920642259e-01],
            [3.8698134427618e00, 1.4168612083433e00, -9.3457292477788e-01],
            [3.8645201728196e00, 1.4314324609914e00, -9.4526873324423e-01],
            [3.8713767371621e00, 1.4125880608155e00, -9.3958339647206e-01],
            [3.8687181728010e00, 1.4199362966407e00, -9.4440082826897e-01],
            [3.8660596084399e00, 1.4272845324660e00, -9.4921826006588e-01],
            [3.7847463501820e00, 1.2869612289286e00, -1.0110875234148e00],
            [3.7882161293470e00, 1.2952473975570e00, -1.0006326084202e00],
            [3.7840036708439e00, 1.3089808408341e00, -9.8189659453120e-01],
            [3.7736944340897e00, 1.3175655146540e00, -9.6829193559890e-01],
            [3.7797912123408e00, 1.3227142841112e00, -9.6316058064216e-01],
            [3.8163322819008e00, 1.1913589544053e00, -9.6740419078720e-01],
            [3.8046827481496e00, 1.2474593204382e00, -9.7922600135387e-01],
            [3.8202228218151e00, 1.1995824283636e00, -9.5733187068101e-01],
            [3.9797161316330e00, 2.5147820926190e-01, -5.1500799817626e-01],
            [3.9831382922541e00, 2.0190980565891e-01, -5.0185526897444e-01],
            [3.9810868976408e00, 2.3910377061737e-01, -5.4962360790281e-01],
            [3.9772930845240e00, 2.8865001362748e-01, -5.6276585706615e-01],
            [3.9816265976187e00, 2.1428739259987e-01, -4.6723916677654e-01],
            [3.9839413943097e00, 1.8949722823843e-01, -5.3648152416530e-01],
            [3.7962006776348e00, 1.2764624207283e00, -9.3931008487698e-01],
            [3.8126101429289e00, 1.2302105573453e00, -9.1545958911180e-01],
            [3.8065408178751e00, 1.2252542025135e00, -9.2029248095042e-01],
            [3.8164164823720e00, 1.2148964928545e00, -9.3639572989640e-01],
            [3.8972892823450e00, 2.7547119775919e-01, -5.6510422311694e-01],
            [3.9015993648189e00, 2.0235606714652e-01, -4.6987255385930e-01],
            [3.9023812010290e00, 1.7705558022279e-01, -5.3881795411458e-01],
            [3.9019902829240e00, 1.8970582368465e-01, -5.0434525398694e-01],
            [3.8998352416870e00, 2.2626338899099e-01, -5.5196108861576e-01],
            [3.8994443235820e00, 2.3891363245285e-01, -5.1748838848812e-01],
            [3.9372911834345e00, 2.8206060569333e-01, -5.6393504009155e-01],
            [3.9416129812188e00, 2.0832172987319e-01, -4.6855586031792e-01],
            [3.9431612976694e00, 1.8327640423061e-01, -5.3764973913994e-01],
            [3.8619577233846e00, 1.4192189812407e00, -9.2587403626770e-01],
            [3.8507167163959e00, 1.4238788373222e00, -9.3661710728291e-01],
            [3.8651039358730e00, 1.4201766685559e00, -9.2771824467570e-01],
            [3.8624692302920e00, 1.4273996853788e00, -9.3359662134515e-01],
            [3.8610467267790e00, 1.4182334490688e00, -9.3810025187748e-01],
            [3.8563372198902e00, 1.4215489092814e00, -9.3124557177530e-01],
            [3.8568487267976e00, 1.4297296134196e00, -9.3896815685275e-01],
            [3.8583881624179e00, 1.4255816848941e00, -9.4291768367439e-01],
            [3.8594834323787e00, 1.4225065965966e00, -9.3308978018331e-01],
        ]
    )

    assert np.allclose(nblock_expected, all_solid_cells_archive.nodes)

    grid = all_solid_cells_archive.grid
    assert grid.number_of_cells == 4
    assert np.unique(grid.celltypes).size == 4
    assert np.all(grid.celltypes > 20)

    assert np.all(all_solid_cells_archive.quality > 0.0)


def test_read_complex_archive_linear(all_solid_cells_archive_linear):
    grid = all_solid_cells_archive_linear.grid
    assert np.all(grid.celltypes < 20)
    assert np.all(all_solid_cells_archive_linear.quality > 0.0)


@pytest.mark.parametrize("celltype", QUADRATIC_CELL_TYPES)
def test_write_quad_complex_archive(tmpdir, celltype, all_solid_cells_archive):
    grid = all_solid_cells_archive.grid
    mask = grid.celltypes == celltype
    assert mask.any()
    grid = grid.extract_cells(mask)

    try:
        tmp_archive_file = str(tmpdir.mkdir("tmpdir").join("tmp.cdb"))
    except:
        tmp_archive_file = "/tmp/nblock.cdb"

    pymapdl_reader.save_as_archive(tmp_archive_file, grid)
    new_archive = pymapdl_reader.Archive(tmp_archive_file)
    assert np.allclose(grid.cells, new_archive.grid.cells)
    assert np.allclose(grid.points, new_archive.grid.points)
    assert (new_archive.quality > 0.0).all()


@pytest.mark.parametrize("celltype", LINEAR_CELL_TYPES)
def test_write_lin_archive(tmpdir, celltype, all_solid_cells_archive_linear):
    linear_grid = all_solid_cells_archive_linear.grid

    mask = linear_grid.celltypes == celltype
    assert mask.any()
    linear_grid = linear_grid.extract_cells(mask)

    tmp_archive_file = str(tmpdir.mkdir("tmpdir").join("tmp.cdb"))

    pymapdl_reader.save_as_archive(tmp_archive_file, linear_grid)
    new_archive = pymapdl_reader.Archive(tmp_archive_file)
    assert new_archive.quality > 0
    assert np.allclose(linear_grid.celltypes, new_archive.grid.celltypes)


def test_write_component(tmpdir):
    items = np.array([1, 20, 50, 51, 52, 53])
    filename = str(tmpdir.mkdir("tmpdir").join("tmp.cdb"))

    comp_name = "TEST"
    pymapdl_reader.write_cmblock(filename, items, comp_name, "node")
    archive = pymapdl_reader.Archive(filename)
    assert np.allclose(archive.node_components[comp_name], items)


def test_write_component_edge_case(tmpdir):
    items = np.arange(2, 34, step=2)
    filename = str(tmpdir.mkdir("tmpdir").join("tmp.cdb"))

    comp_name = "TEST"
    pymapdl_reader.write_cmblock(filename, items, comp_name, "node")
    archive = pymapdl_reader.Archive(filename)
    assert np.allclose(archive.node_components[comp_name], items)


def test_read_parm():
    filename = os.path.join(TESTFILES_PATH, "parm.cdb")
    archive = pymapdl_reader.Archive(filename)
    with pytest.raises(AttributeError):
        archive.parameters

    archive = pymapdl_reader.Archive(filename, read_parameters=True)
    assert len(archive.parameters) == 2
    for parm in archive.parameters:
        assert isinstance(archive.parameters[parm], np.ndarray)


def test_read_wb_nblock():
    expected = np.array(
        [
            [9.89367578e-02, -8.07092192e-04, 8.53764953e00],
            [9.65803244e-02, 2.00906704e-02, 8.53744951e00],
            [9.19243555e-02, 3.98781615e-02, 8.53723652e00],
        ]
    )
    filename = os.path.join(TESTFILES_PATH, "workbench_193.cdb")
    archive = pymapdl_reader.Archive(filename)
    assert np.allclose(archive.nodes, expected)
    assert np.allclose(archive.node_angles, 0)


def test_read_hypermesh():
    expected = np.array(
        [
            [-6.01203, 2.98129, 2.38556],
            [-3.03231, 2.98067, 2.38309],
            [-0.03485, 2.98004, 2.3805],
            [2.98794, 2.97941, 2.37773],
            [5.98956, 2.97878, 2.37488],
            [5.98956, 5.97878, 2.37488],
        ]
    )

    filename = os.path.join(TESTFILES_PATH, "hypermesh.cdb")
    archive = pymapdl_reader.Archive(filename, verbose=True)
    assert np.allclose(archive.nodes[:6], expected)


@pytest.mark.parametrize("angles", [True, False])
def test_cython_write_nblock(hex_archive, tmpdir, angles):
    nblock_filename = str(tmpdir.mkdir("tmpdir").join("nblock.inp"))

    if angles:
        _archive.py_write_nblock(
            nblock_filename,
            hex_archive.nnum,
            hex_archive.nnum[-1],
            hex_archive.nodes,
            hex_archive.node_angles,
        )
    else:
        _archive.py_write_nblock(
            nblock_filename,
            hex_archive.nnum,
            hex_archive.nnum[-1],
            hex_archive.nodes,
            np.empty((0, 0)),
        )

    tmp_archive = pymapdl_reader.Archive(nblock_filename)
    assert np.allclose(hex_archive.nnum, tmp_archive.nnum)
    assert np.allclose(hex_archive.nodes, tmp_archive.nodes)
    if angles:
        assert np.allclose(hex_archive.node_angles, tmp_archive.node_angles)


@pytest.mark.parametrize("has_angles", [True, False])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_write_nblock(hex_archive, tmpdir, dtype, has_angles):
    nblock_filename = str(tmpdir.mkdir("tmpdir").join("nblock.inp"))

    nodes = hex_archive.nodes.astype(dtype)
    if has_angles:
        angles = hex_archive.node_angles
    else:
        angles = None
    archive.write_nblock(nblock_filename, hex_archive.nnum, nodes, angles, mode="w")

    tmp_archive = pymapdl_reader.Archive(nblock_filename)
    assert np.allclose(hex_archive.nnum, tmp_archive.nnum)
    assert np.allclose(hex_archive.nodes, tmp_archive.nodes)
    if has_angles:
        assert np.allclose(hex_archive.node_angles, tmp_archive.node_angles)


def test_cython_write_eblock(hex_archive, tmpdir):
    filename = str(tmpdir.mkdir("tmpdir").join("eblock.inp"))

    etype = np.ones(hex_archive.n_elem, np.int32)
    typenum = hex_archive.etype
    elem_nnodes = np.empty(etype.size, np.int32)
    elem_nnodes[typenum == 181] = 4
    elem_nnodes[typenum == 185] = 8
    elem_nnodes[typenum == 186] = 20
    elem_nnodes[typenum == 187] = 10
    nodenum = hex_archive.nnum

    cells, offset = pymapdl_reader.misc.vtk_cell_info(
        hex_archive.grid,
        force_int64=False,
        shift_offset=False,
    )
    _archive.py_write_eblock(
        filename,
        hex_archive.enum,
        etype,
        hex_archive.material_type,
        np.ones(hex_archive.n_elem, np.int32),
        elem_nnodes,
        cells.astype(np.int32, copy=False),
        offset.astype(np.int32, copy=False),
        hex_archive.grid.celltypes,
        typenum,
        nodenum,
    )


def test_rlblock_prior_to_nblock():
    # test edge case where RLBLOCK is immediately prior to the NBLOCK
    filename = os.path.join(TESTFILES_PATH, "ErnoRadiation.cdb")
    archive = pymapdl_reader.Archive(filename)
    assert archive.n_node == 65
    assert archive.n_elem == 36


class TestPathlibFilename:
    def test_pathlib_filename_property(self, pathlib_archive):
        assert isinstance(pathlib_archive.pathlib_filename, pathlib.Path)

    def test_filename_property_is_string(self, pathlib_archive):
        filename = TESTFILES_PATH_PATHLIB / "ErnoRadiation.cdb"
        a = pymapdl_reader.Archive(filename)
        assert isinstance(a.filename, str)

    def test_filename_setter_pathlib(self, pathlib_archive):
        with pytest.raises(AttributeError):
            pathlib_archive.filename = pathlib.Path("dummy2")

    def test_filename_setter_string(self, pathlib_archive):
        with pytest.raises(AttributeError):
            pathlib_archive.filename = "dummy2"
