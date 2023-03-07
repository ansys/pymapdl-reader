"""
Sample result file generated with

import os
from pyansys import examples
import pyansys

os.environ['I_MPI_SHM_LMT'] = 'shm'
mapdl = pyansys.launch_mapdl(override=True)

mapdl.cdread('db', examples.hexarchivefile)
mapdl.esel('s', 'ELEM', vmin=5, vmax=20)
mapdl.cm('ELEM_COMP', 'ELEM')
mapdl.nsel('s', 'NODE', vmin=5, vmax=20)
mapdl.cm('NODE_COMP', 'NODE')

# boundary conditions
mapdl.allsel()

# dummy steel properties
mapdl.prep7()
mapdl.mp('EX', 1, 200E9)  # Elastic moduli in Pa (kg/(m*s**2))
mapdl.mp('DENS', 1, 7800)  # Density in kg/m3
mapdl.mp('NUXY', 1, 0.3)  # Poissons Ratio
mapdl.emodif('ALL', 'MAT', 1)

# fix one end of the beam
mapdl.nsel('S', 'LOC', 'Z')
mapdl.d('all', 'all')
mapdl.allsel()

mapdl.modal_analysis(nmode=1)

"""
import os
import pathlib
import platform
from shutil import copy

import numpy as np
import pytest
from pyvista.plotting import system_supports_plotting
from pyvista.plotting.renderer import CameraPosition

from ansys.mapdl import reader as pymapdl_reader
from ansys.mapdl.reader import examples
from ansys.mapdl.reader.examples.downloads import _download_and_read
from ansys.mapdl.reader.rst import Result

try:
    vm33 = examples.download_verification_result(33)
except:
    vm33 = None

try:
    vm240 = examples.download_verification_result(240)
except:
    vm240 = None

try:
    vm240_sparse = _download_and_read("vm240_sparse.rst")
except:
    vm240_sparse = None


try:
    pontoon = examples.download_pontoon()
except:
    pontoon = None

IS_MAC = platform.system() == "Darwin"
skip_plotting = pytest.mark.skipif(
    not system_supports_plotting() or IS_MAC, reason="Requires active X Server"
)

test_path = os.path.dirname(os.path.abspath(__file__))
testfiles_path = os.path.join(test_path, "testfiles")

is16_filename = os.path.join(testfiles_path, "is16.rst")
is16_known_result = os.path.join(testfiles_path, "is16.npz")
if os.path.isfile(is16_filename):
    is16 = pymapdl_reader.read_binary(is16_filename)
else:
    is16 = None

temperature_rst = os.path.join(testfiles_path, "temp_v13.rst")
temperature_known_result = os.path.join(testfiles_path, "temp_v13.npz")


@pytest.fixture(scope="module")
def pathlib_result():
    temperature_rst_pathlib = pathlib.Path(temperature_rst)
    return Result(temperature_rst_pathlib)


@pytest.fixture(scope="module")
def hex_pipe_corner():
    filename = os.path.join(testfiles_path, "rst", "cyc_stress.rst")
    return pymapdl_reader.read_binary(filename)


@pytest.fixture(scope="module")
def hex_rst():
    filename = os.path.join(testfiles_path, "hex_201.rst")
    return pymapdl_reader.read_binary(filename)


@pytest.fixture(scope="module")
def volume_rst():
    rst_file = os.path.join(testfiles_path, "vol_test.rst")
    return pymapdl_reader.read_binary(rst_file)


@pytest.fixture(scope="module")
def materials_281_rst():
    rst_file = os.path.join(testfiles_path, "materials", "file.rst")
    return pymapdl_reader.read_binary(rst_file)


@pytest.fixture(scope="module")
def materials_stress_lim_rst():
    rst_file = os.path.join(testfiles_path, "materials", "stress_lim.rst")
    return pymapdl_reader.read_binary(rst_file)


def test_map_flag_section_data():
    # basic 4 element SHELL181 result files
    shell181_2020r2 = os.path.join(testfiles_path, "shell181_2020R2.rst")
    shell181_2021r1 = os.path.join(testfiles_path, "shell181_2021R1.rst")

    rst_2020r2 = pymapdl_reader.read_binary(shell181_2020r2)
    rst_2021r1 = pymapdl_reader.read_binary(shell181_2021r1)

    for key in rst_2020r2.section_data:
        assert np.allclose(rst_2020r2.section_data[key], rst_2021r1.section_data[key])


def test_overwrite(tmpdir):
    tmp_path = str(tmpdir.mkdir("tmpdir"))
    rst = pymapdl_reader.read_binary(copy(examples.rstfile, tmp_path))

    # get the data
    solution_type = "EEL"
    enum, esol_old, _ = rst.element_solution_data(0, solution_type)

    index = 10
    element_id = enum[index]
    old_record = esol_old[index]
    ovr_record = np.random.random(old_record.size)

    rst.overwrite_element_solution_record(ovr_record, 0, solution_type, element_id)

    # verify data has been written
    new_record = rst.element_solution_data(0, solution_type)[1][index]
    assert not np.allclose(new_record, old_record), "nothing overwritten"
    assert np.allclose(ovr_record, new_record)


def test_overwrite_dict(tmpdir):
    tmp_path = str(tmpdir.mkdir("tmpdir"))
    rst = pymapdl_reader.read_binary(copy(examples.rstfile, tmp_path))

    # get the data
    solution_type = "EEL"
    enum, esol_old, _ = rst.element_solution_data(0, solution_type)

    indices = (10, 20)
    old_records = [esol_old[index] for index in indices]

    element_data = {
        enum[indices[0]]: np.random.random(old_records[0].size),
        enum[indices[1]]: np.random.random(old_records[1].size),
    }

    rst.overwrite_element_solution_records(element_data, 0, solution_type)

    # verify data has been written
    new_record = rst.element_solution_data(0, solution_type)[1]

    for i, index in enumerate(indices):
        assert not np.allclose(old_records[i], new_record[index]), "nothing overwritten"
        assert np.allclose(element_data[enum[indices[i]]], new_record[index])


@pytest.mark.skipif(vm33 is None, reason="Requires example files")
def test_write_tables(tmpdir):
    filename = str(tmpdir.mkdir("tmpdir").join("vm33.txt"))
    vm33.write_tables(filename)
    assert os.path.isfile(filename)


@pytest.mark.skipif(vm33 is None, reason="Requires example files")
def test_nodal_displacement():
    nnum, result = vm33.nodal_solution(0)
    nnum_, result_ = vm33.nodal_displacement(0)
    assert np.allclose(nnum, nnum_)
    assert np.allclose(result, result_)


def test_read_volume(volume_rst):
    enum, edata, enode = volume_rst.element_solution_data(0, datatype="ENG")
    edata = np.asarray(edata)
    volume = edata[:, 0]

    enum_vtk = np.sort(volume_rst.grid.cell_data["ansys_elem_num"])
    assert np.allclose(enum, enum_vtk)
    assert np.allclose(volume, 291895460.0)


@pytest.mark.skipif(vm33 is None, reason="Requires example files")
def nodal_displacement():
    nnum, disp = vm33.nodal_displacement(0)
    assert isinstance(nnum, np.ndarray)
    assert isinstance(disp, np.ndarray)


@pytest.mark.skipif(vm33 is None, reason="Requires example files")
def test_nodal_thermal_strain():
    _, tstrain = vm33.nodal_thermal_strain(0)
    assert np.any(tstrain)
    assert tstrain.shape == (vm33.grid.n_points, 8)


@skip_plotting
@pytest.mark.skipif(vm33 is None, reason="Requires example files")
def test_plot_nodal_thermal_strain():
    vm33.plot_nodal_thermal_strain(0, "X", off_screen=True)


@skip_plotting
@pytest.mark.skipif(vm33 is None, reason="Requires example files")
def test_plot_nodal_thermal_strain():
    vm33._animate_time_solution("ENS", off_screen=True)


@pytest.mark.skipif(pontoon is None, reason="Requires example files")
def test_nodal_elastic_strain():
    _, estrain = pontoon.nodal_elastic_strain(0)
    assert np.any(estrain)
    assert estrain.shape == (pontoon.grid.n_points, 7)


@skip_plotting
@pytest.mark.skipif(pontoon is None, reason="Requires example files")
def test_plot_nodal_elastic_strain():
    pontoon.plot_nodal_elastic_strain(0, "X", off_screen=True)


@skip_plotting
@pytest.mark.skipif(pontoon is None, reason="Requires example files")
def test_plot_pontoon():
    pontoon.plot(off_screen=True)


@skip_plotting
@pytest.mark.skipif(pontoon is None, reason="Requires example files")
def test_plot_pontoon_nodal_displacement():
    pontoon.plot_nodal_solution(
        0, show_displacement=True, overlay_wireframe=True, off_screen=True
    )


@skip_plotting
@pytest.mark.skipif(pontoon is None, reason="Requires example files")
def test_plot_pontoon_nodal_displacement_with_text_color():
    pontoon.plot_nodal_solution(
        0,
        show_displacement=True,
        overlay_wireframe=True,
        off_screen=True,
        text_color="k",
    )


@skip_plotting
@pytest.mark.skipif(pontoon is None, reason="Requires example files")
def test_plot_pontoon_nodal_displacement():
    pontoon.plot_nodal_displacement(
        0, show_displacement=True, overlay_wireframe=True, off_screen=True
    )


@pytest.mark.skipif(pontoon is None, reason="Requires example files")
def test_print_pontoon_components():
    assert isinstance(pontoon.node_components, dict)


@pytest.mark.skipif(pontoon is None, reason="Requires example files")
def test_repr():
    assert "Title" in str(pontoon)


@pytest.mark.skipif(vm33 is None, reason="Requires example files")
def test_available_results():
    for etype in vm33.available_results:
        if etype in ["NSL", "VSL", "ASL", "RF"]:  # skip nodal records
            continue
        _, edata, _ = vm33.element_solution_data(0, etype)
        assert edata[0] is not None


@pytest.mark.skipif(vm33 is None, reason="Requires example files")
def test_solution_info():
    info = vm33.solution_info(0)
    assert "omega_a_x" in info


@pytest.mark.skipif(
    vm240 is None or vm240_sparse is None, reason="Requires example files"
)
def test_sparse_nodal_solution():
    nnum, stress = vm240.nodal_stress(0)
    sparse_nnum, sparse_stress = vm240_sparse.nodal_stress(0)
    assert np.allclose(sparse_stress, stress, equal_nan=True)
    assert np.allclose(nnum, sparse_nnum)


@pytest.mark.skipif(is16 is None, reason="Requires example files")
def test_is16():
    npz_rst = np.load(is16_known_result)
    nnum, data = is16.nodal_solution(0)
    # remove NAN data
    mask = ~np.isnan(data[:, 0])
    nnum = nnum[mask]
    data = data[mask]

    assert np.allclose(data, npz_rst["data"], atol=1e-6)
    assert np.allclose(nnum, npz_rst["nnum"])


@pytest.mark.skipif(
    not os.path.isfile(temperature_rst), reason="Requires example files"
)
def test_read_temperature():
    temp_rst = pymapdl_reader.read_binary(temperature_rst)
    nnum, temp = temp_rst.nodal_temperature(0)

    npz_rst = np.load(temperature_known_result)
    assert np.allclose(nnum, npz_rst["nnum"])
    assert np.allclose(temp, npz_rst["temp"])


@skip_plotting
@pytest.mark.skipif(
    not os.path.isfile(temperature_rst), reason="Requires example files"
)
def test_plot_nodal_temperature():
    temp_rst = pymapdl_reader.read_binary(temperature_rst)
    temp_rst.plot_nodal_temperature(0, off_screen=True)


class TestPathlibFilename:
    @skip_plotting
    @pytest.mark.skipif(
        not os.path.isfile(temperature_rst), reason="Requires example files"
    )
    def test_pathlib_filename_property(self, pathlib_result):
        assert isinstance(pathlib_result.pathlib_filename, pathlib.Path)

    @skip_plotting
    @pytest.mark.skipif(
        not os.path.isfile(temperature_rst), reason="Requires example files"
    )
    def test_filename_property_is_string(self, pathlib_result):
        assert isinstance(pathlib_result.filename, str)

    @skip_plotting
    @pytest.mark.skipif(
        not os.path.isfile(temperature_rst), reason="Requires example files"
    )
    def test_filename_setter_pathlib(self, pathlib_result):
        with pytest.raises(AttributeError):
            pathlib_result.filename = pathlib.Path("dummy2")

    @skip_plotting
    @pytest.mark.skipif(
        not os.path.isfile(temperature_rst), reason="Requires example files"
    )
    def test_filename_setter_string(self, pathlib_result):
        with pytest.raises(AttributeError):
            pathlib_result.filename = "dummy2"


def test_rst_node_components(hex_rst):
    assert "ELEM_COMP" not in hex_rst.node_components
    np.allclose(hex_rst.node_components["NODE_COMP"].nonzero()[0], np.arange(4, 20))


def test_rst_node_components(hex_rst):
    assert "NODE_COMP" not in hex_rst.element_components
    np.allclose(hex_rst.element_components["ELEM_COMP"].nonzero()[0], np.arange(4, 20))


def test_rst_beam4_shell63():
    filename = os.path.join(testfiles_path, "shell63_beam4.rst")

    # check geometry can load
    rst = pymapdl_reader.read_binary(filename)
    assert isinstance(rst.grid.points, np.ndarray)
    assert isinstance(rst.grid.cells, np.ndarray)
    assert (rst.grid.cells >= 0).all()
    assert (rst.grid.cells <= 50).all()

    # not a great test, but verifies results load
    assert np.any(rst.nodal_displacement(0)[1] > 0)


def test_cyl_stress(hex_pipe_corner):
    # ANSYS results generated with
    # RSYS, 0
    # PRNSOL, S
    # RSYS, 1
    # PRNSOL, S
    _, my_stress = hex_pipe_corner.cylindrical_nodal_stress(0)
    ans_stress = np.load(os.path.join(testfiles_path, "rst", "cyc_stress.npy"))

    assert np.allclose(my_stress[-114:], ans_stress, atol=1e-7)


@skip_plotting
def test_plot_cyl_stress(hex_pipe_corner):
    with pytest.raises(ValueError):
        cpos = hex_pipe_corner.plot_cylindrical_nodal_stress(0, off_screen=True)
    with pytest.raises(ValueError):
        cpos = hex_pipe_corner.plot_cylindrical_nodal_stress(
            0, comp="X", off_screen=True
        )
    cpos = hex_pipe_corner.plot_cylindrical_nodal_stress(0, comp="R", off_screen=True)
    if cpos is not None:
        assert isinstance(cpos, CameraPosition)


def test_reaction_forces(volume_rst):
    known_result = os.path.join(testfiles_path, "nodal_reaction.npy")
    nnum_known, fz_known = np.load(known_result).T
    nnum_known = nnum_known.astype(np.int32)

    rforces, nnum, dof = volume_rst.nodal_reaction_forces(0)
    assert (np.diff(nnum) >= 0).all()
    assert (np.in1d(dof, [1, 2, 3])).all()
    assert rforces.dtype == np.float64

    fz = rforces[dof == 3]
    # loose tolerance due to table printed from MAPDL
    assert np.allclose(fz_known, fz, rtol=1e-4)
    assert np.allclose(nnum_known, nnum[dof == 1])


@pytest.mark.parametrize("nnum_of_interest", [range(11, 50), "all"])
def test_nnum_of_interest(nnum_of_interest):
    rst = pymapdl_reader.read_binary(examples.rstfile)
    if nnum_of_interest == "all":
        nnum_of_interest = rst.mesh.nnum

    nnum_sel, data_sel = rst._nodal_result(0, "ENS", nnum_of_interest=nnum_of_interest)
    nnum, data = rst._nodal_result(0, "ENS")

    mask = np.in1d(nnum, nnum_of_interest)
    assert np.allclose(nnum[mask], nnum_sel)
    assert np.allclose(data[mask], data_sel, equal_nan=True)


@pytest.mark.parametrize("nodes", [range(11, 50), "NCOMP2", ("NCOMP2", "NODE_COMP")])
def test_nodes_subselection(hex_rst, nodes):
    nnum_sel, data_sel = hex_rst.nodal_solution(0, nodes=nodes)
    nnum, data = hex_rst.nodal_solution(0)

    grid_nnum = hex_rst.grid.point_data["ansys_node_num"]
    if isinstance(nodes, str):
        nnum_of_interest = grid_nnum[hex_rst.grid.point_data[nodes].view(bool)]
    elif isinstance(nodes, tuple):
        mask = np.logical_or(
            hex_rst.grid.point_data[nodes[0]].view(bool),
            hex_rst.grid.point_data[nodes[1]].view(bool),
        )
        nnum_of_interest = grid_nnum[mask]
    else:
        nnum_of_interest = nodes

    mask = ~np.isnan(data_sel[:, 0])

    assert np.allclose(nnum[mask], nnum_sel[mask])
    assert np.allclose(data[mask], data_sel[mask], equal_nan=True)


def test_materials(materials_281_rst):
    # known material results from MPLIST,ALL
    mat = {}
    mat[1] = {
        "EX": 200000.0,
        "NUXY": 0.3000000,
        "ALPX": 0.1200000e-04,
        "DENS": 0.7850000e-08,
        "KXX": 60.50000,
        "RSVX": 0.1700000e-03,
        "C": 0.4340000e09,
        "MURX": 10000.00,
    }

    mat[2] = {
        "EX": 47869.10,
        "EY": 12099.50,
        "EZ": 12099.50,
        "NUXY": 0.7621539e-01,
        "NUYZ": 0.4249800,
        "NUXZ": 0.7621539e-01,
        "GXY": 4118.700,
        "GYZ": 4245.500,
        "GXZ": 4118.700,
        "DENS": 0.1947300e-08,
        "PRXY": 0.3015300,
        "PRYZ": 0.4249800,
        "PRXZ": 0.3015300,
    }

    mat[3] = {
        "EX": 38479.10,
        "EY": 10456.20,
        "EZ": 10456.20,
        "NUXY": 0.8242602e-01,
        "NUYZ": 0.4435900,
        "NUXZ": 0.8242602e-01,
        "GXY": 3548.900,
        "GYZ": 3621.600,
        "GXZ": 3548.900,
        "DENS": 0.1877000e-08,
        "PRXY": 0.3033300,
        "PRYZ": 0.4435900,
        "PRXZ": 0.3033300,
    }

    for mat_type, known_material in mat.items():
        ans_mat = materials_281_rst.materials[mat_type]
        for key, value in known_material.items():
            assert pytest.approx(ans_mat[key]) == known_material[key]


def test_materials_str_lim(materials_stress_lim_rst):
    # output from:
    # from ansys.mapdl.core import launch_mapdl
    # mapdl = launch_mapdl()
    # mapdl.post1()
    # mapdl.set(1, 1)
    # print(mapdl.tblist(mat=1))
    #
    # validate material stress limits can be extracted
    #                   FC LIMIT (FCLI) Table For Material     1
    #               Stress Strength Limits
    #
    #                     1
    # Temps     0.0000000e+00
    # XTEN     5.0000000e-01
    # XCMP     -2.5000000e-01
    # YTEN     5.0000000e-01
    # YCMP     -2.5000000e-01
    # ZTEN     5.0000000e+00
    # ZCMP     -5.0000000e+00
    #   XY     5.0000000e-01
    #   YZ     3.0000000e+00
    #   XZ     3.0000000e+00
    # XYCP     0.0000000e+00
    # YZCP     0.0000000e+00
    # XZCP     0.0000000e+00
    # XZIT     0.0000000e+00
    # XZIC     0.0000000e+00
    # YZIT     0.0000000e+00
    # YZIC     0.0000000e+00
    # GI/GII   0.0000000e+00
    # ETA_L    0.0000000e+00
    # ETA_T    0.0000000e+00
    # ALPHA_0  0.0000000e+00
    #
    #               FC LIMIT (FCLI) Table For Material     2
    #               Stress Strength Limits
    #                     1
    # Temps     0.0000000e+00
    # XTEN     2.0000000e+00
    # XCMP     -2.0000000e+00
    # YTEN     2.0000000e+00
    # YCMP     -2.0000000e+00
    # ZTEN     3.0000000e+00
    # ZCMP     -4.0000000e+00
    #   XY     1.0000000e+00
    #   YZ     2.0000000e+00
    #   XZ     2.0000000e+00
    # XYCP     0.0000000e+00
    # YZCP     0.0000000e+00
    # XZCP     0.0000000e+00
    # XZIT     0.0000000e+00
    # XZIC     0.0000000e+00
    # YZIT     0.0000000e+00
    # YZIC     0.0000000e+00
    # GI/GII   0.0000000e+00
    # ETA_L    0.0000000e+00
    # ETA_T    0.0000000e+00
    # ALPHA_0  0.0000000e+00

    mat = {}
    mat[1] = {
        "XTEN": 5.0000000e-01,
        "XCMP": -2.5000000e-01,
        "YTEN": 5.0000000e-01,
        "YCMP": -2.5000000e-01,
        "ZTEN": 5.0000000e00,
        "ZCMP": -5.0000000e00,
        "XY": 5.0000000e-01,
        "YZ": 3.0000000e00,
        "XZ": 3.0000000e00,
        "XYCP": 0.0000000e00,
        "YZCP": 0.0000000e00,
        "XZCP": 0.0000000e00,
        "XZIT": 0.0000000e00,
        "XZIC": 0.0000000e00,
        "YZIT": 0.0000000e00,
        "YZIC": 0.0000000e00,
    }

    mat[2] = {
        "XTEN": 2.0000000e00,
        "XCMP": -2.0000000e00,
        "YTEN": 2.0000000e00,
        "YCMP": -2.0000000e00,
        "ZTEN": 3.0000000e00,
        "ZCMP": -4.0000000e00,
        "XY": 1.0000000e00,
        "YZ": 2.0000000e00,
        "XZ": 2.0000000e00,
        "XYCP": 0.0000000e00,
        "YZCP": 0.0000000e00,
        "XZCP": 0.0000000e00,
        "XZIT": 0.0000000e00,
        "XZIC": 0.0000000e00,
        "YZIT": 0.0000000e00,
        "YZIC": 0.0000000e00,
    }

    ans_mat = materials_stress_lim_rst.materials[1]["stress_failure_criteria"]

    for mat_type, known_stress_lim in mat.items():
        ans_mat = materials_stress_lim_rst.materials[mat_type][
            "stress_failure_criteria"
        ]
        for key, value in known_stress_lim.items():
            assert pytest.approx(ans_mat[key]) == known_stress_lim[key]


def test_materials_v150():
    """Validate on older result files"""
    rst = pymapdl_reader.read_binary(examples.rstfile)
    mat = {1: {"EX": 16900000.0, "NUXY": 0.31, "DENS": 0.00041408}}

    for mat_type, known_material in mat.items():
        ans_mat = rst.materials[mat_type]
        for key, value in known_material.items():
            assert pytest.approx(ans_mat[key]) == known_material[key]


def test_temperature_dependent_properties():
    path_rth = os.path.join(testfiles_path, "temp_dependent_results.rth")
    rst = pymapdl_reader.read_binary(path_rth)
    mats = rst.materials
    kxx = mats[1]["KXX"]
    expected_value = np.array(
        [[-100.0, 0.0, 100.0, 200.0], [114.0, 144.0, 165.0, 175.0]]
    )

    assert kxx is not None
    assert isinstance(kxx, np.ndarray)
    assert kxx.size > 0
    assert np.allclose(kxx, expected_value)
