# Copyright (C) 2021 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import platform
import shutil

import numpy as np
import pytest
import pyvista as pv
from pyvista.plotting import system_supports_plotting
from pyvista.plotting.renderer import CameraPosition

from ansys.mapdl import reader as pymapdl_reader
from ansys.mapdl.reader import examples
from ansys.mapdl.reader._rst_keys import element_index_table_info

try:
    from ansys.mapdl.core import _HAS_ANSYS
except:
    _HAS_ANSYS = False

if os.environ.get("SKIP_ANSYS", "").upper() == "TRUE":
    _HAS_ANSYS = False


HAS_FFMPEG = True
try:
    import imageio_ffmpeg  # noqa: F401
except ImportError:
    HAS_FFMPEG = False

skip_no_ansys = pytest.mark.skipif(not _HAS_ANSYS, reason="Requires ANSYS installed")


test_path = os.path.dirname(os.path.abspath(__file__))
testfiles_path = os.path.join(test_path, "testfiles")

IS_MAC = platform.system() == "Darwin"
skip_plotting = pytest.mark.skipif(
    not system_supports_plotting() or IS_MAC, reason="Requires active X Server"
)

RSETS = list(zip(range(1, 9), [1] * 8))


@pytest.fixture(scope="module")
def result():
    return pymapdl_reader.read_binary(examples.rstfile)


@pytest.fixture(scope="module")
def static_canteliver_bc():
    filename = os.path.join(testfiles_path, "rst", "beam_static_bc.rst")
    return pymapdl_reader.read_binary(filename)


@pytest.fixture(scope="module")
def thermal_rst():
    filename = os.path.join(testfiles_path, "file.rth")
    return pymapdl_reader.read_binary(filename)


@pytest.fixture(scope="module")
def cyclic_modal(mapdl):
    # build the cyclic model
    mapdl.clear()
    mapdl.prep7()
    mapdl.shpp("off")
    mapdl.cdread("db", pymapdl_reader.examples.sector_archive_file)
    mapdl.prep7()
    mapdl.cyclic()

    # set material properties
    mapdl.mp("NUXY", 1, 0.31)
    mapdl.mp("DENS", 1, 4.1408e-04)
    mapdl.mp("EX", 1, 16900000)
    mapdl.emodif("ALL", "MAT", 1)

    # setup and solve
    mapdl.run("/SOLU")
    mapdl.antype(2, "new")
    mapdl.modopt("lanb", 1, 1, 100000)
    mapdl.eqslv("SPARSE")
    mapdl.lumpm(0)
    mapdl.pstres(0)
    mapdl.bcsoption("INCORE")
    mapdl.mxpand(elcalc="YES")
    mapdl.solve()
    mapdl.finish()

    # setup ansys for output without line breaks
    mapdl.post1()
    mapdl.header("OFF", "OFF", "OFF", "OFF", "OFF", "OFF")
    nsigfig = 10
    mapdl.format("", "E", nsigfig + 9, nsigfig)
    mapdl.page(1e9, "", -1, 240)


@pytest.fixture(scope="module")
def transient_thermal(mapdl):
    mapdl.finish()
    mapdl.clear()
    mapdl.prep7()

    # Material properties-- 1020 steel in imperial
    mapdl.units("BIN")  # U.S. Customary system using inches (in, lbf*s2/in, s, °F).
    mapdl.mp("EX", 1, 30023280.0)
    mapdl.mp("NUXY", 1, 0.290000000)
    mapdl.mp("ALPX", 1, 8.388888889e-06)
    mapdl.mp("DENS", 1, 7.346344000e-04)
    mapdl.mp("KXX", 1, 6.252196000e-04)
    mapdl.mp("C", 1, 38.6334760)

    # use a thermal element type
    mapdl.et(1, "SOLID70")

    # Geometry and Mesh
    mapdl.block(0, 5, 0, 1, 0, 1)
    mapdl.lesize("ALL", 1, layer1=1)
    mapdl.mshape(0, "3D")
    mapdl.mshkey(1)
    mapdl.vmesh(1)
    mapdl.run("/SOLU")
    mapdl.antype(4)  # transient analysis
    mapdl.trnopt("FULL")  # full transient analysis
    mapdl.kbc(0)  # ramp loads up and down

    # Time stepping
    end_time = 500
    mapdl.time(end_time)  # end time for load step
    mapdl.autots("ON")  # use automatic time stepping
    mapdl.deltim(10, 5, 100)  # substep size (seconds)

    # Create a table of convection times and coefficients and transfer it to MAPDL
    my_conv = np.array(
        [
            [0, 0.001],  # start time
            [120, 0.001],  # end of first "flat" zone
            [130, 0.005],  # ramps up in 10 seconds
            [300, 0.005],  # end of second "flat zone
            [400, 0.002],  # ramps down in 10 seconds
            [end_time, 0.002],
        ]
    )  # end of third "flat" zone
    mapdl.load_table("my_conv", my_conv, "TIME")

    # Create a table of bulk temperatures for a given time and transfer to MAPDL
    my_bulk = np.array(
        [
            [0, 100],  # start time
            [120, 100],  # end of first "flat" zone
            [200, 300],  # ramps up in 380 seconds
            [300, 300],  # hold temperature for 200 seconds
            [400, 75],  # temperature ramps down for 200 seconds
            [end_time, 75],
        ]
    )  # end of second "flat" zone
    mapdl.load_table("my_bulk", my_bulk, "TIME")

    # Force transient solve to include the times within the conv and bulk arrays
    # my_tres = np.unique(np.vstack((my_bulk[:, 0], my_conv[:, 0])))[0]  # same as
    mapdl.parameters["my_tsres"] = [120, 130, 300, 400, end_time]
    mapdl.tsres("%my_tsres%")

    mapdl.outres("ERASE")
    mapdl.outres("ALL", "ALL")

    mapdl.eqslv("SPARSE")  # use sparse solver
    mapdl.tunif(75)  # force uniform starting temperature (otherwise zero)

    # apply the convective load (convection coefficient plus bulk temperature)
    # use "%" around table array names
    mapdl.sfa(6, 1, "CONV", "%my_conv%", " %my_bulk%")

    # solve
    mapdl.solve()
    mapdl.post1()
    mapdl.set(1, 1)


@pytest.mark.parametrize("rset", RSETS)
@skip_no_ansys
def test_prnsol_u(mapdl, cyclic_modal, rset):
    mapdl.set(*rset)
    # verify cyclic displacements
    array = mapdl.prnsol("u").to_array()
    ansys_nnum = array[:, 0].astype(np.int_)
    ansys_disp = array[:, 1:-1]

    nnum, disp = mapdl.result.nodal_solution(rset)

    # cyclic model will only output the master sector
    ansys_nnum = ansys_nnum[: nnum.size]
    ansys_disp = ansys_disp[: nnum.size]

    assert np.allclose(ansys_nnum, nnum)
    assert np.allclose(ansys_disp, disp)


@pytest.mark.parametrize("rset", RSETS)
@skip_no_ansys
def test_presol_s(mapdl, cyclic_modal, rset):
    mapdl.set(*rset)

    # verify element stress
    _, element_stress, enode = mapdl.result.element_stress(rset)
    element_stress = np.vstack(element_stress)
    enode = np.hstack(enode)

    # parse ansys result
    array = mapdl.presol("S").to_array()
    ansys_enode = array[:, 0].astype(np.int_)
    ansys_element_stress = array[:, 1:]

    arr_sz = element_stress.shape[0]
    assert np.allclose(element_stress, ansys_element_stress[:arr_sz])
    assert np.allclose(enode, ansys_enode[:arr_sz])


@pytest.mark.parametrize("rset", RSETS)
@skip_no_ansys
def test_prnsol_s(mapdl, cyclic_modal, rset):
    mapdl.set(*rset)

    # verify cyclic displacements
    array = mapdl.prnsol("s").to_array()
    ansys_nnum = array[:, 0].astype(np.int_)
    ansys_stress = array[:, 1:]

    nnum, stress = mapdl.result.nodal_stress(rset)

    # v150 includes nodes in the geometry that aren't in the result
    mask = np.isin(nnum, ansys_nnum)
    nnum = nnum[mask]
    stress = stress[mask]

    arr_sz = nnum.shape[0]
    assert np.allclose(nnum, ansys_nnum[:arr_sz])
    assert np.allclose(stress, ansys_stress[:arr_sz])


@pytest.mark.parametrize("rset", RSETS)
@skip_no_ansys
def test_prnsol_prin(mapdl, cyclic_modal, rset):
    mapdl.set(*rset)

    # verify principal stress
    array = mapdl.prnsol("prin").to_array()
    ansys_nnum = array[:, 0].astype(np.int_)
    ansys_stress = array[:, 1:]

    nnum, stress = mapdl.result.principal_nodal_stress(rset)

    # v150 includes nodes in the geometry that aren't in the result
    mask = np.isin(nnum, ansys_nnum)
    nnum = nnum[mask]
    stress = stress[mask]

    arr_sz = nnum.shape[0]
    assert np.allclose(nnum, ansys_nnum[:arr_sz])
    assert np.allclose(stress, ansys_stress[:arr_sz], atol=1e-5, rtol=1e-3)


def test_loadresult(result):
    # check result is loaded
    assert result.nsets
    assert result.mesh.nnum.size

    # check geometry is generated
    grid = result.grid
    assert grid.points.size
    assert grid.cells.size
    assert "ansys_node_num" in grid.point_data

    # check results can be loaded
    nnum, disp = result.nodal_solution(0)
    assert nnum.size
    assert disp.size

    nnum, disp = result.nodal_solution(0)
    assert nnum.size
    assert disp.size

    nnum, disp = result.principal_nodal_stress(0)
    assert nnum.size
    assert disp.size

    nnum, disp = result.nodal_stress(0)
    assert nnum.size
    assert disp.size

    enum, element_stress, enode = result.element_stress(0)
    assert element_stress[0].size
    assert enum.size
    assert enode[0].size

    enum, element_stress, enode = result.element_stress(0, principal=True)
    assert element_stress[0].size
    assert enum.size
    assert enode[0].size


def test_boundary_conditions(static_canteliver_bc):
    nnum, dof, bc = static_canteliver_bc.nodal_boundary_conditions(0)
    assert np.allclose(nnum, np.array([3, 3, 3, 25, 26, 27]))
    assert np.allclose(dof, np.array([1, 2, 3, 1, 2, 3]))
    assert np.allclose(bc, np.array([0.0, 0.0, 0.0, 0.001, 0.0011, 0.0012]))


def test_force(static_canteliver_bc):
    nnum, dof, force = static_canteliver_bc.nodal_input_force(0)
    assert np.allclose(nnum, np.array([71, 52, 127]))
    assert np.allclose(dof, np.array([2, 1, 3]))
    assert np.allclose(force, np.array([30.0, 20.0, 40.0]))


def test_dof(result):
    assert result.result_dof(0) == ["UX", "UY", "UZ"]


result_types = ["ENS", "EPT", "ETH", "EEL", "ENG"]  # 'ENF']


@pytest.mark.parametrize("result_type", result_types)
def test_save_as_vtk(tmpdir, result, result_type):
    filename = str(tmpdir.mkdir("tmpdir").join("tmp.vtk"))
    result.save_as_vtk(filename, result_types=[result_type])

    grid = pv.UnstructuredGrid(filename)
    for i in range(result.nsets):
        key = "Nodal Solution %d" % i
        assert key in grid.point_data
        arr = grid.point_data[key]
        assert np.allclose(arr, result.nodal_solution(i)[1], atol=1e-5)

        key = "%s %d" % (element_index_table_info[result_type], i)
        assert key in grid.point_data
        arr = grid.point_data[key]
        _, rst_arr = result._nodal_result(i, result_type)
        if rst_arr.shape[1] == 1:
            rst_arr = rst_arr.ravel()
        assert np.allclose(arr, rst_arr, atol=1e-5, equal_nan=True)


@skip_plotting
def test_plot_component():
    """
    # create example file for component plotting
    ansys = launch_mapdl()
    ansys.Cdread('db', examples.hexarchivefile)
    # ansys.open_gui()
    ansys.Esel('S', 'ELEM', vmin=1, vmax=20)
    ansys.Nsle('S', 1)
    ansys.Cm('MY_COMPONENT', 'NODE')

    ansys.Nsel('S', 'NODE', '', 1, 40)
    ansys.Cm('MY_OTHER_COMPONENT', 'NODE')

    ansys.Allsel()

    # Aluminum properties (or something)
    ansys.Prep7()
    ansys.Mp('ex', 1, 10e6)
    ansys.Mp('nuxy', 1, 0.3)
    ansys.Mp('dens', 1, 0.1/386.1)
    ansys.Mp('dens', 2, 0)

    ansys.Run("/SOLU")
    ansys.Antype(2, "NEW")
    ansys.Run("BCSOPT,,INCORE")
    ansys.Modopt("LANB", 1)
    ansys.Mxpand(elcalc='YES')
    ansys.Run("/GOPR")
    ansys.Solve()
    ansys.Finish()
    ansys.exit()
    """

    filename = os.path.join(testfiles_path, "comp_hex_beam.rst")
    result = pymapdl_reader.read_binary(filename)

    components = ["MY_COMPONENT", "MY_OTHER_COMPONENT"]
    result.plot_nodal_solution(
        0, node_components=components, off_screen=True, sel_type_all=False
    )
    result.plot_nodal_stress(0, "x", node_components=components, off_screen=True)
    result.plot_principal_nodal_stress(
        0, "SEQV", node_components=components, off_screen=True
    )


def test_file_close(tmpdir):
    tmpfile = str(tmpdir.mkdir("tmpdir").join("tmp.rst"))
    shutil.copy(examples.rstfile, tmpfile)
    rst = pymapdl_reader.read_binary(tmpfile)
    nnum, stress = rst.nodal_stress(0)
    del rst
    os.remove(tmpfile)  # tests file has been correctly closed


@skip_plotting
@pytest.mark.skipif(not HAS_FFMPEG, reason="requires imageio_ffmpeg")
def test_animate_nodal_solution(tmpdir, result):
    temp_movie = str(tmpdir.mkdir("tmpdir").join("tmp.mp4"))
    result.animate_nodal_solution(
        0, n_frames=20, movie_filename=temp_movie, loop=False, off_screen=True
    )
    assert np.any(result.grid.points)
    assert os.path.isfile(temp_movie)


def test_loadbeam():
    linkresult_path = os.path.join(testfiles_path, "link1.rst")
    linkresult = pymapdl_reader.read_binary(linkresult_path)
    assert np.any(linkresult.grid.cells)


def test_reaction_forces():
    rst = pymapdl_reader.read_binary(os.path.join(testfiles_path, "vm1.rst"))
    nnum, forces = rst.nodal_static_forces(0)
    assert np.allclose(nnum, [1, 2, 3, 4])
    assert np.allclose(forces[:, 1], [-600, 250, 500, -900])


class TestThermalResult:
    def test_nodal_displacement(self, thermal_rst):
        assert thermal_rst._is_thermal
        assert thermal_rst.result_dof(0) == ["TEMP"]
        with pytest.raises(AttributeError):
            thermal_rst.nodal_displacement()

    def test_nodal_velocity(self, thermal_rst):
        assert thermal_rst._is_thermal
        assert thermal_rst.result_dof(0) == ["TEMP"]
        with pytest.raises(AttributeError):
            thermal_rst.nodal_velocity(0)

    def test_nodal_acceleration(self, thermal_rst):
        assert thermal_rst._is_thermal
        assert thermal_rst.result_dof(0) == ["TEMP"]
        with pytest.raises(AttributeError):
            thermal_rst.nodal_acceleration(0)

    def test_nodal_solution(self, thermal_rst):
        assert thermal_rst._is_thermal
        assert thermal_rst.result_dof(0) == ["TEMP"]
        with pytest.raises(AttributeError):
            thermal_rst.plot_nodal_solution(0, "NORM")

    def test_plot_nodal_solution(self, thermal_rst):
        assert thermal_rst._is_thermal
        assert thermal_rst.result_dof(0) == ["TEMP"]
        with pytest.raises(ValueError):
            thermal_rst.plot_nodal_solution(0, "ROTX")


def test_plot_temperature(thermal_rst):
    cpos = thermal_rst.plot_nodal_temperature(0, return_cpos=True)
    if cpos is not None:
        assert isinstance(cpos, CameraPosition)


def test_file_not_found():
    with pytest.raises(FileNotFoundError):
        pymapdl_reader.read_binary("not_a_file.rst")


def test_file_not_supported():
    with pytest.raises(RuntimeError):
        pymapdl_reader.read_binary(os.path.join(testfiles_path, "file.esav"))


@skip_no_ansys
def test_nodal_time_history(mapdl, transient_thermal):
    mapdl._prioritize_thermal = (
        True  # shouldn't be necessary (needed only on mapdl.local = False)
    )
    rst = mapdl.result
    nnum, data = rst.nodal_time_history()
    assert np.allclose(nnum, mapdl.mesh.nnum)
    for i in range(mapdl.post_processing.nsets):
        mapdl.set(1, i + 1)
        assert np.allclose(data[i].ravel(), mapdl.post_processing.nodal_temperature())
