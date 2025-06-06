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
from ansys.mapdl.reader._rst_keys import element_index_table_info
from ansys.mapdl.reader.dis_result import DistributedResult

try:
    from ansys.mapdl.core import _HAS_ANSYS
except:
    _HAS_ANSYS = False

if os.environ.get("SKIP_ANSYS", "").upper() == "TRUE":
    _HAS_ANSYS = False


test_path = os.path.dirname(os.path.abspath(__file__))
testfiles_path = os.path.join(test_path, "testfiles")

IS_MAC = platform.system() == "Darwin"
skip_no_ansys = pytest.mark.skipif(not _HAS_ANSYS, reason="Requires ANSYS installed")
skip_plotting = pytest.mark.skipif(
    not system_supports_plotting() or IS_MAC, reason="Requires active X Server"
)


@pytest.fixture()
def beam_blade():
    filename = os.path.join(testfiles_path, "dist_rst", "blade_stations", "beam3_0.rst")
    return pymapdl_reader.read_binary(filename)


@pytest.fixture(scope="module")
def static_dis():
    filename = os.path.join(testfiles_path, "dist_rst", "static", "file0.rst")
    return pymapdl_reader.read_binary(filename)


@pytest.fixture(scope="module")
def static_rst():
    filename = os.path.join(testfiles_path, "dist_rst", "static", "file.rst")
    return pymapdl_reader.read_binary(filename)


@pytest.fixture(scope="module")
def thermal_solution(mapdl):
    mapdl.clear()

    # create a simple beam and impose some thermal conditions on it
    mapdl.prep7()
    mapdl.mp("kxx", 1, 45)
    mapdl.et(1, 90)
    mapdl.block(-0.3, 0.3, -0.46, 1.34, -0.2, -0.2 + 0.02)
    mapdl.vsweep(1)

    mapdl.finish()
    mapdl.run("/SOLU")
    mapdl.asel("S", vmin=3)
    mapdl.nsla()
    mapdl.d("all", "temp", 5)
    mapdl.asel("S", vmin=4)
    mapdl.nsla()
    mapdl.d("all", "temp", 100)
    mapdl.allsel()

    mapdl.solve()
    mapdl.finish()
    mapdl.post1()
    mapdl.set(1, 1)


def test_not_a_dis_rst(tmpdir):
    filename = os.path.join(testfiles_path, "dist_rst", "static", "file.rst")
    tmp_file = os.path.join(str(tmpdir), "tmp0.rth")
    shutil.copy(filename, tmp_file)
    with pytest.raises(RuntimeError):
        DistributedResult(tmp_file)


@skip_no_ansys
def test_not_all_found(thermal_solution, mapdl, tmpdir):
    if not mapdl._distributed:
        return
    filename = os.path.join(mapdl.directory, "file0.rth")
    tmp_file = os.path.join(mapdl.directory, "tmp0.rth")
    shutil.copy(filename, tmp_file)
    with pytest.raises(FileNotFoundError):
        dist_rst = pymapdl_reader.read_binary(tmp_file)


@skip_no_ansys
def test_temperature(thermal_solution, mapdl, tmpdir):
    if not mapdl._distributed:
        return

    ans_temp = mapdl.post_processing.nodal_temperature
    dist_rst = pymapdl_reader.read_binary(os.path.join(mapdl.directory, "file0.rth"))

    # normal result should match
    rst = mapdl.result  # normally not distributed
    nnum, rst_temp = rst.nodal_temperature(0)
    assert np.allclose(rst_temp, ans_temp)

    # distributed result should match
    dist_nnum, dist_temp = dist_rst.nodal_temperature(0)
    assert np.allclose(dist_temp, ans_temp)
    assert np.allclose(dist_nnum, nnum)


@skip_no_ansys
def test_plot_temperature(thermal_solution, mapdl):
    if not mapdl._distributed:
        return

    dist_rst = pymapdl_reader.read_binary(os.path.join(mapdl.directory, "file0.rth"))
    if cpos is not None:
        cpos = dist_rst.plot_nodal_temperature(0)
    assert isinstance(cpos, CameraPosition)


def test_blade_result(beam_blade):
    filename = os.path.join(
        testfiles_path, "dist_rst", "blade_stations", "ans_xdisp.npz"
    )
    ans_rst = np.load(filename)

    nnum, disp = beam_blade.nodal_displacement(0)
    mask = np.isin(ans_rst["nnum"], nnum)
    assert np.allclose(disp[:, 0], ans_rst["x_disp"][mask])


def test_plot_blade_result(beam_blade):
    cpos = beam_blade.plot_nodal_displacement(0)
    if cpos is not None:
        assert isinstance(cpos, CameraPosition)


def test_nodal_stress(static_dis, static_rst):
    nnum_dis, stress_dis = static_dis.nodal_stress(0)
    nnum, stress = static_rst.nodal_stress(0)
    assert np.allclose(nnum_dis, nnum)
    assert np.allclose(stress_dis, stress, equal_nan=True)


def test_plot_nodal_stress(static_dis, static_rst):
    cpos = static_dis.plot_nodal_stress(0, "x")
    if cpos is not None:
        assert isinstance(cpos, CameraPosition)


def test_nodal_reaction(static_dis, static_rst):
    nnum_dis, data_dis = static_dis.nodal_static_forces(0)
    nnum, data = static_rst.nodal_static_forces(0)
    assert np.allclose(nnum_dis, nnum)
    assert np.allclose(data_dis, data, equal_nan=True)


def test_nodal_principal_stress(static_dis, static_rst):
    nnum_dis, data_dis = static_dis.principal_nodal_stress(0)
    nnum, data = static_rst.principal_nodal_stress(0)
    assert np.allclose(nnum_dis, nnum)
    assert np.allclose(data_dis, data, equal_nan=True)


def test_plot_nodal_principal_stress(static_dis):
    cpos = static_dis.plot_principal_nodal_stress(0, "SEQV")
    if cpos is not None:
        assert isinstance(cpos, CameraPosition)


rtypes = ["EMS", "ENF", "ENS", "ENG", "EEL", "ETH", "EUL", "EMN"]


@pytest.mark.parametrize("rtype", rtypes)
def test_element_solution_data(static_dis, static_rst, rtype):
    elemnum_dis, element_data_dis, enode_dis = static_dis.element_solution_data(
        0, rtype
    )
    elemnum, element_data, enode = static_rst.element_solution_data(0, rtype)

    # assert len(element_data_dis) == len(element_data)
    assert np.allclose(elemnum_dis, elemnum)

    for i in range(static_dis.mesh.n_elem):
        if element_data_dis[i] is None:
            assert element_data[i] is None
        else:
            assert np.allclose(element_data_dis[i], element_data[i], equal_nan=True)

        assert np.allclose(enode_dis[i], enode[i])


def test_element_stress(static_dis, static_rst):
    elemnum_dis, element_stress_dis, enode_dis = static_dis.element_stress(0)
    elemnum, element_stress, enode = static_rst.element_stress(0)

    assert len(element_stress_dis) == len(element_stress)
    assert np.allclose(elemnum_dis, elemnum)

    for i in range(static_dis.mesh.n_elem):
        if element_stress_dis[i] is None:
            assert element_stress[i] is None
        else:
            assert np.allclose(element_stress_dis[i], element_stress[i], equal_nan=True)

        assert np.allclose(enode_dis[i], enode[i])


rtypes = ["ENS", "EPT", "ETH", "EEL", "ENG"]


@pytest.mark.parametrize("rtype", rtypes)
def test_save_as_vtk(tmpdir, static_dis, static_rst, rtype):
    filename = str(tmpdir.mkdir("tmpdir").join("tmp.vtk"))
    static_dis.save_as_vtk(filename, result_types=[rtype])

    grid = pv.read(filename)
    key = "Nodal Solution %d" % 0
    assert key in grid.point_data
    arr = grid.point_data[key]
    assert np.allclose(arr, static_dis.nodal_solution(0)[1], atol=1e-5)

    key = "%s %d" % (element_index_table_info[rtype], 0)
    assert key in grid.point_data
    arr = grid.point_data[key]
    _, rst_arr = static_dis._nodal_result(0, rtype)
    if rst_arr.shape[1] == 1:
        rst_arr = rst_arr.ravel()
    assert np.allclose(arr, rst_arr, atol=1e-5, equal_nan=True)


def test_cylindrical_nodal_stress(static_dis):
    nnum_dis, stress_dis = static_dis.cylindrical_nodal_stress(0)
    nnum, stress = static_dis.cylindrical_nodal_stress(0)
    assert np.allclose(nnum_dis, nnum)
    assert np.allclose(stress_dis, stress, equal_nan=True)


def test_plot_cylindrical_nodal_stress(static_dis):
    cpos = static_dis.plot_cylindrical_nodal_stress(0, "Z")
    if cpos is not None:
        assert isinstance(cpos, CameraPosition)


def test_nodal_thermal_strain(static_dis):
    nnum_dis, data_dis = static_dis.nodal_thermal_strain(0)
    nnum, data = static_dis.nodal_thermal_strain(0)
    assert np.allclose(nnum_dis, nnum)
    assert np.allclose(data_dis, data, equal_nan=True)


def test_plot_nodal_thermal_strain(static_dis):
    cpos = static_dis.plot_nodal_thermal_strain(0, "Z")
    if cpos is not None:
        assert isinstance(cpos, CameraPosition)


def test_nodal_elastic_strain(static_dis):
    nnum_dis, data_dis = static_dis.nodal_elastic_strain(0)
    nnum, data = static_dis.nodal_elastic_strain(0)
    assert np.allclose(nnum_dis, nnum)
    assert np.allclose(data_dis, data, equal_nan=True)


def test_plot_nodal_elastic_strain(static_dis):
    cpos = static_dis.plot_nodal_elastic_strain(0, "Z")
    if cpos is not None:
        assert isinstance(cpos, CameraPosition)


def test_result_does_not_exist(static_dis):
    with pytest.raises(ValueError):
        static_dis.nodal_plastic_strain(0)


@skip_plotting
def test_animate_nodal_solution(static_dis):
    static_dis.animate_nodal_solution(0, loop=False)


# maybe add this...
# nnum_dis, data_dis = static_dis.plot_nodal_static_forces(0)
