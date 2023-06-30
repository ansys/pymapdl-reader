import os

from ansys.mapdl.reader import common

TEST_PATH = os.path.dirname(os.path.abspath(__file__))
TESTFILES_PATH = os.path.join(TEST_PATH, "testfiles")
BEAM44_RST = os.path.join(TESTFILES_PATH, "beam44.rst")


def test_beam44():
    result = common.read_binary(BEAM44_RST)
    assert result.grid.n_cells
    assert result.grid.n_points
    nnum, disp = result.nodal_solution(0)
    assert nnum.size == result.grid.n_points
    assert disp.any()
