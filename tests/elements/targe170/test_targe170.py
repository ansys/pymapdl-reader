import os

import numpy as np
import pytest

from ansys.mapdl import reader as pymapdl_reader

path = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope="module")
def rst():
    rst_file = os.path.join(path, "file.rst")
    return pymapdl_reader.read_binary(rst_file)


def test_tshape(rst):
    assert isinstance(rst.mesh.tshape, np.ndarray)
    assert 99 in rst.mesh.tshape
    assert 7 in rst.mesh.tshape
    assert len(rst.mesh.tshape) == 207


def test_tshape_keyopt(rst):
    assert isinstance(rst.mesh.tshape_key, dict)
    assert 99 in rst.mesh.tshape_key.values()
    assert 7 in rst.mesh.tshape_key.values()
    assert len(rst.mesh.tshape_key) == 7
    assert 1 in rst.mesh.tshape_key
    assert 2 in rst.mesh.tshape_key
    assert 7 in rst.mesh.tshape_key

    assert rst.mesh.tshape_key[1] == 0
    assert rst.mesh.tshape_key[2] == 99
    assert rst.mesh.tshape_key[3] == 0
    assert rst.mesh.tshape_key[4] == 7
    assert rst.mesh.tshape_key[7] == 0


def test_et_id(rst):
    assert isinstance(rst.mesh.et_id, np.ndarray)
    assert 1 in rst.mesh.et_id
    assert 2 in rst.mesh.et_id
    assert len(rst.mesh.et_id) == 207
