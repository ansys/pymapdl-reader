import platform
import os
import warnings

import pytest
from pyvista.plotting import system_supports_plotting

from ansys.mapdl.reader import examples


HAS_IMAGEIO = True
try:
    import imageio_ffmpeg
except ImportError:
    HAS_IMAGEIO = False

try:
    shaft = examples.download_shaft_modal()
except:
    warnings.warn('Unable to execute ``examples.download_shaft_modal``')
    shaft = None


IS_MAC = platform.system() == 'Darwin'
skip_plotting = pytest.mark.skipif(not system_supports_plotting() or IS_MAC,
                                   reason="Requires active X Server")
skip_no_shaft = pytest.mark.skipif(shaft is None, reason="Requires example file")


def test_load_verif():
    for filename in examples.vmfiles.values():
        assert os.path.isfile(filename)


@skip_plotting
def test_show_hex_archive():
    examples.show_hex_archive(off_screen=True)


def test_load_result():
    examples.load_result()


@skip_plotting
def test_show_displacement():
    examples.show_displacement(off_screen=True)


@skip_plotting
def test_show_stress():
    examples.show_stress(off_screen=True)


def test_load_km():
    examples.load_km()


@skip_plotting
def test_show_cell_qual():
    examples.show_cell_qual(meshtype='tet', off_screen=True)
    examples.show_cell_qual(meshtype='hex', off_screen=True)


@skip_plotting
@skip_no_shaft
@pytest.mark.skipif(not HAS_IMAGEIO, reason='Requires imageio_ffmpeg')
def test_shaft_animate(tmpdir):
    filename = str(tmpdir.mkdir("tmpdir").join('tmp.mp4'))
    shaft.animate_nodal_solution(5, element_components='SHAFT_MESH',
                                 comp='norm', loop=False,
                                 nangles=10, show_edges=True,
                                 off_screen=True,
                                 movie_filename=filename)

@skip_plotting
@skip_no_shaft
def test_shaft_nodal_solution_ncomp(tmpdir):
    filename = str(tmpdir.mkdir("tmpdir").join('tmp.mp4'))
    shaft.plot_nodal_solution(5, node_components='N_AREA_BC1', sel_type_all=False,
                              off_screen=True)


@skip_plotting
@skip_no_shaft
def test_shaft_plot_screenshot(tmpdir):
    filename = str(tmpdir.mkdir("tmpdir").join('tmp.png'))
    shaft.plot_nodal_solution(0, off_screen=True, screenshot=filename)
    assert os.path.isfile(filename)
